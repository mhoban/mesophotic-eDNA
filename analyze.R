# load libraries quietly

library(phyloseq)
library(insect)
library(httr)

################################################################################################
# functions to convert OTU tables & sample data to phyloseq objects, as well as do other things
################################################################################################
vdist <- function(comm, method="jaccard", pa=TRUE,...) {
  if (method == "beta.sim" | method == "sim") {
    comm <- if(pa) decostand(comm,"pa") else comm
    beta.pair(comm,...)$beta.sim
  } else {
    vegdist(comm,method=method,...)
  }
}

# make a matrix out of a data frame, with optional row names
make_matrix <- function(df,row_names=NULL) {
  m <- as.matrix(df)
  if (!is.null(row_names)) {
    rownames(m) <- row_names
  }
  return(m) 
}

# assemble an OTU table and sample data into a phyloseq object
as_ps <- function(otus,sample_pattern,sample_data) {
  # create matrices for OTUs and taxonomy
  otutable <- make_matrix(otus %>% select(matches(sample_pattern)), otus$OTU)
  taxtable <- make_matrix(otus %>% select(domain:species), otus$OTU)
  
  # get sample data
  sd <- sample_data %>% 
    filter(id %in% colnames(otutable)) %>%
    column_to_rownames('id')
  
  #construct phyloseq object
  phyloseq(
    otu_table(otutable, taxa_are_rows = TRUE),
    tax_table(taxtable),
    sample_data(sd)
  )
}

as_ps2 <- function(otutable,taxtable,sample_data) {
  sd <- sample_data %>% 
    filter(id %in% rownames(otutable)) %>%
    column_to_rownames('id')
  # taxtable <- make_matrix(taxtable %>% select(domain:species), taxtable$OTU)
  taxtable <- make_matrix(taxtable %>% select(-OTU), taxtable$OTU)
  phyloseq(
    otu_table(otutable,taxa_are_rows=FALSE),
    tax_table(taxtable),
    sample_data(sd)
  )
}

summarize_communities <- function(comm,ctype='relative') {
  cat(ctype,"\n")
  comm %>% walk2(names(.),~{
    cat(.y)
    sites <- .x[[ctype]] %>%
      pivot_longer(matches(sample_pattern),names_to="site",values_to="read_count") %>%
      distinct(site) %>%
      pull(site)

        sd <- sample_data %>% filter(id %in% sites)
    print(table(sd$station_grouping,sd$depth))
  })
}

is_blank <- function(x) {
  is.na(x) | x == ""
}

map_pal <- function(cg,ci,rev=FALSE) {
  map2_chr(ifelse(!rev,cg,rev(cg)),ci,~brewer.pal(9,.x)[9-.y])
}

sample_tibble <- function(ps, sid = 'sample') {
  as_tibble(as(sample_data(ps),'data.frame'),rownames = sid)
}

taxa_tibble <- function(ps, otus = 'zotu') {
  as_tibble(as(tax_table(ps),'matrix'),rownames=otus)
}

otu_tibble <- function(ps, rownames = 'sample') {
  as_tibble(as(otu_table(ps),'matrix'),rownames=rownames)
}

plot_betadisp <- function(ps, group, method="jaccard", list=FALSE, usable_groups=TRUE,expand=FALSE,binary=TRUE,...) {
  # dd <- vdist(otu_table(ps),method=method)
  if (usable_groups & expand) {
    stop("expand and usable_groups are mutually exclusive")
  }
  
  dd <- distance(ps,method=method,binary=binary)
  sd <- sample_tibble(ps,sid="sample")
  bds <- betadisper(dd,group = sd %>% pull(all_of(group)),...)
  sco <- scores(bds)
  
  centroids <- sco$centroids %>%
    as_tibble(rownames=group) %>%
    rename(x=PCoA1,y=PCoA2)
  sd <- sd %>%
    inner_join(as_tibble(sco$sites,rownames="sample"), by="sample") %>%
    rename(x=PCoA1,y=PCoA2)
  
  # x_var <- bds$eig[1]/sum(bds$eig)
  # y_var <- bds$eig[2]/sum(bds$eig)
  varx <- map(bds$eig,~.x/sum(bds$eig))
  
  hull <- sd
  if (expand) {
    hull <- hull %>%
      count(!!sym(group)) %>%
      filter(n == 2) %>%
      pull(all_of(group)) %>% 
      map_dfr(~{
        first <- hull %>%
          filter(!!sym(group) == .x) %>%
          # slice(1) %>%
          mutate(
            x = jitter(x,factor=0.02),
            y = jitter(y,amount=0.02)
          )
      }) %>%
      bind_rows(hull)
    change <- hull %>%
      group_by(sample) %>%
      filter(n() > 1) %>%
      summarise(x=mean(x),y=mean(y)) 
    sd <- sd %>%
      full_join(change,by="sample") %>%
      mutate(
        x = case_when(
          !is.na(x.y) ~ x.y,
          TRUE ~ x.x
        ),
        y = case_when(
          !is.na(y.y) ~ y.y,
          TRUE ~ y.x
        )
      ) %>% select(-x.x,-x.y,-y.x,-y.y)
  }
  
  # make these things factors so when/if they're filtered out,
  # we still see them in the plot legend (if we so desire)
  hull <- hull %>%
    mutate(across(all_of(group),as.factor))
  centroids <- centroids %>%
    mutate(across(all_of(group),as.factor))
  
  if (usable_groups) {
    to_remove <- hull %>%
      count(!!sym(group)) %>%
      filter(n <= 2) %>%
      pull(!!sym(group))
    
    hull <- hull %>%
      group_by(!!sym(group)) %>%
      filter(n() > 2) %>%
      ungroup()
    centroids <- centroids %>%
      filter(!(!!sym(group) %in% to_remove))
  }
  
  hull <- hull %>%
    group_by(!!sym(group)) %>%
    slice(chull(x,y))

  p <- ggplot(sd) + 
    geom_polygon(data=hull, aes_string(x=quote(x),y=quote(y),fill=group),alpha=0.7) + 
    geom_label_repel(data=centroids, aes_string(x=quote(x), y=quote(y), fill=group, label=group), show.legend = FALSE, segment.color=NA)
  if (list) {
    list(plot = p, betadisp= bds, x_var = varx[[1]], y_var = varx[[2]])
  } else {
    p
  }
}

tukey <- function(bds) {
  if (inherits(bds,'betadisper')) {
    tk <- TukeyHSD(bds)
    as_tibble(tk$group, rownames='comparison') %>%
      separate(comparison,into=c("left","right"),sep="-") %>%
      rename(lower=lwr,upper=upr,p_adj=`p adj`)
  } else {
    return(NULL) 
  }
}

print_anova <- function(model,name="") {
  cat(name,"\n")
  if (inherits(model,"adonis")) {
    tab <- as_tibble(model$aov.tab,rownames = "term") %>%
      rename(pval=`Pr(>F)`) %>% mutate(
        ppval = case_when(
          pval < 0.0001 ~ "*p* < 0.0001",
          pval < 0.001 ~ "*p* < 0.001",
          pval < 0.01 ~ "*p* < 0.01",
          TRUE ~ as.character(str_glue("*p* = {round(pval,2)}"))
        )
      )
    pwalk(tab,~{
      row <- list(...)
      if (!is.na(row$F.Model) & !is.na(row$pval)) {
        p <- row$pval
        s <- str_glue("pseudo-F = {round(row$F.Model,2)}, {row$ppval}")
        ss <- str_glue("{row$term}:\t\t{s}")
        f <- readline(ss)
        clipr::write_clip(s) 
      }
    })
  }
}

pairwise_adonis <- function(comm, factors, permutations = 1000, correction = "fdr", method = "bray") {
  factor_combos <- combn(levels(factors), 2)
  model_output <- map_dfr(array_branch(factor_combos,2), ~{
    fact <- factors[factors %in% .x]
    if (inherits(comm,'dist')) {
      dd <- as.matrix(comm)[factors %in% .x, factors %in% .x]
    } else {
      comm <- as(comm,'matrix')
      dd <- vdist(comm[factors %in% .x,], method = method)
    }
    
    model <- adonis(dd ~ fact, permutations = permutations)
    df <- model$aov.tab$Df[1]
    ss <- model$aov.tab$SumsOfSqs[1]
    f_mod <- model$aov.tab$F.Model[1]
    p_val <- model$aov.tab$`Pr(>F)`[1]
    
    list(
      "left" = .x[1],
      "right" = .x[2],
      "df" = df,
      "ss" = ss,
      "f_model" = f_mod,
      "p_val" = p_val
    )
  })
  
  model_output %>%
    mutate(
      p_adj = p.adjust(p_val, method = correction)
    )
}

gnu_sed <- function() {
	system("sed --version > /dev/null 2>&1",ignore.stderr = TRUE, ignore.stdout = TRUE, intern = FALSE) == 0  
}

get_ncbi_taxonomy <- function() {
  nlf <- here("data","lineage.txt")
  if (!file_exists(nlf)) {
    url <- "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip"
    taxdir <- dir_create(here("data","taxonomy"),recurse=TRUE)
    taxdump <- path(taxdir,"taxdump.zip")
    resp <- GET(url,write_disk(taxdump,overwrite = T),progress())
    if (resp$status_code == 200) {
      unzip(taxdump,exdir=taxdir)       
      lf <- path(taxdir,"rankedlineage.dmp")
      sed_cmd <- "'s!\\t\\|\\t!\\t!g;s!\\t\\|$!!g'"
      if (!gnu_sed()) {
        sed_cmd <- str_c("$",sed_cmd)
      }   
      system2("sed",args=c("-E",sed_cmd,lf),stdout=nlf)
    } else {
      return(NULL)
    }
  }
  read_delim(nlf,delim="\t",col_names = c("taxid","name","species","genus","family","order","class","phylum","kingdom","domain"),col_types=cols()) %>%
    select(taxid,name,domain,kingdom,phylum,class,order,family,genus,species) %>%
    return()
}

# helper function to do purrr-style keep call with list element names
keep_names <- function(.x,.names,...) {
  if (is.character(.names)) {
    idx <- intersect(names(.x),.names)
  } else if (is.function(.names) || is_formula(.names)) {
    .names <- rlang::as_function(.names)
    idx <- .names(names(.x))
    if (is.logical(idx)) {
      idx[is.na(idx)] <- FALSE
    } else if (is.character(idx)) {
      idx <- intersect(names(.x),idx)
    } else if (!is.integer(idx)) {
      abort("If `.names` is a function, it must return an logical, integer, or character vector")
    }
  }
  .x[idx]
}

# filter a phyloseq object so that there are >min representatives per group
ps_min_group <- function(ps,group,min) {
  filter_groups <- ps %>%
    sample_tibble() %>%
    group_by(!!enquo(group)) %>%
    filter(n() >= min) %>%
    pull(sample)
  prune_samples(filter_groups,ps)
}

# standardize phyloseq object
ps_standardize <- function(ps, method="total", ...) {
  vegan_methods <- c("total", "max", "frequency", "normalize", "range", "rank", "standardize", "pa", "chi.square", "hellinger", "log")
  other_methods <- c("wisconsin","sqrt")
  
  trows <- phyloseq::taxa_are_rows(ps)
  if(trows == TRUE) { 
    marg <- 2 
  } else {
    marg <- 1 
  }
  
  otus <- ps %>%
    otu_table() %>%
    as("matrix")
  
  if (method == "relative") {
    method <- "total"
  }
  
  if (method %in% vegan_methods) {
    otus_std <- decostand(otus,method=method,MARGIN=marg,...)
    if (method == "chi.square") {
      otus_std <- t(otus_std)
    }
  } else if (method %in% other_methods) {
    if (method == "wisconsin") {
      if (trows) {
        otus_std <- wisconsin(t(otus))
        otus_std <- t(otus_std)
      } else {
        otus_std <- wisconsin(otus)
      }
    } else if (method == "sqrt") {
      otus_std <- sqrt(otus)
    }
  } else {
    otus_std <- NULL
  }
  
  if (!is.null(otus_std)) {
    otu_table(ps) <- otu_table(otus_std,taxa_are_rows = trows)
  }
  return(ps) 
}


# start here --------------------------------------------------------------
# shorten calls to suppressMessages (readr outputs all sorts of nonsense)
sm <- suppressMessages            
# read sample metadata (this is project-universal)
sample_data <- read_csv(here("data","sample_data.csv"),col_types=cols()) %>%
  filter(project %in% active_project)

# get marker info
markers <- read_csv(here::here("data","marker_config.csv"),col_types=cols())

# make a little lookup table of marker descriptions
marker_descriptions <- markers %>%
  filter(skip == "N") %>%
  pull(gene) %>%
  set_names %>%
  map_chr(~markers %>% filter(gene == .x) %>% pull(description))

# gene targets
targets <- names(marker_descriptions)

# load ranked lineage taxonomy from NCBI
message("loading NCBI taxonomy")
ncbi <- get_ncbi_taxonomy()

communities <- targets %>%
  set_names(marker_descriptions[.]) %>%
  map(~{
    # read the OTU table and do some pre-filtration
    # we assume that rows are taxa and columns are samples
    # use a regex to determine which columns are samples, blanks or everything else
    # first, we get rid of taxa that show up in extraction blanks
    # next, we get rid of samples that have zero reads across all taxa
    # note that you need to use '&&' inside of 'where', rather than '&'
    tab_file <- path(project_dir,str_glue("{.x}_collapsed_taxa.tab"))
    if (file_exists(tab_file)) {
      
      message("loading OTU table and taxonomy data")
      otu_table <- read_tsv(tab_file,col_types=cols()) #%>%
     
      unassigned_file <- path(project_dir,str_glue("{.x}_all_taxa.tab"))
      if (include_unassigned & file_exists(unassigned_file)) {
        unassigned <- read_tsv(unassigned_file,col_types=cols()) %>% 
          filter(!(OTU %in% otu_table$OTU)) 
        classified_file <- path(project_dir,str_glue("{.x}_classified.csv"))    
        # if (insect_classify) {
        #   if (!file_exists(classified_file)) {
        #     classifier_file <- here("classifiers",str_glue("{.x}_classifier.rds"))
        #     zotu_file <- path(project_dir,str_glue("{.x}_zotus.fasta"))
        #     if (all(file_exists(c(classifier_file,zotu_file)))) {
        #       classifier <- readRDS(classifier_file)
        #       zotu <- readFASTA(zotu_file)
        #       cores <- parallel::detectCores() 
        #       cat("running insect classifier model....\n")
        #       classified <- classify(zotu, classifier, threshold = 0.5, metadata = TRUE, offset = -2, mincount = 1, cores=cores)
        #       write_csv(classified,classified_file) 
        #     }
        #   }
        # }
        
        if (file_exists(classified_file)) {
          message("loading insect classified taxonomy")
          insect_taxa <- read_csv(classified_file,col_types=cols()) %>%
            rename(taxid=taxID) %>%
            left_join(ncbi %>% select(taxid,domain),by="taxid") %>%
            mutate(
              domain = case_when(
                taxon == "Eukaryota" ~ "Eukaryota",
                TRUE ~ domain
              ) 
            ) %>%
            select(domain,kingdom:species,OTU=representative,id_method)
          
          unassigned <- unassigned %>%
            left_join(insect_taxa,by="OTU") %>%
            mutate(numberOfUnq_BlastHits=0) %>%
            select(domain:species,OTU,id_method,numberOfUnq_BlastHits,everything())
           
        } else {
          unassigned <- unassigned %>%
            mutate(domain="Unassigned",kingdom="Unassigned",phylum="Unassigned",class="Unassigned",order="Unassigned",family="Unassigned",genus="Unassigned",species="Unassigned",numberOfUnq_BlastHits=0) %>%
            select(domain:species,numberOfUnq_BlastHits,everything()) %>%
            select(any_of(names(otu_table)))
        }
        
        otu_table <- bind_rows(otu_table,unassigned %>% select(all_of(names(otu_table)))) 
      }
      
      otu_table <- otu_table %>%
          rowwise() %>%
          mutate(
            blanks=sum(c_across(matches(blank_pattern))),
          ) %>%
          ungroup() 
      
      message("filtering OTUs found in extraction blanks")
      filtered_blanks <- otu_table %>%
        filter(blanks >= max_blank)
      unfiltered <- otu_table %>% select(-blanks)
      otu_table <- otu_table %>% 
        filter(blanks < max_blank) %>%
        select(-blanks,-matches(blank_pattern)) %>%
        select(!matches(sample_pattern) | where(~(is.numeric(.x) && sum(.x,na.rm=T) >= min_total)))
      
      totals <- otu_table %>% 
        select(matches(sample_pattern)) %>%
        colSums()
      
      # separate otu table and taxonomy table
      tax_data <- otu_table %>%
        select(domain:species,OTU,id_method) %>%
        mutate(across(domain:species,~na_if(.x,"dropped")))
      
      
      curated_file <- path(project_dir,str_glue("{.x}_curated.csv"))
      if (file_exists(curated_file)) {
        message("loading manually curated taxonomy")
        curated <- read_csv(curated_file,col_types=cols())
        c2 <- curated %>%
          left_join(ncbi,by=c("new_name" = "name")) %>%
          pmap_dfr(~{
            row <- list(...)
            row[[row$level]] <- row$new_name
            if (is_blank(row$qualifier)) {
              return(row)
            } else {
              if (row[[row$qualifier_level]] == row$qualifier) {
                return(row)
              } else {
                return(NULL)
              }
            }
          }) %>%
          select(domain:species,OTU) %>%
          mutate(id_method="curated")
        tax_data <- tax_data %>%
          rows_upsert(c2,by="OTU")
      }
      
      # filter out entries where no taxonomy is assigned
      unknown_filtered <- 0
      if (drop_unknown) {
        message("dropping unidentified taxa")
        unid <- tax_data  %>%
          filter(across(domain:species,is.na)) %>% 
          pull(OTU)
        
        unknown_filtered <- length(unid)
        
        otu_table <- otu_table %>%
          filter(!OTU %in% unid)
        tax_data <- tax_data %>%
          filter(!OTU %in% unid)
      }
      
      otu_table <- otu_table %>%
        select(OTU,matches(sample_pattern))
      
      message("filtering by minimum reads and dumping OTUs with zero reads")
      otu_rel <- otu_table %>%
        mutate(
          # across(matches(sample_pattern),~.x/sum(.x)),
          across(
            matches(sample_pattern),
            ~case_when(
              .x < abundance_threshold ~ 0,
              TRUE ~ .x
            )
          )
        ) %>%
        rowwise() %>%
        mutate( asv_reads = sum(c_across(matches(sample_pattern))) ) %>%
        ungroup() %>%
        filter( asv_reads > 0 ) %>% 
        select(-asv_reads) #%>%
        # mutate(
        #   across(
        #     matches(sample_pattern),
        #     ~ .x * totals[cur_column()]
        #   )
        # ) 
      list('raw'=otu_table,'relative'=otu_rel,'tax_data'=tax_data, 'blanks'=filtered_blanks,'unfiltered'=unfiltered,'unknown_filtered'=unknown_filtered)
    } else {
      list('raw'=NULL,'relative'=NULL)
    }
  }) 


