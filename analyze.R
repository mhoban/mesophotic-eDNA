library(tidyverse)
library(here)
library(fs)
library(phyloseq)

################################################################################################
# functions to convert OTU tables & sample data to phyloseq objects, as well as do other things
################################################################################################
vdist <- function(comm, method="jaccard", pa=TRUE,...) {
  if (method == "beta.sim") {
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
  taxtable <- make_matrix(taxtable %>% select(domain:species), taxtable$OTU)
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

plot_betadisp <- function(ps, group, method="jaccard", list=FALSE) {
  dd <- vdist(otu_table(ps),method=method)
  sd <- sample_tibble(ps,sid="sample")
  bds <- betadisper(dd,group = sd %>% pull(all_of(group)))
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
  
  hull <- sd %>%
    group_by(across(all_of(group))) %>%
    slice(chull(x,y))
  p <- ggplot(sd) + 
    geom_polygon(data=hull, aes_string(x=quote(x),y=quote(y),fill=group),alpha=0.7) + 
    geom_label(data=centroids, aes_string(x=quote(x), y=quote(y), fill=group, label=group), show.legend = FALSE)
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


# shorten calls to suppressMessages (readr outputs all sorts of nonsense)
sm <- suppressMessages            
# read sample metadata (this is project-universal)
sample_data <- sm(read_csv(here("data","sample_data.csv")))
markers <- sm(read_csv(here::here("data","marker_config.csv")))
targets <- markers %>%
  filter(skip == "N") %>%
  pull(gene)

communities <- targets %>%
  map(~{
    # read the OTU table and do some pre-filtration
    # we assume that rows are taxa and columns are samples
    # use a regex to determine which columns are samples, blanks or everything else
    # first, we get rid of taxa that show up in extraction blanks
    # next, we get rid of samples that have zero reads across all taxa
    # note that you need to use '&&' inside of 'where', rather than '&'
    tab_file <- path(project_dir,str_glue("{.x}_collapsed_taxa.tab"))
    if (file_exists(tab_file)) {
      
      otu_table <- sm(read_tsv(tab_file)) %>%
        mutate(
          across(domain:species,~replace_na(.x,"unspecified"))
        )
      
      unassigned_file <- path(project_dir,str_glue("{.x}_all_taxa.tab"))
      if (include_unassigned & file_exists(unassigned_file)) {
        unassigned <- read_tsv(unassigned_file) %>% 
          filter(!(OTU %in% otu_table$OTU)) %>%
          mutate(domain="Unassigned",kingdom="Unassigned",phylum="Unassigned",class="Unassigned",order="Unassigned",family="Unassigned",genus="Unassigned",species="Unassigned",numberOfUnq_BlastHits=0) %>%
          select(domain:species,numberOfUnq_BlastHits,everything()) %>%
          select(any_of(names(otu_table)))
        
        otu_table <- bind_rows(otu_table,unassigned) 
      }
      otu_table <- otu_table %>%
          rowwise() %>%
          mutate(
            blanks=sum(c_across(matches(blank_pattern))),
          ) %>%
          ungroup() 
      
      filtered_blanks <- otu_table %>%
        filter(blanks >= max_blank)
      otu_table <- otu_table %>% 
        filter(blanks < max_blank) %>%
        select(-blanks,-matches(blank_pattern)) %>%
        select(!matches(sample_pattern) | where(~(is.numeric(.x) && sum(.x,na.rm=T) > min_total)))
      
      totals <- otu_table %>% 
        select(matches(sample_pattern)) %>%
        colSums()
      tax_data <- otu_table %>%
        select(-matches(sample_pattern),-numberOfUnq_BlastHits) %>%
        pmap_dfr(~{
          items <- list(...)
          dropped <- which(items == "dropped" | is.na(items))
          if (length(dropped) > 0) {
            first_dropped <- min(dropped)
            if (first_dropped > 1) {
              dropped_name <- str_c("Unidentified ",items[first_dropped-1])
            } else {
              dropped_name <- "Unidentified"
            }
            items[dropped] <- dropped_name
          }
          items
        }) %>%
        mutate(
          spp = str_detect(species,"Unidentified") & !str_detect(genus,"Unidentified"),
          species = case_when(
            spp ~ str_c(genus,"spp.",sep = " "),
            TRUE ~ species
          )  
        ) %>%
        select(-spp)
      otu_table <- otu_table %>%
        select(OTU,matches(sample_pattern))
      # convert read count to relative read abundance
      message("WE CHANGED THE RELATIVE THING TO ABSOLUTE, REMEMBER THAT","\n")
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
      list('raw'=otu_table,'relative'=otu_rel,'tax_data'=tax_data, 'blanks'=filtered_blanks)
    } else {
      list('raw'=NULL,'relative'=NULL)
    }
  }) %>%
  set_names(str_c("m",targets))

sample_data <- sample_data %>% 
  filter(str_detect(id,sample_pattern))

