# TODO: investigate inverting clustering: do it by species and see if there is any secondary depth effect
#       similarly, ordinate species by depth and see which groups are driven in which direction by which depth

# load libraries quietly
lib <- function(...) suppressPackageStartupMessages(library(...))

lib(here)
lib(EcolUtils)
lib(beyonce)
lib(plotly)
lib(betapart)
lib(patchwork)
lib(dendextend)
lib(rfishbase)
lib(fs)
lib(vegan)
lib(viridis)
lib(ggrepel)

lib(sf)
lib(mapdata)
lib(marmap)
lib(ggspatial)
lib(cowplot)
lib(tidyverse)


# set random seed for reproducibility
set.seed(31337)

# Community setup and initial analysis ------------------------------------

# shorten calls to suppressMessages (readr outputs all sorts of nonsense)
sm <- suppressMessages            

# set the project, directory, and markers for this analysis
active_project <- "mesophotic"
project_dir <- here::here("pipeline_output",active_project)
sample_pattern <- "^BMSP[0-9]+$"  # pattern to match sample IDs
blank_pattern <- "^B[0-9]+$"      # pattern to match extraction blank IDs
abundance_threshold <- 0.1     # minimum threshold for relative abundance
abundance_threshold <- 5     # minimum threshold for relative abundance

remove_families <- c("Hominidae","Bovidae","Felidae","Salmonidae")

# read count transformation, can be "sqrt", "wisconsin", "relative", or anything else (like "none")
read_transform <- "none"
# whether to drop zotus with no assigned taxonomy (NA's across the board)
drop_unknown <- FALSE

# format to save figures
fig_format <- "svg"

max_blank <- 5      # maximum reads in blank to believe it
min_total <- 2e4    # minimum total reads in sample 
# TODO: maybe get rid of min_total

insect_classify <- FALSE
include_unassigned <- TRUE        # include zOTUs that didn't blast to anything

distance_methods <- c("sim","jaccard","bray")
distance_method <- "sim"
plot_theme <- "light"

# save a ggplot figure to preferred format
save_fig <- function(plotz,basename,format="pdf",...) {
  device = switch(
    format,
    pdf = cairo_pdf,
    svg = svg,
    eps = cairo_ps,
    NULL
  )
  filename <- path_ext_set(basename,format)
  if (!is.null(device)) {
    ggsave(filename,plotz,device=device,...)
  }
}

# convert an adonis result to a tibble in a useful way
as_tibble_adonis <- function(a) {
  a %>%
    as_tibble(rownames="term") %>%
    rename(df=Df,ss=SumOfSqs,r2=R2,pseudo_f=F,p_value=`Pr(>F)`)
}

plotz_theme <- function(which="dark") {
  which <- ifelse(which %in% c("dark","light","transparent"),which,"dark")
  thm <- switch(
    which,
    dark = theme_bw() + 
    theme(
      plot.background = element_rect(fill="#1a1a1a"),
      text = element_text(color="grey80"),
      panel.background = element_rect(fill="#000000"),
      legend.background = element_rect(fill="grey27"),
      legend.key = element_rect(fill="grey76", color = "black")
    ),
    light = theme_bw() + 
    theme(
      plot.background = element_rect(fill="white"),
      text = element_text(color="black"),
      panel.background = element_rect(fill="white"),
      legend.background = element_rect(fill="grey95"),
      legend.key = element_rect(fill="grey95", color = "black")
    ),
    transparent = ggfun::theme_transparent() + 
    theme(
      text = element_text(color="grey80"),
      axis.title = element_text(color="grey80"),
      axis.line = element_line(color="grey80"),
      axis.ticks = element_line(color="grey80"),
      axis.text = element_text(color="grey80")
    )
  )
  thm + 
    theme(
      panel.grid = element_blank(),
      legend.key.size = unit(1.4,"line"),
      legend.position = c(0.08,0.94)
    )
    
}

# load our OTU tables and sample data
source("analyze.R")

blank_summary <- communities %>%
  map(~{
    .x$blanks %>%
      mutate(across(
        domain:species,
        ~case_when(
          .x == "dropped" | is.na(.x) ~ "Unidentified",
          TRUE ~ .x
        )
      )) %>%
      arrange(domain,kingdom,phylum,class,order,family,genus,species)
  })

# cut depths into discrete 'deep' and 'shallow' zones
sample_data <- sample_data %>%
  mutate(
    depth=depth_zone,
    depth_zone = cut(depth,c(-Inf,30,Inf),labels=c("Shallow","Deep")),
    depth_zone45 = cut(depth,c(-Inf,45,Inf),labels=c("Shallow","Deep")),
    depth_f=recode_factor(factor(str_c(depth,'m')), "0m" = "Surface")
  ) %>%
  mutate(dive = factor(dive))

# dive computer temperature data
temp_data <- read_csv(here("data","dives.csv"),col_names = c("diveno","date","time","sample_time","sample_depth","sample_temp","sample_pressure","heart"),skip=1) %>%
  # separate(sample_time,into=c("min","sec"),sep=":",convert=T) %>%
  # mutate(sample_time=min*60+sec,diff=c(0,diff(sample_temp)),diveno=factor(diveno)) %>%
  mutate(diveno = factor(diveno)) %>%
  drop_na(sample_temp) %>%
  left_join(sample_data,by=c("diveno" = "dive")) %>%
  filter(sample_depth >= depth_min-5 & sample_depth <= depth_max+4) %>%
  group_by(id) %>%
  summarise(temp=mean(sample_temp))

sample_data <- sample_data %>%
  left_join(temp_data,by="id")

# figure out which samples constitute negative controls 
# so we can use them to prune out taxa that occur in them
# (fortunately for us there aren't any, I don't think)
negative_controls <- sample_data %>%
  filter(substrate %in% c("extraction blank","bleach solution control")) %>%
  pull(id)


# function to construct a phyloseq object from raw data
# including possible taxa filtration
make_ps <- function(otus,taxa,sample_data,minimals=0,negative_controls=character(0),remove_families=character(0),read_transform="",rarefy=FALSE,rarefy_perm=99,filtr) {
  # filter taxa as necessary
  new_tax <- taxa %>%
    filter({{filtr}}) %>%
    filter(!family %in% remove_families)
  keep_otus <- new_tax %>% pull(OTU)
  
  
  # get filtered OTU table
  rarefied <- otus %>%
    filter(OTU %in% keep_otus) %>%
    pivot_longer(matches(sample_pattern),names_to="site",values_to="reads") %>%
    pivot_wider(names_from="OTU",values_from="reads") %>%
    # filter out by minimum read count
    rowwise() %>%
    mutate(total = sum(c_across(-site))) %>% 
    filter(total >= minimals) %>%
    select(-total) %>%
    column_to_rownames("site")
  
  new_tax <- new_tax %>%
    filter(OTU %in% colnames(rarefied))
  
  # rarefy if we want to
  if (rarefy) {
    rarefied <- rarefied %>%
      rrarefy.perm(n=rarefy_perm)
  }
  # create phyloseq object
  ps <- as_ps2(rarefied,new_tax,sample_data)
  
  
  # prune taxa with zero reads or reads occurring in negative controls
  pruno <- filter_taxa(ps,function(x) {
    sum(x[negative_controls],na.rm=T) == 0 & sum(x) > 0
  })
  ps <- prune_taxa(pruno,ps)
  
  
  # get rid of insects, because we know they don't live there
  # we have to include is.na because NA != anything == NA
  ps <- ps %>%
    subset_taxa(is.na(class) | class != "Insecta")
  
  # finally, do transformation if one has been requested
  ps <- ps_standardize(ps,read_transform)
  
  return(ps)
}

minimals <- 300
# construct our phyloseq objects
message("constructing phyloseq objects")
datasets <- c("complete","rarefied") %>%
  set_names() %>%
  map(~{
    rarefy <- .x == "rarefied"
    rarefy_perm <- 99
    
    # full dataset
    msg <- str_glue("making all taxa {.x} dataset") 
    message(msg)
    all_taxa <- communities %>%
      map(
        ~make_ps(
          .x$raw,
          .x$tax_data,
          sample_data,
          negative_controls=negative_controls,
          read_transform=read_transform,
          remove_families=remove_families,
          rarefy=rarefy,
          rarefy_perm=rarefy_perm
        )
      )
    
    # animals subset
    msg <- str_glue("making animals {.x} dataset") 
    message(msg)
    animals <- communities %>%
      keep_names(~.x != "fish") %>% 
      map(
        ~make_ps(
          .x$raw,
          .x$tax_data,
          sample_data,
          minimals=minimals,
          negative_controls=negative_controls,
          remove_families=remove_families,
          read_transform=read_transform,
          filtr=kingdom == "Metazoa",
          rarefy=rarefy,
          rarefy_perm = rarefy_perm
        )
      )
    
    # benthic subset
    msg <- str_glue("making benthic {.x} dataset") 
    message(msg)
    remove_classes <- c("Hexanauplia","Appendicularia","Thaliacea",NA_character_)
    benthic <- communities %>%
      keep_names(~.x != "fish") %>%
      map(
        ~make_ps(
          .x$raw,
          .x$tax_data,
          sample_data,
          minimals=minimals,
          negative_controls=negative_controls,
          remove_families=remove_families,
          read_transform=read_transform,
          filtr=kingdom == "Metazoa" & phylum != "Ctenophora" & !class %in% remove_classes,
          rarefy=rarefy,
          rarefy_perm = rarefy_perm
          # filtr=kingdom == "Metazoa" & phylum != "Ctenophora" & !class %in% remove_classes & !is.na(class)
          # filtr=kingdom == "Metazoa" & !class %in% remove_classes 
        )
      )
    
    # smash datasets together
    list(all=all_taxa,animals=animals,benthic=benthic)
  })


# text map for plots
plot_text <- c(
  'fish' = 'Fishes (16S rRNA)',
  'inverts' = '"Eukaryotes" (18S rRNA)',
  'metazoans' = '"Metazoans" (COI)'
)
plot_text2 <- c(
  'fish' = '16S-fishes',
  'inverts' = '18S-eukaryotes',
  'metazoans' = 'COI-metazoans'
)

# set up depth zone color palette
pal <- rev(beyonce_palette(75))[5:12]
pal <- pal[c(1:3,5,4,6:7)]
pal[4] <- '#827222'
names(pal) <- levels(sample_data$depth_f)

# beta diversity calculations ---------------------------------------------

beta_wrapper <- function(x,method) {
  if (method %in% c("sim","sorensen","jaccard")) {
    x %>%
      ps_standardize("pa") %>%
      otu_table() %>%
      beta.pair(if_else(method == "sim","sorensen","jaccard")) %>%
      return()
  } else if (method == "bray") {
    dist <- x %>%
      distance(method="bray")
    list(beta.bray = dist) %>%
      return()
  } else {
    return(NULL)
  }
}

# pairwise (by depth zone) beta diversity metrics
message("caluclating mean pairwise dissimilarity")
beta_pairs <- datasets %>%
  map(~{
    .x %>%
      map(~{
        dataset <- .x
        distance_methods %>%
          set_names() %>%
          map(~{
            dm <- .x
            dataset %>%
              map(~{
                sd <- sample_tibble(.x)
                .x %>%
                  beta_wrapper(method=dm) %>%
                  imap_dfr(~{
                    .x %>%
                      as("matrix") %>%
                      as_tibble(m,rownames="row") %>%
                      pivot_longer(-row,names_to="col",values_to="dist") %>%
                      left_join(sd %>% select(sample,depth_f),by=c("row" = "sample")) %>%
                      left_join(sd %>% select(sample,depth_f),by=c("col" = "sample"),suffix = c("_s1","_s2")) %>%
                      mutate(row=factor(row),col=factor(col)) %>%
                      filter(as.numeric(depth_f_s1) > as.numeric(depth_f_s2)) %>%
                      select(sample1=row,sample2=col,depth1=depth_f_s1,depth2=depth_f_s2,dist) %>%
                      group_by(depth1,depth2) %>%
                      summarise(sd=sd(dist),dist=mean(dist)) %>%
                      mutate(measurement = .y) %>%
                      select(depth1,depth2,dist,sd,measurement)
                  })
              })
          })
      })
  })
  


# overall beta diversity metrics
message("caluclating overall beta diversity statistics")
beta_diversity <- datasets %>%
  map(~{
    .x %>%
      map(~{
        dataset <- .x
        distance_methods %>%
          set_names() %>%
          map(~{
            dm <- .x
            dataset %>%
              imap(~{
                .x %>%
                  ps_standardize("pa") %>%
                  otu_table() %>%
                  beta.multi(if_else(dm == "sim","sorensen","jaccard"))
              })
          })
      })
  })

# start here --------------------------------------------------------------



# permanova analyses ------------------------------------------------------
anovas <- datasets %>%
  map(~{ # rarefaction
    .x %>%
      map(~{ # dataset
        dataset <- .x
        distance_methods %>%
          set_names() %>%
          map(~{ # distance method
            dm <- .x
            dataset %>%
              imap(~{
                ps <- .x
                min_group <- 2
                # go through permanova terms and make sure we're only using groups with >=3 members
                c("depth_zone","depth_zone45","depth_f","station_grouping") %>%
                  set_names() %>%
                  map(~{
                    perm <- 9999
                    group <- .x
                    pps <- ps_min_group(ps,!!sym(group),2)
                    sd <- sample_tibble(pps)
                    dd <- distance(pps,method=dm,binary=dm != "bray")
                    formula <- as.formula(str_glue("dd ~ {.x}"))
                    adonis2(formula, data=sd, permutations=perm)
                  })
              }) 
          })
      })
  })


# make a nice little table of anova results
# this assumes the name of the list entry is the same as the
# name of the variable being examined
anova_table <- anovas %>%
  imap_dfr(~{ # rarefaction
    .x %>%
      map2_dfr(c("Complete dataset","Metazoans","Benthic metazoans"),~{ # dataset
        .x %>%
          keep_names(~.x != "bray") %>%
          imap_dfr(~{ # distance method
            # name_map <- c("sim" = "Simpson ($\\sim$)", "jaccard" = "Jaccard ($\\jac$)","bray" = "Bray-Curtis")
            name_map <- c("sim" = "Simpson ($\\sim$)", "jaccard" = "Jaccard ($\\jac$)")
            method <- name_map[.y]
            .x %>%
              imap_dfr(~{ # taxa
                marker_name <- plot_text2[.y]
                .x %>%
                  imap_dfr(~{ # analysis
                    # as_tibble(.x,rownames="term") %>%
                    .x %>%
                      as_tibble_adonis() %>%
                      filter(!term %in% c("Residual","Total")) %>%
                      select(term,pseudo_f,p_value) %>%
                      mutate(analysis=.y)
                  }) %>%
                  mutate(marker=marker_name)
              }) %>%
              mutate(index=method)
          }) %>%
          mutate(dataset=.y) %>%
          select(dataset,index,marker,term,everything())   
      }) %>%
      mutate(rarefied = .y == "rarefied")
  })


# manuscript figures ------------------------------------------------------
figure_dir <- dir_create(here("output/figures"),recurse=TRUE)


#### Figure: beta diversity heatmaps (everything)
# map distance types to pretty names
beta_title_map <- list(
  "sim" = c(
    beta.sim = "Turnover",
    beta.sne = "Nestedness",
    beta.sor = "Overall\n(Sørensen)"
  ), 
  "jaccard" = c(
    beta.jtu = "Turnover",
    beta.jne = "Nestedness",
    beta.jac = "Overall\n(Jaccard)"
  )
)
beta_label_map <- list(
  all = list(
    fish = c("A","B","C"),
    inverts = c("D","E","F"),
    metazoans = c("G","H","I")
  ),
  animals = list(
    inverts = c("A","B","C"),
    metazoans = c("D","E","F")
  ),
  benthic = list(
    inverts = c("A","B","C"),
    metazoans = c("D","E","F")
  )
)

beta_map <- list(
  factor_levels = list(
    "jaccard" = c("beta.jac","beta.jtu","beta.jne"),
    "sim" = c("beta.sor","beta.sim","beta.sne")
  ),
  measurement = c(
    "sim" = "beta.sim",
    "jaccard" = "beta.jac",
    "bray" = "beta.bray"
  )
)

method_map <- c(
  'sim' = 'Simpson',
  'jaccard' = 'Jaccard',
  'bray' = 'Bray-Curtis'
)

dataset_map <- c(
  'all' = 'Complete dataset',
  'animals' = 'Metazoans',
  'benthic' = 'Benthic metazoans'
)


# create composite figure of all beta diversity heatmaps
beta_pairs_composite <- beta_pairs %>%
  map(~{
    .x %>%
      imap(~{
        dataset_name <- .y 
        label_map = beta_label_map[[dataset_name]]
        .x %>% imap(~{
          dm <- .y 
          .x %>%
            imap(~{
              marker <- .y
              wrap_elements(
                .x %>%
                  mutate(
                    measurement = factor(measurement,levels=beta_map$factor_levels[[dm]]),
                    depth1 = fct_rev(depth1)
                  ) %>%
                  group_by(measurement) %>%
                  group_map(~{
                    measurement <- as.character(.y$measurement)
                    ggplot(.x) +
                      geom_tile(aes(x=depth2,y=depth1,fill=dist)) +
                      scale_fill_viridis(option="inferno",name=beta_title_map[[dm]][measurement]) +
                      scale_x_discrete(position="top") +
                      theme_bw() +
                      theme(
                        axis.text.x = element_text(angle=25,vjust=1,hjust=0),
                        panel.grid = element_blank()
                      ) +
                      labs(x="Depth zone",y="Depth zone")
                  }) %>%
                  reduce(`+`) +
                  plot_annotation(tag_levels = list(label_map[[marker]])) &
                  theme(plot.tag = element_text(face="bold"))
              ) 
            }) %>%
            reduce(`/`) 
        })
      })
  })

beta_pairs_composite$complete$animals$sim

# save all the different versions to PDF
beta_pairs_composite %>%
  keep_names(~.x != "rarefied") %>%
  iwalk(~{ # rarefaction
    rare <- if_else(.y == "rarefied","_rarefied_","")
    .x %>%
      iwalk(~{ # dataset
        dataset <- .y
        .x %>% 
          keep_names(~.x != "bray") %>%
          iwalk(~{ # distance method
            dm <- .y
            fn <- path(figure_dir,str_glue("mesophotic_beta_pairs_{dataset}_{dm}{rare}"))
            ht = if_else(dataset == "all",9,6)
            save_fig(.x,fn,fig_format,width=12,height=ht,units="in")
          })
      })
  })

#### Figure: cluster plots
cluster_composite <- beta_pairs %>%
  map(~{ # rarefaction
    .x %>%
      imap(~{ # dataset
        .x %>%
          imap(~{ # distance method
            dm <- .y 
            .x %>%
              imap(~{ # marker
                dd <- .x %>%
                  filter(measurement == beta_map$measurement[[dm]]) %>%
                  bind_rows(
                    tibble(
                      depth1=factor(levels(.$depth1),levels=levels(.$depth1)),
                      depth2=factor(levels(.$depth1),levels=levels(.$depth1)),
                      dist = 0,
                      sd = 0,
                      measurement = beta_map$measurement[[dm]]
                    ) 
                  ) %>%
                  arrange(depth1,depth2) %>%
                  pivot_wider(-c(sd,measurement),names_from=depth2,values_from=dist) %>%
                  column_to_rownames("depth1") %>% 
                  as.matrix() %>% 
                  as.dist()
                title <- plot_text2[.y]
                labels(dd) <- str_c(labels(dd)," ")
                expando <- if_else(.y == "fish",0.6,0.3)
                ggplot(as.dendrogram(hclust(dd),hang=0.5)) +
                  theme(
                    plot.caption = element_text(hjust=0.5,size=14)
                  ) + 
                  expand_limits(y=-0.59)
              }) %>%
              reduce(`+`) + 
              plot_annotation(tag_levels="A") &
              theme(plot.tag = element_text(face="bold"))
          })
      })
  })
  
# put all the cluster diagrams in one big figure
# simpson version
big_cluster_sim <- cluster_composite %>%
  map(~{
    .x %>%
      map(~.x$sim) %>%
      reduce(`/`) + 
      plot_annotation(tag_levels="A") &
      theme(plot.tag = element_text(face="bold"))
  })
big_cluster_sim$complete

# save the simpson version of the figure
big_cluster_sim %>%
  iwalk(~{
    rare <- if_else(.y == "rarefied","_rarefied","")
    fn <- path(figure_dir,str_glue("mesophotic_cluster_all_sim{rare}"))
    save_fig(.x,fn,fig_format,width=12,height=9,units="in")
  })

# jaccard version
big_cluster_jac <- cluster_composite %>%
  map(~{
    .x %>%
      map(~.x$jaccard) %>%
      reduce(`/`) + 
      plot_annotation(tag_levels="A") &
      theme(plot.tag = element_text(face="bold"))
  })

# save the jaccard version of the figure
big_cluster_jac %>%
  iwalk(~{
    rare <- if_else(.y == "rarefied","_rarefied","")
    fn <- path(figure_dir,str_glue("mesophotic_cluster_all_jaccard{rare}"))
    save_fig(.x,fn,fig_format,width=12,height=9,units="in")
  })


#### Figure: shallow vs deep ordinations
# generate depth zone ordinations
depth_zones <- c("depth_zone","depth_zone45")
ord_zone_composite <- datasets %>%
  map(~{ # rarefaction?
    .x %>%
      imap(~{ # datasets
        dataset <- .x
        distance_methods %>%
          set_names() %>%
          map(~{ # distance methods
            dm <- .x
            depth_zones %>%
              set_names() %>%
              map(~{ # depth zones
                zone <- .x
                dataset %>%
                  imap(~{ # marker
                    title <- plot_text2[.y]
                    dpal <- pal[c(1,length(pal))]
                    names(dpal) <- c("Shallow","Deep")
                    p <- plot_betadisp(.x, group=zone, method=dm, list=TRUE,binary=dm != "bray")
                    p$plot <- p$plot +
                      scale_fill_manual(values=dpal,name="Depth Zone", drop=FALSE) +
                      plotz_theme("light") +
                      xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
                      ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
                      theme(legend.position = "right")
                    p$plot
                  }) %>%
                  reduce(`+`) +
                  plot_layout(guides="collect") +
                  plot_annotation(tag_levels = "A") &
                  theme(plot.tag = element_text(face="bold"))
              })
          })
      })
  })

# Save composite figures of ordinations by depth zone
# split by 30/45 and distance method
# return the plots so we can see them
zone_plotz <- depth_zones %>%
  set_names() %>%
  map(~{
    dz <- .x
    distance_methods %>%
      set_names() %>%
      map(~{
        dm <- .x
        zone <- if_else(dz == "depth_zone","shallowdeep30","shallowdeep45")
        plotz <- ord_zone_composite$complete$all[[dm]][[dz]] /
          ord_zone_composite$complete$animals[[dm]][[dz]] /
          ord_zone_composite$complete$benthic[[dm]][[dz]] +
          plot_annotation(tag_levels = "A") + 
          plot_layout(guides="collect") &
          theme(plot.tag = element_text(face="bold"))
        fn <- path(figure_dir,str_glue("mesophotic_ord_{zone}_{dm}"))
        save_fig(plotz,fn,fig_format,width=12,height=12,units="in")
        return(plotz)
      })
  })


#### Figure: depth zone ordinations
# generate composite figures of ordinations by individual depth zones
ord_composite <- datasets %>%
  map(~{ # rarefaction
    .x %>%
      map(~{ # dataset
        dataset <- .x
        distance_methods %>%
          set_names() %>%
          map(~{ # distance method
            dm <- .x
            dataset %>%
              imap(~{ # marker
                title <- plot_text2[.y]
                to_plot <- .x
                 
                p <- plot_betadisp(to_plot, group="depth_f", method=dm, list=TRUE, usable_groups=TRUE,binary=dm != "bray")
                p$plot <- p$plot +
                  scale_fill_manual(values=pal,name="Depth",drop=FALSE) +
                  plotz_theme("light") +
                  xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
                  ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
                  theme(legend.position = "right") +
                  guides(color="none")
                p$plot
              }) %>%
              reduce(`+`) +
              plot_layout(guides="collect") +
              plot_annotation(tag_levels = "A") &
              theme(plot.tag = element_text(face="bold"))
          })
      })
  })

# Save composite figures of ordinations by depth zone
# split by distance method
# return the plots so we can see them
ord_plotz <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    plotz <- ord_composite$complete$all[[dm]] /
      ord_composite$complete$animals[[dm]] /
      ord_composite$complete$benthic[[dm]] +
      plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(face="bold"))
    fn <- path(figure_dir,str_glue("mesophotic_ord_{dm}"))
    save_fig(plotz,fn,fig_format,width=12,height=12,units="in")
    return(plotz)
  })

#### Figure: ordinations by site
# generate composite ordinations figure
ord_sites_composite <- datasets %>%
  map(~{
    .x %>%
      map(~{
        dataset <- .x
        distance_methods %>%
          set_names() %>%
          map(~{
            dm <- .x
            dataset %>%
              imap(~{
                title <- plot_text2[.y]
                spal <- pal[c(1,4,7)]
                names(spal) <- .x %>%
                  sample_tibble() %>%
                  pull(station_grouping) %>%
                  unique()
                p <- plot_betadisp(.x, group="station_grouping", method=dm, list=TRUE,binary=dm != "bray")
                p$plot <- p$plot +
                  scale_fill_manual(values=spal,name="Site", drop=FALSE) +
                  plotz_theme("light") + 
                  xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
                  ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
                  theme(legend.position = "right")
                p$plot
              }) %>%
              reduce(`+`) +
              plot_layout(guides="collect") +
              plot_annotation(tag_levels = "A") &
              theme(plot.tag = element_text(face="bold"))
          })
      })
  })

# save the simpson version, un-rarefied
fn <- path(figure_dir,"mesophotic_ord_sites_all_sim")
save_fig(ord_sites_composite$complete$all$sim,fn,fig_format,width=12,height=4,units="in")

# save the jaccard version, un-rarefied
fn <- path(figure_dir,"mesophotic_ord_sites_all_jaccard")
save_fig(ord_sites_composite$complete$all$jaccard,fn,fig_format,width=12,height=4,units="in")

#### Figure: fish depth distributions vs eDNA detections

# get composite fish data / taxonomy grid
# use the raw data so we don't miss anything
fish <- communities$fish$raw %>%
  inner_join(communities$fish$tax_data,by="OTU") %>%
  select(domain:species,OTU,matches(sample_pattern)) %>%
  filter(species != "dropped" & !is.na(species))

# pull species names
fish_species <- fish %>%
  distinct(species) %>%
  pull(species)

# join fish detection depths with fishbase depth ranges
detected_depth <- fish %>%
  filter(species %in% fish_species) %>%
  select(family,species,matches(sample_pattern)) %>%
  pivot_longer(matches(sample_pattern),names_to="sample",values_to="reads") %>%
  inner_join(sample_data %>% select(id,depth),by=c("sample" = "id")) %>%
  left_join(
    fb_tbl("species") %>%
      mutate(sci_name = paste(Genus, Species)) %>%
      select(sci_name, contains('Depth')) ,
    by = c("species" = "sci_name")
  ) %>%
  select(family,species,depth,reads,fb_shallow=DepthRangeShallow,fb_deep=DepthRangeDeep) %>%
  filter(reads > 0) %>%
  mutate(deeper = depth > fb_deep)

# let's do the little statistical tests that Muff et al 2022 did

## t-test to compare number of reads inside and outside of published depth zone
## this was to test for "abundant center" effect. we want the reads to be the same
## inside and outside of their zone, so we're not just picking up a stray, is the idea

# get the data for the t-test
deeper <- detected_depth %>%
  filter(deeper == TRUE) %>%
  pull(reads)
shallower <- detected_depth %>%
  filter(deeper == FALSE) %>%
  pull(reads)

# do the t-test as a welch two sample kind
t.test(deeper,shallower,paired=FALSE,alternative="two.sided",var.equal = FALSE)

## here we do a kendall rank correlation test
## testing the relationship between read count and depth detection difference
## according to Muff & pals, a significant relationship might indicate a niche shift

# assemble data for the correlation test
dd <- detected_depth %>%
  filter(depth > fb_deep | depth < fb_shallow) %>%
  mutate(
    depth_diff = case_when(
      depth < fb_shallow ~ abs(depth-fb_shallow),
      TRUE ~ abs(depth-fb_deep)
    )
  )

# kendall rank correlation test
cor.test(dd$depth_diff,dd$reads,method="kendall",exact = F)

# finish making the depth detection summary
# taking only detections with more than 10 reads
# keeping only species with deepest depth <150m
depth_data <- detected_depth %>%
  filter(reads >= 10) %>%
  group_by(family,species,depth,fb_deep,fb_shallow) %>%
  summarise(n=n(),reads=sum(reads)) %>%
  ungroup() %>%
  filter(fb_deep < 150 | fb_shallow > 90)  %>%
  ungroup() %>%
  group_by(species) %>%
  mutate(depth_diff = depth - max(fb_deep)) %>%
  ungroup() %>%
  mutate(
    depth_range = fb_deep-fb_shallow,
    species = fct_reorder(species,depth_range,.desc = TRUE)
  ) 
  
# what is the most number of observations for a species?
# we need this to plot the depth ranges elegantly
# since it'll plot the line range as many time as there
# are observations, we just duplicate the observations
# for whoever doesn't have the max number
most <- depth_data %>%
  count(species) %>% 
  pull(n) %>% 
  max()

# this is a little hack to make the depth ranges look good
# we basically make them be drawn an equal number of times per
# species, so the lines are all equally thick.
depth_data <- depth_data %>%
  group_by(species) %>%
  group_modify(~{
    grp <- .x
    to_add <- most-nrow(grp)
    if (to_add > 0) {
      # here we duplicate the row *to_add* times but blank out the reads 
      add <- map_dfr(seq(to_add),~grp %>% slice(1) %>% mutate(reads=NA))
      grp <- grp %>%
        bind_rows(add)
    }
    return(grp)
  }) 

# make the fish depth detections plot
depth_plotz <- ggplot(depth_data) + 
  # this "error bar" is the reported depth range
  geom_errorbar(aes(x=species,ymin=fb_shallow,ymax=fb_deep),width=0.5,color="dodgerblue4") +
  # these points are where we've actually detected them
  geom_point(aes(x=species,y=depth,size=reads,color=recode_factor(factor(str_c(depth,'m')), "0m" = "Surface")),alpha=0.7) +
  # put a black bar at the mean detection depth
  stat_summary(aes(x=species,y=depth),geom="point",fun="mean",col="black",size=10,shape="-") + 
  # color them by depth zone
  scale_color_manual(values=pal,name="Depth zone") +
  # size them by read count
  scale_size(name="Read count",labels=scales::comma,range=c(2,12),breaks=c(1e2,1e3,1e4,5e4,1e5)) +
  # flip it upside down because it's depth
  scale_y_reverse(breaks=seq(0,300,by=30),minor_breaks=seq(0,300,by=15)) +
  # put the x axis on top
  scale_x_discrete(position="top") + 
  plotz_theme("light") +
  theme(
    # rotate the species names so they're legible
    axis.text.x = element_text(face="italic",angle=20,vjust=1,hjust=0),
    legend.position = "right",
    legend.key = element_rect(color=NA),
    panel.grid.major = element_line(color="grey95"),
    panel.grid.minor = element_line(color="grey97"),
    panel.grid.major.x = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size=7))) +
  ylab("Detection depth (m)") + 
  xlab("Species") 

depth_plotz
fn <- path(figure_dir,"mesophotic_fish_depth")
save_fig(depth_plotz,fn,fig_format,width=12,height=10,units="in")

#### Figure: species accumulations
# test species accumulation against various models
all_models <- datasets$complete$all %>%
  imap(~{
    # get the sample data
    sd <- sample_tibble(.x)
    specs <- otu_tibble(.x) %>%
      # join depth category to otu table
      inner_join(sd %>% select(sample,depth_f),by="sample") %>%
      select(sample,depth_f,everything()) %>%
      # smash the otu table data into a nested sub-table
      nest(otu_table=-depth_f) %>%
      # mutate the nested
      mutate(otu_table = map(otu_table,~{
        # here we use specaccum to permute a bunch of species accumulations
        # for this subset of the data (depth zone)
        # .x is the otu table for this subset
        # we'll likely get wraning messages, so suppress them
        sm(
          curve <- .x %>%
            column_to_rownames("sample") %>% 
            as("matrix") %>%
            # convert the data to presence/absence first
            decostand("pa",margin=NULL) %>%               
            # do the accumulation with 1000 permutations
            specaccum(method="random",permutations=1000)
        ) 
        # make the permutation results into a data frame
        # and flip it to long format
        as_tibble(curve$perm) %>%
          mutate(sites=row_number()) %>%
          select(sites,everything()) %>%
          pivot_longer(-sites,names_to = "permutation", values_to = "richness") %>%
          select(sites,richness)
      })) %>%
      unnest(otu_table)
    
    # test each model against the permuted accumulation data
    modls <- list(
      asymp = nls(richness ~ SSasymp(sites, Asym, R0, lrc),  data=specs),
      gompertz = nls(richness ~ SSgompertz(sites, Asym,  xmid, scal), data=specs),
      `michaelis-menten` = nls(richness ~  SSmicmen(sites, Vm, K), data=specs),
      logis = nls(richness ~ SSlogis(sites,  Asym, xmid, scal), data=specs) 
    )
    # get the best AIC value
    AICs <- modls %>%
      map_dbl(~{
        if (!is_null(.x))
          AIC(.x)
        else
          Inf
      }) %>%
      sort()
    
    # best model's name
    best_mod <- names(AICs)[1]
    # the best model itself
    pred <- modls[[best_mod]]
    # generate some data with the model
    predictions <- tibble(
      x = seq(1:2000),
      y = predict(pred,list(sites=seq(1:2000)))
    )
    
    # find where the slope goes below 1
    slopes <- diff(predictions$y)
    slopes <- c(1,slopes)
    # this is our predicted "asymptote", but really it's just
    # where fewer than one species gets added per additional replicate
    asymp <- predictions$x[slopes < 1][1]
    # predicated number of taxa at asymptote
    asymp_taxa <- predictions$y[slopes < 1][1]
    # how many would we expect with 6 replicates?
    taxa_6 <- predictions$y[6]
    # what's the proportion of total at "asymptote"?
    percent_taxa <- taxa_6/asymp_taxa
    
    # provide some output
    cat("Community:",.y,"\n")
    cat("Best model:",best_mod,"\n")
    cat("Asymptote at:",asymp,"replicates\n")
    cat("Taxa at asymptote:",asymp_taxa,"\n")
    walk(c(3,4,6),~{
      tt <- predictions$y[.x]
      pc <- tt/asymp_taxa
      cat(str_glue("Taxa at {.x} replicates: {tt}\n\n"))
      cat(str_glue("Percent of taxa represented by {.x} replicates: {scales::percent(pc)}\n\n"))
    })
    cat("\n\n")
    # return all the data as a list
    list(models = modls, best_model = best_mod, predictions = predictions, asymptote=asymp, asymptote_taxa=asymp_taxa, our_taxa = taxa_6, accum_data=specs)
  })

# make composite plot from the accumulation data
accum_composite <- all_models %>%
  imap(~{
    d <- .x$accum_data %>%
      group_by(depth_f,sites) %>%
      summarise(sd=sd(richness),richness=mean(richness))
    p <- .x$predictions
    tit <- plot_text2[.y]
    ggplot() + 
      # plot lines with error bars rather than every permuted point
      # these are the accumulation curves for each depth zone
      geom_line(data=d,aes(x=sites,y=richness,color=depth_f)) + 
      geom_errorbar(data=d,aes(x=sites,ymin=richness-sd,ymax=richness+sd),width=0.5) + 
      # this is the best-fit model line
      geom_line(data=p,aes(x=x,y=y),color="blue",size=1.1) + 
      # this marks the "ideal" replication level
      geom_vline(xintercept = .x$asymptote,color="firebrick",size=0.7,linetype="dashed") +
      scale_color_manual(values=pal,name="Depth Zone") + 
      # keep us within a reasonable range
      xlim(1,25) +
      xlab("Replicates") + 
      ylab("zOTUs") + 
      theme_bw() + 
      theme(panel.grid = element_blank())
  }) %>%
  # smash the plots together
  reduce(`+`) + 
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face="bold"))

accum_composite
# ggsave(path(figure_dir,"mesophotic_accum.pdf"),accum_composite,device=cairo_pdf,width=12,height=4,units="in")
fn <- path(figure_dir,"mesophotic_accum")
save_fig(accum_composite,fn,fig_format,width=12,height=4,units="in")

#### Figure: sample rarefaction curves
rarecurves <- datasets$complete$all %>%
  map(~{
    curvedata <- .x %>%
      otu_table()  %>%
      as("matrix") %>% 
      rarecurve(step=500,tidy=TRUE) 
    minreads <- curvedata %>%
      group_by(Site) %>%
      summarise(maxreads=max(Sample)) %>%
      pull(maxreads) %>%
      min()
    
    curvedata %>%
      ggplot() + 
      geom_line(aes(x=Sample,y=Species,color=Site)) + 
      scale_color_viridis(option="magma",discrete=TRUE,drop=FALSE,name="Sample",guide="none") + 
      geom_vline(xintercept = minreads,color="firebrick",size=0.7,linetype="dashed") +
      scale_x_continuous(labels = scales::label_comma())  + 
      plotz_theme("light") + 
      # theme(
      #   legend.background = element_rect(fill="white"),
      #   legend.key = element_rect( color = NA)
      # ) + 
      labs(x = "Sequencing depth",y="zOTUs")
  }) 

curve_composite <- rarecurves %>% 
  reduce(`/`) + 
  plot_layout(guides="collect") + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face="bold")) 
curve_composite

# ggsave(path(figure_dir,"mesophotic_rarefactions.pdf"),curve_composite,device=cairo_pdf,width=6,height=6,units="in")
fn <- path(figure_dir,"mesophotic_rarefactions")
save_fig(curve_composite,fn,fig_format,width=6,height=6,units="in")

#### Figure: site map

# make an sf object for what we want to display
spots <- sample_data %>% 
  filter(
    project == "mesophotic",
    !is.na(lon) & !is.na(lat)
  ) %>%
  distinct(station_grouping,lat,lon,.keep_all = T) %>% 
  # let's cheat the point over a bit to make it easier to see
  mutate(
    lat = case_when(
      str_detect(station,"Outhouse") ~  19.636009,
      TRUE ~ lat
    ),
    lon = case_when(
      str_detect(station,"Outhouse") ~  -156.005123,
      TRUE ~ lon
    ),
    # get the spelling right
    station_grouping = str_replace(station_grouping,"Honaunau","Hōnaunau")
  ) %>%
  st_as_sf(coords=c("lon","lat"),crs=4326,remove=F)

# make up a bounding box
box <- c(ymax=20.337941429661853, xmin=-156.62139587608522, ymin=18.84548284767367, xmax=-155.4117057257134)

# set up some geometries
pbm <- cbind(c(box['xmin'],box['xmax'],box['xmax'],box['xmin'],box['xmin']),
             c(box['ymin'],box['ymin'],box['ymax'],box['ymax'],box['ymin']))
boxpoly <- st_polygon(list(pbm))
boxsfc <- st_sfc(boxpoly,crs=4326)
plotbox <- st_sf(name="study area",geometry=boxsfc)

# this is the hawaii box
hbox <- c("ymin"=18.63314470051002, "xmax"=-154.6469259428272, "ymax"=22.434801224940056, "xmin"=-160.76622913125863)

# get the world map and filter it to just hawaii
hawaii <- st_as_sf(maps::map("worldHires",plot=F,fill=T),crs=4326) %>% 
  filter(ID == "Hawaii")

# a little cheat value
x <- 0.01

# get bathymetry data
bath <- getNOAA.bathy(lon1=hbox['xmin'],lon2=hbox['xmax'],lat1=hbox['ymin'],lat2=hbox['ymax'],resolution = 1)

# convert it into something we can use
bm <- bath %>%
  as("matrix") %>%
  as_tibble(rownames = "lon") %>%
  mutate(lon = as.numeric(lon)) %>%
  pivot_longer(-lon,names_to = "lat",values_to = "value") %>%
  mutate(across(everything(),as.numeric)) %>%
  filter(value < 0)

# plot our main map
mainmap <- ggplot() + 
  # show depth contours
  geom_contour(data=bm,mapping=aes(x = lon, y = lat, z = value), breaks=c(-10,-20,seq(-50,-10000,by=-50)), color="#8FCFEA") + 
  # show hawaii island
  geom_sf(data=hawaii, fill="#8BCA7B",color="#4f4f4f",size=0.5) +
  # show sampling stations
  geom_sf(data=spots, shape=25,size=6,color="black",fill="firebrick") +
  coord_sf(xlim = c(box['xmin']-x,box['xmax']+x),
           ylim = c(box['ymin']-x,box['ymax']+x)) +
  # show sampling station names
  # do this twice so we can have a transparent
  # background but opaque foreground
  geom_label_repel(
    data=spots,
    aes(label=station_grouping,geometry=geometry),
    point.padding = 12,
    stat = "sf_coordinates",
    min.segment.length = 0,
    alpha=0.4,
    seed = 31337
  ) +
  geom_label_repel(
    data=spots,
    aes(label=station_grouping,geometry=geometry),
    point.padding = 12,
    stat = "sf_coordinates",
    min.segment.length = 0,
    alpha=1,
    fill = NA,
    seed = 31337
  ) +
  # make it look nice
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  theme(
    axis.text = element_text(color="black", size=10),
    axis.title = element_blank(),
    panel.background = element_rect(fill="#D1EDFB") ,
    plot.background = element_rect(fill="white",color=NA),
    text=element_text(color="black")      
  ) +
  # add a scale bar
  annotation_scale(location="br", width_hint=0.5) 

# make the inset plot of the hawaiian islands
inset <- ggplot() + 
  # plot hawaii
  geom_sf(data=hawaii, fill="#8BCA7B",color="#4f4f4f",size=0.1) +
  # plot the box to show where the big map is
  geom_sf(data=plotbox, color="red", fill="transparent") + 
  coord_sf(xlim = c(hbox['xmin']-x,hbox['xmax']+x),
           ylim = c(hbox['ymin']-x,hbox['ymax']+x)) +
  theme_void() +
  theme(
    panel.background = element_rect(fill="#d1edfb", color="black") 
  ) +
  # add a scale bar
  annotation_scale(location="bl",width_hint=0.5)


# smash them together
sitemap <- ggdraw(mainmap) + 
  draw_plot(
    inset,
    x = -0.130,
    y = 0.35,
    scale=0.4
  ) 
# and save it
save_fig(sitemap,path(figure_dir,"mesophotic_map"),"pdf",width=8.5,height=11,units="in")


# manuscript tables -------------------------------------------------------
table_dir <- dir_create(here("output/tables"),recurse=TRUE)

write_ms_table <- function(tbl,file,caption="",bold_header=TRUE,...) {
  if (bold_header) {
    tbl <- tbl %>% rename_with(~str_c("**",.x,"**"))
  }
  if (caption != "") {
    write_lines(str_c(": ",caption),file)
    write_csv(tbl,file,append=TRUE,col_names = TRUE,...)
  } else {
    write_csv(tbl,file,...)    
  }
}

#### Table (supplemental?): PERMANOVA results
term_map = c(
  "depth_zone" = "Depth zone\n(shallow = 0--30m)",
  "depth_zone45" = "Depth zone\n(shallow = 0--45m)",
  "depth_f" = "15m depth interval",
  "station_grouping" = "Sites"
)

anova_table %>%
  mutate(pseudo_f=round(pseudo_f,3), p_value=round(p_value,4)) %>%
  mutate( term = term_map[term] ) %>%
  mutate(
    p_value = case_when(
      p_value < 0.05 ~ str_glue("**{p_value}**"),
      TRUE ~ as.character(p_value)
    ),
    rarefied = if_else(rarefied,"Rarefied to minimum read depth","Un-rarefied")
  ) %>%
  select(-analysis) %>%
  rename(Dataset=dataset,`Dissimilarity Index`=index,`Assay`=marker,`Term`=term,`Pseudo-F`=pseudo_f,`p-value`=p_value,Rarefaction=rarefied) %>%
  select(Dataset,Rarefaction,everything()) %>%
  write_ms_table(path(table_dir,"mesophotic_anovas.csv"),"Results of PERMANOVA analyses",bold_header = TRUE)

#### Supplemental table: eDNA reads summary
all_samples <- sample_data %>%
  select(`Sample ID` = id,substrate,station_grouping,depth_f)

# normalized data
normed <- datasets$complete$all %>%
  imap_dfr(~{
    otu_tibble(.x) %>%
      column_to_rownames("sample") %>%
      rowSums() %>%
      enframe(name = "Sample ID", value = "Normalized sequence reads") %>%
      mutate(Assay=plot_text[.y])
  })

# build summary table of eDNA reads across samples, etc.
reads_summary <- 
  communities %>%
  imap_dfr(~{
    gene <- markers %>%
      filter(description == .y) %>%
      pull(gene)
    
    counts <- .x$unfiltered %>%
      select(where(is.numeric),-matches("Blast")) %>%
      colSums() %>%
      enframe(name = "Sample ID", value = "Sequence reads")
    
    s <- all_samples %>%
      left_join(counts,by="Sample ID")
    
    rs <- s %>% 
      mutate(
        "Site" = station_grouping,
        "Depth" = depth_f,
        "Assay" = plot_text[.y],
        "Amplification" = !is.na(`Sequence reads`),
        "Type" = case_when(
          str_detect(substrate,"control") ~ "field control",
          str_detect(`Sample ID`,sample_pattern) ~ "eDNA Sample",
          str_detect(`Sample ID`,blank_pattern) ~ "extraction blank",
          TRUE ~ ""
        ),
        "Passing" = !is.na(`Sequence reads`) & `Sequence reads` > min_total
      ) %>%
      select(`Sample ID`,Assay,Site,Depth,`Sequence reads`,Type,Passing,Amplification) %>%
      arrange(`Sample ID`)
  }) %>%
  left_join(normed,by=c('Sample ID','Assay'))

rs <- reads_summary %>%
  replace_na(list(`Sequence reads`=0,`Normalized sequence reads`=0)) %>%
  mutate(
    Passing = case_when(
      Passing & str_detect(Assay,"16S") ~ "16S",
      Passing & str_detect(Assay,"18S") ~ "18S",
      Passing & str_detect(Assay,"COI") ~ "COI",
      TRUE ~ NA_character_
    ),
    Amplification = case_when(
      Amplification & str_detect(Assay,"16S") ~ "16S",
      Amplification & str_detect(Assay,"18S") ~ "18S",
      Amplification & str_detect(Assay,"COI") ~ "COI",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(`Sample ID`,Site,Depth) %>%
  summarise(
    `Sequence reads (16S; 18S; COI)`=str_c(scales::comma(na.omit(`Sequence reads`),accuracy = 1),collapse="; "),
    `Normalized sequence reads (16S; 18S; COI)`=str_c(scales::comma(na.omit(`Normalized sequence reads`),accuracy = 1),collapse="; "),
    `Sample type`=unique(Type),
    `Successful amplification`=str_c(na.omit(Amplification),collapse="; "),
    `Passing QA/QC`=str_c(na.omit(Passing),collapse="; ")
  ) %>%
  mutate(
    Site = factor(Site),
    Site = fct_relevel(Site,"none",after=Inf),
    order = case_when(
      str_detect(`Sample ID`,sample_pattern) ~ 0,
      str_detect(`Sample ID`,blank_pattern) ~ 1,
      TRUE ~ 999
    )
  ) %>%
  arrange(order,`Sample ID`) %>%
  select(-order,-starts_with("Normalized"))

# write reads summary to manuscript tables directory
# write_ms_table(rs,path(table_dir,"mesophotic_reads_summary.csv"),"Sequencing results for mesophotic eDNA samples from West Hawai&#x02BB;i across assay types",na="",bold_header = TRUE)
write_ms_table(rs,path(table_dir,"mesophotic_reads_summary.csv"),"Sequence read summary for West Hawai&#x02BB;i eDNA samples",na="",bold_header = TRUE)

#### Supplemental table: mean pairwise beta diversity data
bd <- beta_pairs %>%
  imap_dfr(~{ # rarefaction
    .x %>% imap_dfr(~{
      .x %>% 
        keep_names(~.x != "bray") %>%
        imap_dfr(~{
          .x %>%
            imap_dfr(~{
              .x %>%
                group_by(measurement) %>%
                summarise(mean=round(mean(dist),2),sd=round(sd(dist),2)) %>%
                ungroup() %>%
                mutate(marker=.y)
            }) %>%
            mutate(method=.y)
        }) %>%
        mutate(dataset=.y)
    }) %>%
    mutate(
      measurement = case_when(
        method == "sim" ~ beta_title_map$sim[measurement],
        method == "jaccard" ~ beta_title_map$jaccard[measurement]
      ),
      method = method_map[method],
      marker = plot_text2[marker]
    ) %>%
    arrange(dataset,desc(method),marker,match(measurement,c("Overall (Jaccard)","Overall (Sørensen)","Turnover","Nestedness"))) %>%
    mutate(
      dataset = dataset_map[dataset],
      rarefied = if_else(.y == "rarefied","Rarefied to minimum read depth","Un-rarefied")
    ) %>%
    select(Dataset=dataset, Rarefaction=rarefied, `Dissimilarity Index`=method, Assay=marker,  
           Measurement=measurement, Mean=mean, `Std. Dev.`=sd)
  })
write_ms_table(bd,path(table_dir,"mesophotic_beta_diversity.csv"),
               caption="Beta diversity summary statistics",
               bold_header = TRUE)

reps <- datasets$complete$all %>%
  imap_dfr(~{
    dataset <- .y
    modl <- all_models[[dataset]]
    zotus <- psmelt(.x) %>% 
      filter(Abundance > 0) %>%
      group_by(depth_f) %>%
      summarise(zotus = n_distinct(OTU))
    rep <- .x %>%
      sample_tibble() %>%
      count(depth_f,name = "passed") %>%
      inner_join(zotus,by="depth_f") %>%
      mutate(dataset=dataset) %>%
      mutate(
        percent = scales::percent(zotus/modl$asymptote_taxa),
        zotus = scales::comma(zotus,accuracy = 1),
        status = str_glue("{passed} ({zotus} zOTUs, ~{percent})")
      ) %>%
      select(depth_f,dataset,status)
  }) %>%
  mutate(dataset = plot_text2[dataset]) %>%
  pivot_wider(names_from="dataset",values_from = "status") %>%
  arrange(depth_f) %>%
  rename(`Depth Zone`=depth_f)
write_ms_table(reps,path(table_dir,"mesophotic_sample_reps.csv"),caption="Sample replicates passing QA/QC by depth zone with recovered zOTUs and percentage of total predicted diversity",bold_header=T)

# manuscript include file -------------------------------------------------

# write PERMANOVA and beta diversity results to an include file for the
# manuscript so I don't have to keep changing the numbers when they change here
resource_dir <- dir_create(here("output/resources"),recurse=TRUE)
include_file <- path(resource_dir,"mesophotic_include.m4")


file_delete(include_file)

# write permanova results
anovas$complete %>%
  iwalk(~{ # dataset
    dataset <- .y
    .x %>%
      # keep_names(~.x != "bray") %>%
      iwalk(~{ # distance method
        method <- .y
        .x %>%
          iwalk(~{ # marker
            marker <- .y
            lines <- .x %>%
              imap_dfr(~{
                as_tibble(.x,rownames="term") %>%
                  filter(term == .y) %>%
                  select(term,pseudo_f=`F`,p_value=`Pr(>F)`) %>%
                  mutate(
                    pseudo_f = str_glue("pseudo-F = {round(pseudo_f,2)}"),
                    p_value = case_when(
                      p_value < 0.001 ~ "*p* \\< 0.001",
                      p_value < 0.01 ~ "*p* \\< 0.01",
                      TRUE ~ as.character(str_glue("*p* = {round(p_value,2)}"))
                    )
                  )
              }) 
            pseudo_f <- str_glue("define({{{{{dataset}_{method}_{marker}_{lines$term}_f}}}},{{{{{lines$pseudo_f}}}}})")
            p_val <- str_glue("define({{{{{dataset}_{method}_{marker}_{lines$term}_p}}}},{{{{{lines$p_value}}}}})")
            lines <- c(pseudo_f,p_val)
            write_lines(lines,include_file,append=TRUE)
          }) 
      }) 
  })

stat_map <- c(
  "beta.sor" = "overall",
  "beta.sim" = "turnover",
  "beta.sne" = "nestedness",
  "beta.jac" = "overall",
  "beta.jtu" = "turnover",
  "beta.jne" = "nestedness",
  "beta.bray" = "overall"
)

# write beta diversity stats
beta_diversity$complete %>%
  iwalk(~{
    dataset_name <- .y
    .x %>%
      keep_names(~.x != "bray") %>%
      iwalk(~{
        method <- .y
        .x %>% iwalk(~{
          marker <- .y
          .x %>% iwalk(~{
            stat <- stat_map[str_to_lower(.y)]
            .x <- str_pad(round(.x,2),4,side="right",pad="0")
            line <- str_glue("define({{{{{dataset_name}_{method}_{marker}_{stat}}}}},{{{{{.x}}}}})") 
            write_lines(line,include_file,append=TRUE)
          })
        })
      })
  })

# write mean/sd beta diversity stats for depth zone comparisons 
with(
  beta_pairs$complete %>%
    imap_dfr(~{
      .x %>% 
        keep_names(~.x != "bray") %>%
        imap_dfr(~{
          .x %>%
            imap_dfr(~{
              .x %>%
                group_by(measurement) %>%
                summarise(mean=round(mean(dist),2),sd=round(sd(dist),2)) %>%
                ungroup() %>%
                mutate(marker=.y)
            }) %>%
            mutate(method=.y)
        }) %>%
        mutate(dataset=.y)
    }),
  write_lines(
    c(
      str_glue("define({{{{{dataset}_{method}_{marker}_{stat_map[measurement]}_mean}}}},{{{{{mean}}}}})"),
      str_glue("define({{{{{dataset}_{method}_{marker}_{stat_map[measurement]}_sd}}}},{{{{{sd}}}}})")
    ),
    include_file,
    append=TRUE
  )
)

category_map = c(
  "sample" = "sample",
  "station" = "site",
  "depth_f" = "depth"
)

#### eDNA read counts for different categories squished together
datasets$complete %>%
  imap(~{
    dataset <- .y
    .x %>% imap(~{
      ps <- .x
      marker <- .y 
      c("sample","station","depth_f") %>%
        map(~{
          if (.x != "sample") {
            ps <- merge_samples(ps,.x)
          } 
          
          cat <- category_map[.x]
          ss <- sample_sums(ps)
          m <- scales::number(mean(ss),big.mark = ",")
          s <- scales::number(sd(ss),big.mark = ",")
          ms <- str_glue("define({{{{{dataset}_{marker}_reads_{cat}_mean}}}},{{{{{m}}}}})")
          ss <- str_glue("define({{{{{dataset}_{marker}_reads_{cat}_sd}}}},{{{{{s}}}}})")
          c(ms,ss)
        })
    })
  }) %>% 
  unlist() %>%
  write_lines(include_file,append=TRUE)
