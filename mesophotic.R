library(here)
library(EcolUtils)
library(beyonce)
library(plotly)
library(indicspecies)
library(betapart)
library(patchwork)
library(dendextend)
library(rfishbase)
library(fs)
library(vegan)
library(viridis)
library(ggVennDiagram)


# set random seed for reproduceability
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

# whether to rarefy samples to minimum read depth
rarefy <- FALSE
# whether to transform read counts to eDNA index (wisconsin double standardized)
wisco <- FALSE
# whether to drop zotus with no assigned taxonomy (NA's across the board)
drop_unknown <- TRUE

rarefy_perm <- 99

summary_grouping <- c("station_grouping","station")  # how to group samples for replication report

max_blank <- 5                    # maximum reads in blank to believe it
min_total <- 2e4                   # minimum total reads in sample

insect_classify <- FALSE
include_unassigned <- TRUE        # include zOTUs that didn't blast to anything
filter_unknown <- FALSE

distance_methods <- c("sim","jaccard")
distance_method <- "sim"
plot_theme <- "light"

plotz_theme <- function(which="dark") {
  which <- ifelse(which %in% c("dark","light"),which,"dark")
  thm <- switch(
    which,
    dark = theme_bw() + 
    theme(
      plot.background = element_rect(fill="#1a1a1a"),
      text = element_text(color="grey80"),
      panel.background = element_rect(fill="#000000"),
      legend.background = element_rect(fill="grey27"),
      legend.key = element_rect(fill="grey76", color = "black"),
    ),
    light = theme_bw() + 
    theme(
      plot.background = element_rect(fill="white"),
      text = element_text(color="black"),
      panel.background = element_rect(fill="white"),
      legend.background = element_rect(fill="grey95"),
      legend.key = element_rect(fill="grey95", color = "black"),
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
  )

# figure out which samples constitute negative controls 
# so we can use them to prune out taxa that occur in them
# (fortunately for us there aren't any, I don't think)
negative_controls <- sample_data %>%
  filter(substrate %in% c("extraction blank","bleach solution control")) %>%
  pull(id)

# construct our phyloseq object
comm_ps <- communities %>% 
  map(~{
    rarefied <- .x$raw %>%
      pivot_longer(matches(sample_pattern),names_to="site",values_to="reads") %>%
      pivot_wider(names_from="OTU",values_from="reads") %>%
      column_to_rownames("site") 
    if (rarefy) {
      rarefied <- rarefied %>%
        rrarefy.perm(n=rarefy_perm)
    }
    if (wisco) {
      rarefied <- wisconsin(rarefied)
    }
    ps <- as_ps2(rarefied,.x$tax_data,sample_data)
    pruno <- filter_taxa(ps,function(x) {
      sum(x[negative_controls],na.rm=T) == 0 & sum(x) > 0
    })
    ps <- prune_taxa(pruno,ps)
    ps <- subset_taxa(ps, !family %in% c("Hominidae","Bovidae","Felidae","Salmonidae"))
    
    prune_samples(sample_sums(ps) > 0, ps)
  }) 

# create the metazoan subset of our three communities
minimals <- 1000 # minimum reads to retain a sample
# if we wisconsin'd the data, we have to reconstruct the animals from
# the raw data and re-wisconsin it
if (wisco) {
  animals <- communities %>%
    map(~{
      new_tax <- .x$tax_data %>%
        filter(kingdom == "Metazoa")
      keep_otus <- new_tax$OTU
      rarefied <- .x$raw %>%
        filter(OTU %in% keep_otus) %>%
        pivot_longer(matches(sample_pattern),names_to="site",values_to="reads") %>%
        pivot_wider(names_from="OTU",values_from="reads") %>%
        column_to_rownames("site") 
      if (rarefy) {
        rarefied <- rarefied %>%
          rrarefy.perm(n=rarefy_perm)
      }
      if (wisco) {
        rarefied <- wisconsin(rarefied)
      }
      ps <- as_ps2(rarefied,new_tax,sample_data)
      prune_samples(sample_sums(ps) > 0, ps)
    })
} else { # otherwise just subset
  animals <- comm_ps %>%
    map(~{
      f <- subset_taxa(.x,kingdom == "Metazoa")         # filter by metazoans
      f <- prune_samples(sample_sums(f) >= minimals,f)  # filter by reads >= minimum
      new_otu <- otu_table(f) %>%                       
        as("matrix") 
      if (rarefy) {
        # re-rarefy new otu table to equal depth
        new_otu <- new_otu %>%
          rrarefy.perm(n=rarefy_perm)
        otu_table(f) <- otu_table(new_otu,taxa_are_rows = FALSE) # reassign new rarefied otu table
      }
      return(f)
    })
}

animals <- animals[c("inverts","metazoans")]

# text map for plots
plot_text <- c(
  'fish' = 'Fishes (16S rRNA)',
  'inverts' = '"Eukaryotes" (18S rRNA)',
  'metazoans' = '"Metazoans" (COI)'
)
plot_text2 <- c(
  'fish' = '16S-Fish',
  'inverts' = '18S-Eukaryote',
  'metazoans' = 'COI-Leray'
)

# set up depth zone color palette
pal <- rev(beyonce_palette(75))[5:12]
pal <- pal[c(1:3,5,4,6:7)]
pal[4] <- '#827222'

# beta diversity calculations ---------------------------------------------

# pairwise (by depth zone) beta diversity metrics
beta_pairs <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    comm_ps %>%
      map(~{
        sd <- sample_tibble(.x)
        .x %>%
          otu_table() %>%
          as("matrix") %>%
          decostand("pa") %>%
          beta.pair(if_else(dm == "sim","sorensen","jaccard")) %>%
          map2_dfr(names(.),~{
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


# pairwise (by depth zone) beta diversity metrics (animals only)
beta_pairs_animals <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    animals %>%
      map(~{
        sd <- sample_tibble(.x)
        .x %>%
          otu_table() %>%
          as("matrix") %>%
          decostand("pa") %>%
          beta.pair(if_else(dm == "sim","sorensen","jaccard")) %>%
          map2_dfr(names(.),~{
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

# overall beta diversity metrics
beta_diversity <- comm_ps %>%
  map(~{
    comm <- .x %>%
      otu_table() %>%
      as("matrix") %>%
      decostand("pa")
    comm_merged <- .x %>%
      merge_samples("depth_f") %>%
      otu_table() %>%
      as("matrix") %>%
      decostand("pa")
    bp <- betapart.core(comm)
    bp_merged <- betapart.core(comm_merged)
    sorensen <- comm %>%
      betapart.core() %>%
      beta.multi("sorensen")
    sorensen_merged <- comm %>%
      betapart.core() %>%
      beta.multi("sorensen")
    return(list(unmerged=sorensen,merged=sorensen_merged))
  }) 

beta_diversity_animals <- animals %>%
  map(~{
    .x <- merge_samples(.x,"depth_f")
    comm <- as(otu_table(.x),"matrix")
    comm <- decostand(comm,"pa")
    # now the presence-absence one
    bp <- betapart.core(comm)
    
    
    sorensen <- beta.multi(bp,"sorensen")
    jaccard <- beta.multi(bp,"jaccard")
    
    list(sorensen=sorensen,jaccard=jaccard)
  }) 

# start here --------------------------------------------------------------



# permanova analyses ------------------------------------------------------
anovas <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    comm_ps %>%
      map2(names(.),~{
        perm <- 9999
        sd <- sample_tibble(.x)
        dd <- distance(.x,method=dm)
        list(
          depth_zone = adonis(dd ~ depth_zone, data=sd, permutations=perm),
          depth_zone45 = adonis(dd ~ depth_zone45, data=sd, permutations=perm),
          depth_f = adonis(dd ~ depth_f, data=sd, permutations = perm), 
          station_grouping = adonis(dd ~ station_grouping, data=sd, permutations = perm)
        )
      }) 
  })

animal_anovas <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    animals %>%
      map2(names(.),~{
        perm <- 9999
        sd <- sample_tibble(.x)
        dd <- distance(.x,method=dm)
        list(
          depth_zone = adonis(dd ~ depth_zone, data=sd, permutations=perm),
          depth_zone45 = adonis(dd ~ depth_zone45, data=sd, permutations=perm),
          depth_f = adonis(dd ~ depth_f, data=sd, permutations = perm), 
          station_grouping = adonis(dd ~ station_grouping, data=sd, permutations = perm)
        )
      }) 
  })

# make a nice little table of anova results
# this assumes the name of the list entry is the same as the
# name of the variable being examined
anova_table <- list(anovas,animal_anovas) %>%
  map2_dfr(c("Complete dataset","Metazoans"),~{
    .x %>%
      map2_dfr(names(.),~{
        name_map <- c("sim" = "Simpson ($\\sim$)", "jaccard" = "Jaccard ($\\jac$)")
        method <- name_map[.y]
        .x %>%
          map2_dfr(names(.),~{
            marker_name <- plot_text2[.y]
            .x %>%
              map2_dfr(names(.),~{
                as_tibble(.x$aov.tab,rownames="term") %>%
                  filter(term == .y) %>%
                  select(term,pseudo_f=`F.Model`,p_value=`Pr(>F)`)
              }) %>%
              mutate(marker=marker_name)
          }) %>%
          mutate(index=method)
      }) %>%
      mutate(dataset=.y) %>%
      select(dataset,index,marker,term,everything())   
  })


# manuscript figures ------------------------------------------------------
figure_dir <- "~/projects/dissertation/manuscript/figures"

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
  fish = c("A","B","C"),
  inverts = c("D","E","F"),
  metazoans = c("G","H","I")
)

beta_map <- list(
  factor_levels = list(
    "jaccard" = c("beta.jac","beta.jtu","beta.jne"),
    "sim" = c("beta.sor","beta.sim","beta.sne")
  )
)

# create composite figure of all beta diversity heatmaps
beta_pairs_composite <- beta_pairs %>%
  map2(names(.),~{
    dm <- .y 
    .x %>%
      map2(names(.),~{
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
            plot_annotation(tag_levels = list(beta_label_map[[marker]])) & theme(plot.tag = element_text(face="bold"))
        ) 
      }) %>%
      reduce(`/`) 
  })

ggsave(path(figure_dir,"mesophotic_beta_pairs.pdf"),beta_pairs_composite$sim,device=cairo_pdf,width=12,height=9,units="in")

#### Figure: beta diversity heatmaps (animals)
beta_label_map <- list(
  inverts = c("A","B","C"),
  metazoans = c("D","E","F")
)
beta_pairs_composite_animals <- beta_pairs_animals %>%
  map2(names(.),~{
    dm <- .y 
    .x %>%
      map2(names(.),~{
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
            plot_annotation(tag_levels = list(beta_label_map[[marker]])) & theme(plot.tag = element_text(face="bold"))
        ) 
      }) %>%
      reduce(`/`) 
  })
beta_pairs_composite_animals
ggsave(path(figure_dir,"mesophotic_beta_pairs_animals.pdf"),beta_pairs_composite_animals$sim,device=cairo_pdf,width=12,height=6,units="in")


#### Figure: cluster plots
cluster_plotz <- beta_pairs$sim %>%
  map2(names(.),~{
    dd <- .x %>%
      filter(measurement == "beta.sim") %>%
      bind_rows(
        tibble(
          depth1=factor(levels(.$depth1),levels=levels(.$depth1)),
          depth2=factor(levels(.$depth1),levels=levels(.$depth1)),
          dist = 0,
          sd = 0,
          measurement = "beta.sim"
        ) 
      ) %>%
      arrange(depth1,depth2) %>%
      pivot_wider(-c(sd,measurement),names_from=depth2,values_from=dist) %>%
      column_to_rownames("depth1") %>% 
      as.matrix() %>% 
      as.dist()
    
    title <- plot_text2[.y]
    labels(dd) <- str_c(labels(dd)," ")
    ggplot(as.dendrogram(hclust(dd),hang=0.5)) +
      theme(
        plot.caption = element_text(hjust=0.5,size=14)
      ) + 
      expand_limits(y=-0.25)
  })

cluster_composite <- (cluster_plotz$fish + cluster_plotz$inverts + cluster_plotz$metazoans) +
  plot_annotation(tag_levels="A") & theme(plot.tag = element_text(face="bold"))
cluster_composite
ggsave(path(figure_dir,"mesophotic_cluster.pdf"),cluster_composite,device=cairo_pdf,width=12,height=4,units="in")

#### Figure: cluster plots (animals)
aminal_cluster_plotz <- beta_pairs_animals$sim %>%
  map2(names(.),~{
    dd <- .x %>%
      filter(measurement == "beta.sim") %>%
      bind_rows(
        tibble(
          depth1=factor(levels(.$depth1),levels=levels(.$depth1)),
          depth2=factor(levels(.$depth1),levels=levels(.$depth1)),
          dist = 0,
          sd = 0,
          measurement = "beta.sim"
        ) 
      ) %>%
      arrange(depth1,depth2) %>%
      pivot_wider(-c(sd,measurement),names_from=depth2,values_from=dist) %>%
      column_to_rownames("depth1") %>% 
      as.matrix() %>% 
      as.dist()
    
    title <- plot_text2[.y]
    labels(dd) <- str_c(labels(dd)," ")
    ggplot(as.dendrogram(hclust(dd),hang=0.5)) +
      theme(
        plot.caption = element_text(hjust=0.5,size=14)
      ) + 
      expand_limits(y=-0.25)
  })

aminal_cluster_composite <- (aminal_cluster_plotz$inverts + aminal_cluster_plotz$metazoans) +
  plot_annotation(tag_levels="A") & theme(plot.tag = element_text(face="bold"))
aminal_cluster_composite
ggsave(path(figure_dir,"mesophotic_cluster_animals.pdf"),aminal_cluster_composite,device=cairo_pdf,width=8,height=4,units="in")

#### Figure: shallow vs deep ordinations
ord_zone_plotz <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    c("depth_zone","depth_zone45") %>%
      set_names() %>%
      map(~{
        zone <- .x
        comm_ps %>%
          map2(names(.),~{
            title <- plot_text2[.y]
            p <- plot_betadisp(.x, group=zone, method=dm, list=TRUE)
            
            p$plot <- p$plot +
              scale_fill_manual(values=pal[c(1,length(pal))],name="Depth Zone") + 
              plotz_theme("light") +
              xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
              ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
              theme(legend.position = "right")
            p$plot
          })
      })
  })

# main text figure (shallow = 0-45m)
ord_zone_composite <- 
  (ord_zone_plotz$sim$depth_zone45$fish + ord_zone_plotz$sim$depth_zone45$inverts + ord_zone_plotz$sim$depth_zone45$metazoans) +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
ord_zone_composite
ggsave(path(figure_dir,"mesophotic_ord_shallowdeep45.pdf"),ord_zone_composite,device=cairo_pdf,width=12,height=4,units="in")

# supplemental figure (shallow = 0-30m)
ord_zone_composite <- 
  (ord_zone_plotz$sim$depth_zone$fish + ord_zone_plotz$sim$depth_zone$inverts + ord_zone_plotz$sim$depth_zone$metazoans) +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
ord_zone_composite
ggsave(path(figure_dir,"mesophotic_ord_shallowdeep30.pdf"),ord_zone_composite,device=cairo_pdf,width=12,height=4,units="in")

#### supplemental figure: shallow vs deep ordinations for animals
ord_zone_plotz_animals <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    c("depth_zone","depth_zone45") %>%
      set_names() %>%
      map(~{
        zone <- .x
        animals %>%
          map2(names(.),~{
            title <- plot_text2[.y]
            p <- plot_betadisp(.x, group=zone, method=dm, list=TRUE)
            
            p$plot <- p$plot +
              scale_fill_manual(values=pal[c(1,length(pal))],name="Depth Zone") + 
              plotz_theme("light") +
              xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
              ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
              theme(legend.position = "right")
            p$plot
          })
      })
  })

# composite figure
ord_zone_composite_animals <- 
  (ord_zone_plotz_animals$sim$depth_zone45$inverts + ord_zone_plotz_animals$sim$depth_zone45$metazoans) +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
ord_zone_composite_animals
ggsave(path(figure_dir,"mesophotic_ord_shallowdeep45_animals.pdf"),ord_zone_composite_animals,device=cairo_pdf,width=8,height=4,units="in")


#### Figure: depth zone ordinations
ord_plotz <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    comm_ps %>%
      map2(names(.),~{
        title <- plot_text2[.y]
        p <- plot_betadisp(.x, group="depth_f", method=dm, list=TRUE, expand=TRUE)
        p$plot <- p$plot +
          scale_fill_manual(values=pal,name="Depth") +
          plotz_theme("light") +
          xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
          ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
          theme(legend.position = "right") +
          guides(color="none")
        p$plot
      })
  })
ord_composite <- 
  (ord_plotz$sim$fish + ord_plotz$sim$inverts + ord_plotz$sim$metazoans) +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
ord_composite

ggsave(path(figure_dir,"mesophotic_ord_samples.pdf"),ord_composite,device=cairo_pdf,width=12,height=4,units="in")

### Supplemental figure: depth zone ordinations for animals
ord_plotz_animals <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    animals %>%
      map2(names(.),~{
        title <- plot_text2[.y]
        p <- plot_betadisp(.x, group="depth_f", method=dm, list=TRUE, expand=TRUE)
        p$plot <- p$plot +
          scale_fill_manual(values=pal,name="Depth") +
          plotz_theme("light") +
          xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
          ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
          theme(legend.position = "right") +
          guides(color="none")
        p$plot
      })
  })
ord_composite_animals <- 
  (ord_plotz_animals$sim$inverts + ord_plotz_animals$sim$metazoans) +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
ord_composite_animals

ggsave(path(figure_dir,"mesophotic_ord_samples_animals.pdf"),ord_composite_animals,device=cairo_pdf,width=8,height=4,units="in")

#### Figure: fish depth distributions vs eDNA detections
fish <- communities$fish$raw %>%
  inner_join(communities$fish$tax_data,by="OTU") %>%
  select(domain:species,OTU,matches(sample_pattern)) %>%
  filter(species != "dropped" & !is.na(species))

fish_species <- fish %>%
  distinct(species) %>%
  pull(species)

# make a depth detection summary
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

# get the data for the t-test
deeper <- detected_depth %>%
  filter(deeper == TRUE) %>%
  pull(reads)
shallower <- detected_depth %>%
  filter(deeper == FALSE) %>%
  pull(reads)

# do the t-test
# let's do the little statistical tests that Muff et al 2022 did
t.test(deeper,shallower,paired=FALSE,alternative="two.sided",var.equal = FALSE)

dd <- detected_depth %>%
  filter(depth > fb_deep | depth < fb_shallow) %>%
  mutate(
    depth_diff = case_when(
      depth < fb_shallow ~ abs(depth-fb_shallow),
      TRUE ~ abs(depth-fb_deep)
    )
  )

cor.test(dd$depth_diff,dd$reads,method="kendall",exact = F)

# finish making the depth detection summary
detected_depth <- detected_depth %>%
  filter(reads >= 10) %>%
  group_by(family,species,depth,fb_deep,fb_shallow) %>%
  summarise(n=n(),reads=sum(reads)) %>%
  ungroup() %>%
  filter(fb_deep < 150)  %>%
  ungroup() %>%
  group_by(species) %>%
  mutate(depth_diff = depth - max(fb_deep)) %>%
  ungroup() %>%
  mutate(
    depth_range = fb_deep-fb_shallow,
    species = fct_reorder(species,depth_range,.desc = TRUE)
    # species = fct_reorder(species,as.numeric(factor(family)))
  ) #%>%
  # select(all_of(names(dd)))
  
# what is the most number of observations for a species?
most <- detected_depth %>%
  count(species) %>% 
  pull(n) %>% 
  max()

# this is a little hack to make the depth ranges look good
# we basically make them be drawn an equal number of times per
# species, so the lines are all equally thick.
detected_depth <- detected_depth %>%
  group_by(species) %>%
  group_modify(~{
    grp <- .x
    to_add <- most-nrow(grp)
    if (to_add > 0) {
      add <- map_dfr(seq(to_add),~grp %>% slice(1) %>% mutate(reads=NA))
      grp <- grp %>%
        bind_rows(add)
    }
    return(grp)
  }) 

depth_plotz <- ggplot(detected_depth) + 
  geom_errorbar(aes(x=species,ymin=fb_shallow,ymax=fb_deep),width=0.5,color="dodgerblue4") +
  geom_point(aes(x=species,y=depth,size=reads,color=recode_factor(factor(str_c(depth,'m')), "0m" = "Surface")),alpha=0.7) +
  stat_summary(aes(x=species,y=depth),geom="point",fun="mean",col="black",size=10,shape="-") + 
  scale_color_manual(values=pal,name="Depth zone") +
  scale_size(name="Read count",labels=scales::comma,range=c(2,12),breaks=c(1e2,1e3,1e4,5e4,1e5)) +
  scale_y_reverse(breaks=seq(0,300,by=30),minor_breaks=seq(0,300,by=15)) +
  scale_x_discrete(position="top") + 
  plotz_theme("light") +
  theme(
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
ggsave(path(figure_dir,"mesophotic_fish_depth.pdf"),depth_plotz,device=cairo_pdf,width=12,height=10,units="in")

#### Supplemental Figure: ordinations by site
ord_sites <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    comm_ps %>%
      map2(names(.),~{
        title <- plot_text2[.y]
        p <- plot_betadisp(.x, group="station_grouping", method=dm, list=TRUE)
        p$plot <- p$plot +
          scale_fill_manual(values=pal[c(1,4,7)],name="Site") +
          plotz_theme("light") + 
          xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
          ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
          theme(legend.position = "right")
        p$plot
      })
  })

ord_site_composite <- 
  (ord_sites$sim$fish + ord_sites$sim$inverts + ord_sites$sim$metazoans) +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
ord_site_composite
ggsave(path(figure_dir,"mesophotic_ord_sites.pdf"),ord_site_composite,device=cairo_pdf,width=12,height=4,units="in")

ord_site_composite_jaccard <- 
  (ord_sites$jaccard$fish + ord_sites$jaccard$inverts + ord_sites$jaccard$metazoans) +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
ord_site_composite_jaccard
ggsave(path(figure_dir,"mesophotic_ord_sites_jaccard.pdf"),ord_site_composite_jaccard,device=cairo_pdf,width=12,height=4,units="in")

#### Supplemental Figure: species accumulations
sup <- function(...) suppressWarnings(suppressMessages(...))
all_models <- comm_ps %>%
  map2(names(.),~{
    sd <- sample_tibble(.x)
    specs <- otu_tibble(.x) %>%
      inner_join(sd %>% select(sample,depth_f),by="sample") %>%
      select(sample,depth_f,everything()) %>%
      nest(otu_table=-depth_f) %>%
      mutate(otu_table = map(otu_table,~{
        sup(curve <- .x %>%
              column_to_rownames("sample") %>% 
              as("matrix") %>%
              decostand("pa",margin=NULL) %>%
              specaccum(method="random",permutations=1000))
        as_tibble(curve$perm) %>%
          mutate(sites=row_number()) %>%
          select(sites,everything()) %>%
          pivot_longer(-sites,names_to = "permutation", values_to = "richness") %>%
          select(sites,richness)
      })) %>%
      unnest(otu_table)
    
    try_mod <- function(e) {
      tryCatch(
        e,
        error = function(e) NULL
      )
    }
    
    modls <- list(
      asymp = nls(richness ~ SSasymp(sites, Asym, R0, lrc),  data=specs),
      gompertz = nls(richness ~ SSgompertz(sites, Asym,  xmid, scal), data=specs),
      `michaelis-menten` = nls(richness ~  SSmicmen(sites, Vm, K), data=specs),
      logis = nls(richness ~ SSlogis(sites,  Asym, xmid, scal), data=specs) 
    )
    AICs <- sort(map_dbl(modls,~{
      if (!is_null(.x))
        AIC(.x)
      else
        Inf
    }))
    best_mod <- names(AICs)[1]
    pred <- modls[[best_mod]]
    predictions <- tibble(
      x = seq(1:2000),
      y = predict(pred,list(sites=seq(1:2000)))
    )
    
    slopes <- diff(predictions$y)
    slopes <- c(1,slopes)
    asymp <- predictions$x[slopes < 1][1]
    asymp_taxa <- predictions$y[slopes < 1][1]
    taxa_6 <- predictions$y[6]
    percent_taxa <- taxa_6/asymp_taxa
    
    
    cat("Community:",.y,"\n")
    cat("Best model:",best_mod,"\n")
    cat("Asymptote at:",asymp,"replicates\n")
    cat("Taxa at asymptote:",asymp_taxa,"\n")
    cat("Taxa at 6 replicates:",taxa_6,"\n")
    cat("Percent of taxa represented by 6 replicates:",scales::percent(percent_taxa),"\n")
    cat("\n\n")
    list(models = modls, best_model = best_mod, predictions = predictions, asymptote=asymp, asymptote_taxa=asymp_taxa, our_taxa = taxa_6, accum_data=specs)
  })

accum_plotz <- all_models %>%
  map2(names(.),~{
    d <- .x$accum_data %>%
      group_by(depth_f,sites) %>%
      summarise(sd=sd(richness),richness=mean(richness))
    p <- .x$predictions
    tit <- plot_text2[.y]
    ggplot() + 
      geom_line(data=d,aes(x=sites,y=richness,color=depth_f)) + 
      geom_errorbar(data=d,aes(x=sites,ymin=richness-sd,ymax=richness+sd),width=0.5) + 
      geom_line(data=p,aes(x=x,y=y),color="blue",size=1.1) + 
      geom_vline(xintercept=.x$asymptote,color="red") + 
      scale_color_manual(values=pal,name="Depth Zone") + 
      xlim(1,25) +
      xlab("Replicates") + 
      ylab("zOTUs") + 
      theme_bw() 
  })
accum_composite <- 
  accum_plotz$fish + accum_plotz$inverts + accum_plotz$metazoans + 
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
accum_composite
ggsave(path(figure_dir,"mesophotic_accum.pdf"),accum_composite,device=cairo_pdf,width=12,height=4,units="in")

#### Figure (supplemental?) shallow/deep shared zOTU venn diagrams
venn_data <- comm_ps %>%
  map(~{
    .x %>%
      psmelt() %>%
      # mutate(depth_zone45 = cut(depth,c(-Inf,31,61,Inf),labels=c("Shallow","Mid","Deep"))) %>%
      group_by(OTU,depth_zone45) %>%
      summarise(present = as.integer(sum(Abundance) > 0)) %>%
      pivot_wider(names_from="depth_zone45",values_from="present") %>%
      ungroup()
  })
venn_pal <- c(pal[1],pal[length(pal)],pal[(length(pal)+1)/2])
venn_plotz <- venn_data %>%
  map(~{
    shallow <- .x %>% filter(Shallow == 1) %>% pull(OTU)
    deep <- .x %>% filter(Deep == 1) %>% pull(OTU)
    # mid <- .x %>% filter(Mid == 1) %>% pull(OTU)
    f <- list("Shallow\n(0–30m)"=shallow,"Deep\n(75–90m)"=deep)
    v <- Venn(f)
    p <- process_data(v)
    names(venn_pal) <- p@region$name
    p@region <- p@region %>%
      mutate(count_label = str_glue("{count}\n({scales::percent(count/sum(count))})"))
    ggplot() + 
      geom_sf(aes(fill=name),data = p@region, alpha=0.7) + 
      geom_sf(color="black", size = 1,data = p@setEdge, show.legend = F) +
      geom_sf_text(aes(label = name), data = p@setLabel,size=4, nudge_x=c(-50,50)) + 
      # geom_sf_label(aes(label = name), data = p@setLabel,size=5) + 
      geom_sf_text(aes(label=count_label), fontface = "bold", family = "serif", size = 6, color = "grey3",  data = p@region) +
      scale_fill_manual(values=venn_pal) + 
      expand_limits(y=c(225,815)) +
      theme_void() +
      guides(fill="none")  #+ 
      # theme(panel.border = element_rect(fill=NA))
    # ggvenn(list("Shallow\n(0–45m)"=shallow,"Deep\n(60–90m)"=deep)) + 
    #   theme(text = element_text(size=2))
  })

venn_composite <- venn_plotz %>%
  reduce(`+`) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides="collect") &
  theme(plot.tag = element_text(face="bold"), plot.caption = element_text(hjust=0.5,size=16))
venn_composite
ggsave(path(figure_dir,"mesophotic_venn.pdf"),venn_composite,device=cairo_pdf,width=12,height=3.5,units="in")

### supplemental figure: 16S zotu upset plot
`
# function to make an upset plot from a presence-absen`ce matrix
# of category occurrences 
plot_upset <- function(dataset,
                       name_column,
                       data_columns,
                       dot_size = 6,
                       line_size=2,
                       prefilter=F,
                       intersects=NA,
                       min_intersects=0,
                       bar_lab = "intersections",
                       sidebar_lab = "number in category",
                       label_top_bars = FALSE,
                       label_side_bars = FALSE,
                       group_palette = NULL) {
  if (prefilter) {
    dataset <- dataset %>%
      mutate(
        across({{data_columns}},~if_else(.x > 0,1,0,missing=0))
      )
  }
  
  # get category names
  sets <- dataset %>%
    select({{data_columns}}) %>%
    names()
  
  dataset <- dataset %>%
    unite("code",{{data_columns}},sep="",remove = FALSE) %>%
    mutate(
      code = str_replace_all(code,"0","n"),
      code = str_replace_all(code,"1","p")
    )
  dataset_long <- dataset %>%
    pivot_longer({{data_columns}},names_to = "metric_name", values_to = "metric")
  
  data1 <- dataset_long %>%
    group_by(code,metric_name,metric) %>%
    summarise(n = n_distinct(!!sym(name_column))) %>%
    arrange(n) %>%
    ungroup()
  data2 <- dataset_long %>%
    group_by(metric_name) %>%
    summarise(n=sum(metric))
  
  x_breaks <- dataset_long %>%
    group_by(code) %>%
    summarise(n=n_distinct(!!sym(name_column))) %>%
    arrange(desc(n)) %>%
    filter(n > min_intersects) %>%
    pull(code) %>%
    as.character()
  
  if (!is.na(intersects)) {
    x_breaks <- x_breaks[1:intersects]
  }
  
  intersects <- length(x_breaks)
  
  top_bars <-
    ggplot(data1, aes(x=reorder(factor(code),-n), y=n)) +
    geom_col(fill="grey5", position="dodge") +
    scale_x_discrete(limits=x_breaks) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    theme_bw() +
    ylab(bar_lab) +
    theme(legend.position = "none", 
          # axis.title = element_blank(),
          axis.line.y.left = element_line(),
          axis.ticks.x = element_blank(),
          panel.border = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(0,0,0,0,"cm"))
  
  if (!is.null(group_palette)) {
    fill_vals <- group_palette
    names(fill_vals) <- sets
  } else {
    fill_vals <- rep("grey5",length(sets))
    names(fill_vals) <- sets
  }
  
  side_bars <-
    ggplot(data2, aes(x=metric_name, y=n)) +
    geom_col(aes(fill=metric_name), position="dodge") +
    scale_fill_manual(values=fill_vals) + 
    scale_x_discrete(
      position="top",
      limits=rev(sets)
    ) +
    scale_y_reverse(
      labels=scales::format_format(big.mark = ",", decimal.mark = ".", scientific = FALSE, digits=0),
      expand = expansion(mult = c(0.6, 0))
    ) +
    coord_flip() +
    theme_bw() +
    theme(
      axis.line.x = element_line(),
      panel.border = element_blank(),
      legend.position = "none",
      axis.title = element_blank(),
      axis.text.x = element_text(size=6),
      axis.text.y = element_text(size=12),
      panel.grid = element_blank(),
      axis.title.x = element_text(),
      plot.margin = margin(0,0,0,0,"cm")
    ) +
    ylab(sidebar_lab) 
  
  if (label_side_bars) { 
    side_bars <- side_bars + 
      geom_text(aes(label=n), position = position_dodge(0.9),  hjust=1.1, vjust=0.5) 
  }
  
  if (label_top_bars) {
    top_bars <- top_bars +  
      geom_text(aes(label=n), position = position_dodge(0.9), hjust=0.5, vjust=-0.25)
  }

 
  
  dot_lines <- data1 %>%
    filter(metric == 1) %>%
    group_by(code) %>%
    summarise(
      f = list(factor(metric_name,levels=sets)),
      n = unique(n)
    ) %>%
    mutate(
      start = map_chr(f,~sets[min(as.integer(.x))]),
      end = map_chr(f,~sets[max(as.integer(.x))])
    ) %>%
    select(-f)
  
  cols <- c("0" = "grey77", "1" = "grey2")  
  data1 <- data1 %>%
      mutate(color_group = str_glue("{metric_name}_{metric}"))
  cols <- data1 %>%
    ungroup() %>%
    mutate(
      color_value = fill_vals[metric_name],
      color_value = case_when(
        metric == 0 ~ "grey95",
        TRUE ~ color_value
      )
    ) %>%
    distinct(color_group,color_value) %>%
    deframe()
  dots <-
    ggplot(data1, aes(y=metric_name, x=reorder(factor(code),-n))) +
    geom_point(shape=21, size=dot_size, colour="black", aes(fill=color_group)) +
    geom_point(data=data1 %>% filter(metric == 1),shape=19,size=dot_size/2,color="black") + 
    geom_segment(data=dot_lines,mapping=aes(x=reorder(factor(code),-n),xend=reorder(factor(code),-n),y=start,yend=end),size=line_size) +
    scale_fill_manual(values = cols) +
    scale_color_manual(values=fill_vals) + 
    scale_x_discrete(limits=x_breaks) +
    scale_y_discrete(limits=rev(sets)) +
    theme_minimal() +
    labs(x="", y="") +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank(),
          plot.margin = margin(0,0,0,0,"cm"))
  layout <- "
    #1111
    #1111
    #1111
    23333
    23333
  "
  return(top_bars + side_bars + dots + plot_layout(design=layout))
}

upset_data <- comm_ps %>%
  map(~{
    .x %>%
      psmelt() %>%
      group_by(OTU,depth_f) %>%
      summarise(present = as.integer(sum(Abundance) > 0)) %>%
      pivot_wider(names_from="depth_f",values_from="present") %>%
      ungroup()
  })

upset_plotz <- upset_data %>%
  map(~{
    .x %>%
      plot_upset("OTU",
                 Surface:`90m`,
                 bar_lab="Shared zOTUs",
                 sidebar_lab="zOTUs detected",
                 label_side_bars = TRUE,
                 # label_top_bars = TRUE,
                 group_palette = pal,
                 intersects = 30,
                 dot_size = 6,
                 line_size=2)
  })


upset_composite <- upset_plotz %>%
  map(wrap_elements) %>%
  reduce(`/`) + 
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
upset_composite 
# upset_plotz$fish

ggsave(path(figure_dir,"mesophotic_upset.pdf"),upset_composite,device=cairo_pdf,width=9,height=14,units="in")

# manuscript tables -------------------------------------------------------
table_dir <- "~/projects/dissertation/manuscript/tables"

write_ms_table <- function(tbl,file,caption="",bold_header=TRUE,...) {
  if (bold_header) {
    tbl <- tbl %>%rename_with(~str_c("**",.x,"**"))
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
  mutate(pseudo_f=round(pseudo_f,2), p_value=round(p_value,2)) %>%
  mutate( term = term_map[term] ) %>%
  rename(Dataaset=dataset,`Dissimilarity Index`=index,`Assay`=marker,`Term`=term,`Pseudo-F`=pseudo_f,`p-value`=p_value) %>%
  write_ms_table(path(table_dir,"mesophotic_anovas.csv"),"Results of PERMANOVA analyses",bold_header = TRUE)

#### Supplemental table: eDNA reads summary
all_samples <- sample_data %>%
  select(`Sample ID` = id,substrate,station_grouping,depth_f)

normed <- comm_ps %>%
  map2_dfr(names(.),~{
    otu_tibble(.x) %>%
      column_to_rownames("sample") %>%
      rowSums() %>%
      enframe(name = "Sample ID", value = "Normalized sequence reads") %>%
      mutate(Assay=plot_text[.y])
  })

reads_summary <- 
  communities %>%
  map2_dfr(names(.),~{
    gene <- markers %>%
      filter(description == .y) %>%
      pull(gene)
    
    counts <- enframe(
        .x$unfiltered %>%
          select(where(is.numeric),-matches("Blast")) %>%
          colSums(),
        name = "Sample ID",
        value = "Sequence reads"
      )
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

write_ms_table(rs,path(table_dir,"mesophotic_reads_summary.csv"),"Sequencing results for mesophotic eDNA samples across assay types",na="",bold_header = TRUE)

# numbers reported in the text --------------------------------------------

#### Beta diversity summaries
# pairwise beta diversity summary
beta_pairs$sim %>%
  map(~{
    .x %>%
      group_by(measurement) %>%
      summarise(mean=round(mean(dist),2),sd=round(sd(dist),2))
  })

# pairwise beta diversity summary (animals)
beta_pairs_animals$sim %>%
  map(~{
    .x %>%
      group_by(measurement) %>%
      summarise(mean=round(mean(dist),2),sd=round(sd(dist),2))
  })

#### eDNA read counts for different categories squished together
cat("sample:\n")
walk2(comm_ps,names(comm_ps),~{
  ss <- sample_sums(.x)
  cat(.y,": mean: ",scales::number(mean(ss),big.mark = ",")," sd: ",scales::number(sd(ss),big.mark = ","),"\n")
})
cat("site:\n")
walk2(comm_ps,names(comm_ps),~{
  f <- merge_samples(.x,"station")
  ss <- sample_sums(f)
  cat(.y,": mean: ",scales::number(mean(ss),big.mark = ",")," sd: ",scales::number(sd(ss),big.mark = ","),"\n")
})

cat("depth zone:\n")
walk2(comm_ps,names(comm_ps),~{
  f <- merge_samples(.x,"depth_f")
  ss <- sample_sums(f)
  cat(.y,": mean: ",scales::number(mean(ss),big.mark = ",")," sd: ",scales::number(sd(ss),big.mark = ","),"\n")
})


#### Indicator species analyses

# helper function: does IndVal analysis for a phyloseq object across some factor
get_indicators <- function(community,variable,fdr=0.10,permutations=999,list=TRUE) {
  # get sample data
  sd <- sample_tibble(community)
  # get presence/absence-transformed community matrix
  mat <- community %>%
    otu_table() %>%
    as("matrix") %>%
    decostand("pa")
  
  # yank analysis factor
  categories <- sd %>% pull(variable)
  # do the IndVal analysis
  iv <- multipatt(mat,categories,control=how(nperm=permutations))
  
  # multipatt returns a thing with columns called s.<variable_category>
  # here we're just enumerating those
  sign_cats <- str_c("s.",unique(categories))
  # construct a table with taxa, Indval, p value, and category
  d <- as_tibble(iv$sign,rownames="zotu") %>% 
    pivot_longer(all_of(sign_cats),names_to="group",values_to="active") %>%
    mutate(group = str_replace(group,"^s\\.","")) %>% # strip off the "s." at the beginning of the group name
    filter(active == 1) %>% # get get only values that apply to the group in question
    drop_na(p.value) %>%  # drop anything with an NA p value
    rename(p_value=p.value) %>%
    # do p value correction for multiple comparisons using the false discovery rate method
    mutate(p_value = suppressWarnings(fdrtool(p_value,statistic="pvalue",plot=FALSE,verbose=FALSE,cutoff.method="pct0",pct0=fdr)$pval)) %>%
    select(zotu,group,stat,p_value) %>%
    arrange(zotu) %>%
    left_join(community %>% taxa_tibble(),by="zotu")
  # return the results of the multipatt call and the table we've made
  # to make it easier to interpret
  if (list) {
    return(list(iv=iv,table=d))
  } else {
    return(d)
  }
}

## all taxa
# do indicator analysis for deep v shallow
indicators <- animals %>% 
  map(~{
    get_indicators(.x,variable="depth_zone45",list=FALSE)
  })

itable <- indicators %>%
  enframe(name = "assay", value="indval") %>%
  unnest(indval) %>%
  filter(p_value < 0.05) %>%
  select(assay,group,stat,p_value,domain:species) %>%
  mutate(assay = plot_text2[assay]) %>%
  rename(Assay=assay,Depth=group,IndVal=stat,`p-value`=p_value) %>%
  rename_with(str_to_title,domain:species) %>%
  arrange(Assay,desc(Depth),Domain,Kingdom,Phylum,Class,Order,Family,Genus,Species) %>%
  mutate(
    Depth = case_when(
      Depth == "Shallow" ~ "Shallow\n(0--45m)",
      Depth == "Deep" ~ "Deep\n(60--90m)"
    )
  ) %>%
  mutate(Species = str_c("*",Species,"*"))
itable
write_ms_table(itable,path(table_dir,"mesophotic_indicators.csv"),caption = "Significant indicator species (*IndVal*) analysis results",na="",bold_header = TRUE)
# indval_shallowdeep <- comm_ps %>%
#   map(get_indicators,variable="depth_zone45")
# # do indicator analysis for individual depth zones
# indval_depthzone <- comm_ps %>%
#   map(get_indicators,variable="depth_f")
# ## animals
# # do indicator analysis for deep v shallow
# indval_shallowdeep_animals <- animals %>%
#   map(get_indicators,variable="depth_zone45")
# # do indicator analysis for individual depth zones
# indval_depthzone_animals <- animals %>%
#   map(get_indicators,variable="depth_f")
