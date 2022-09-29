library(here)
library(EcolUtils)
library(beyonce)
library(plotly)
library(indicspecies)
library(betapart)
library(Polychrome)
library(patchwork)
library(dendextend)
library(rfishbase)
library(fs)
library(vegan)


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
wisco <- FALSE

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
    # table(b$domain,b$phylum)
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
minimals <- 1000                                      # minimum reads to retain a sample
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
} else {
  animals <- comm_ps %>%
    map(~{
      f <- subset_taxa(.x,kingdom == "Metazoa")         # filter by metazoans
      f <- prune_samples(sample_sums(f) >= minimals,f)  # filter by reads >= minimum
      new_otu <- otu_table(f) %>%                       # re-rarefy new otu table to equal depth
        as("matrix") 
      if (rarefy) {
        new_otu <- new_otu %>%
          rrarefy.perm(n=rarefy_perm)
        otu_table(f) <- otu_table(new_otu,taxa_are_rows = FALSE) # reassign new rarefied otu table
      }
      return(f)
    })
}
# animals <- comm_ps %>%
#   map(~{
#     f <- subset_taxa(.x,kingdom == "Metazoa")         # filter by metazoans
#     f <- prune_samples(sample_sums(f) >= minimals,f)  # filter by reads >= minimum
#     new_otu <- otu_table(f) %>%                       # re-rarefy new otu table to equal depth
#       as("matrix") 
#     if (rarefy) {
#       new_otu <- new_otu %>%
#         rrarefy.perm(n=rarefy_perm)
#       otu_table(f) <- otu_table(new_otu,taxa_are_rows = FALSE) # reassign new rarefied otu table
#     }
#     return(f)
#   })

to_plot <- comm_ps

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
beta_pairs <- comm_ps %>%
  map(~{
    sd <- sample_tibble(.x)
    .x %>%
      otu_table() %>%
      as("matrix") %>%
      decostand("pa") %>%
      beta.pair("sorensen") %>%
      map2_dfr(names(.),~{
        .x %>%
          as("matrix") %>%
          as_tibble(m,rownames="row") %>%
          pivot_longer(-row,names_to="col",values_to="dist") %>%
          left_join(sd %>% select(sample,depth_f),by=c("row" = "sample")) %>%
          left_join(sd %>% select(sample,depth_f),by=c("col" = "sample"),suffix = c("_s1","_s2")) %>%
          mutate(row=factor(row),col=factor(col)) %>%
          # filter(as.numeric(row) > as.numeric(col)) %>%
          filter(as.numeric(depth_f_s1) > as.numeric(depth_f_s2)) %>%
          select(sample1=row,sample2=col,depth1=depth_f_s1,depth2=depth_f_s2,dist) %>%
          group_by(depth1,depth2) %>%
          summarise(sd=sd(dist),dist=mean(dist)) %>%
          mutate(measurement = .y) %>%
          select(depth1,depth2,dist,sd,measurement)
      })
  })

# pairwise (by depth zone) beta diversity metrics (animals only)
beta_pairs_animals <- animals %>%
  map(~{
    .x %>%
      otu_table() %>%
      as("matrix") %>%
      decostand("pa") %>%
      beta.pair("sorensen") %>%
      map2_dfr(names(.),~{
        .x %>%
          as("matrix") %>%
          as_tibble(m,rownames="row") %>%
          pivot_longer(-row,names_to="col",values_to="dist") %>%
          left_join(sample_data %>% select(id,depth_f),by=c("row" = "id")) %>%
          left_join(sample_data %>% select(id,depth_f),by=c("col" = "id"),suffix = c("_s1","_s2")) %>%
          mutate(row=factor(row),col=factor(col)) %>%
          # filter(as.numeric(row) > as.numeric(col)) %>%
          filter(as.numeric(depth_f_s1) > as.numeric(depth_f_s2)) %>%
          select(sample1=row,sample2=col,depth1=depth_f_s1,depth2=depth_f_s2,dist) %>%
          group_by(depth1,depth2) %>%
          summarise(sd=sd(dist),dist=mean(dist)) %>%
          mutate(measurement = .y) %>%
          select(depth1,depth2,dist,sd,measurement)
      })
  })

# overall beta diversity metrics
beta_diversity <- map(comm_ps,~{
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
}) %>% set_names(names(to_plot))

beta_diversity_animals <- map(animals,~{
  .x <- merge_samples(.x,"depth_f")
  comm <- as(otu_table(.x),"matrix")
  # do the abundance-based beta partitioning
  # bc <- betapart.core.abund(comm)
  comm <- decostand(comm,"pa")
  # now the presence-absence one
  bp <- betapart.core(comm)


  sorensen <- beta.multi(bp,"sorensen")
  jaccard <- beta.multi(bp,"jaccard")

  list(sorensen=sorensen,jaccard=jaccard)
}) %>% set_names(names(to_plot))

# start here --------------------------------------------------------------


# beta diversity summaries ------------------------------------------------

# pairwise beta diversity summary
beta_pairs %>%
  map(~{
    .x %>%
      group_by(measurement) %>%
      summarise(mean=round(mean(dist),2),sd=round(sd(dist),2))
  })

# pairwise beta diversity summary (animals)
beta_pairs_animals %>%
  map(~{
    .x %>%
      group_by(measurement) %>%
      summarise(mean=round(mean(dist),2),sd=round(sd(dist),2))
  })

# permanova analyses ------------------------------------------------------
anovas <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    map2(to_plot,names(to_plot),~{
      perm <- 8000
      sd <- sample_tibble(.x)
      dd <-distance(.x,method=dm)
      list(
        # all = adonis(dd~ depth_zone + depth_zone45 + depth_f + station_grouping, data=sd, permutations=perm),
        depth_zone = adonis(dd ~ depth_zone, data=sd, permutations=perm),
        # adonis_shallowdeep_st = adonis(dd ~ depth_zone + station_grouping, data=sd, permutations=perm),
        depth_zone45 = adonis(dd ~ depth_zone45, data=sd, permutations=perm),
        # adonis_shallowdeep45_st = adonis(dd ~ depth_zone45 + station_grouping, data=sd, permutations=perm),
        # adonis_sites = adonis(dd ~ depth_f*station_grouping, data=sd, permutations = perm), 
        # adonis_shallowdeep_sites = adonis(dd ~ depth_zone*station_grouping, data=sd, permutations=perm),
        adonis_depth = adonis(dd ~ depth_f, data=sd, permutations = perm), 
        # adonis_depth_st = adonis(dd ~ depth_f + station_grouping, data=sd, permutations = perm), 
        adonis_sites = adonis(dd ~ station_grouping, data=sd, permutations = perm)
        # adonis_everything = adonis(dd ~ depth_f + station_grouping + depth_zone + depth_zone45, data=sd, permutations=perm),
        # anosim_shallowdeep = anosim(dd,sd$depth_zone, permutations = perm), 
        # anosim = anosim(dd, sd$depth_f, permutations = perm),
        # indicators_shallowdeep = multipatt(community, sd$depth_zone, control = how(nperm=perm)),
        # indicators = multipatt(community, sd$depth_f, control = how(nperm=perm)),
        # simper = simper(community, sd$depth_f, permutations = perm),
        # pairwise = adonis.pair(dd, sd$depth_f, nper = perm, corr.method = "bonferroni")
        # pairwise = pairwise_adonis(dd, sd$depth_f, permutations = perm,correction = 'BH')
      )
    }) %>% set_names(names(to_plot))
  })






# manuscript figures ------------------------------------------------------
figure_dir <- "~/projects/dissertation/manuscript/figures"

#### Figure: beta diversity heatmaps (everything)
# map distance types to pretty names
beta_title_map <- c(
  beta.sim = "Turnover",
  beta.sne = "Nestedness",
  beta.sor = "Overall\n(Sørensen)"
)
beta_label_map <- list(
  fish = c("A","B","C"),
  inverts = c("D","E","F"),
  metazoans = c("G","H","I")
)

# TODO: figure out how to make plot labels look right
# (it's already been done for the figure in the MS folder so I did it at least once)
beta_pairs_composite <- beta_pairs %>%
  map2(names(.),~{
    # .x %>%
    #   mutate(
    #     measurement = factor(measurement,levels=c("beta.sor","beta.sim","beta.sne")),
    #     depth1 = fct_rev(depth1)
    #   ) %>%
    #   group_by(measurement) %>%
    #   group_map(~{
    #     measurement <- as.character(.y$measurement)
    #     ggplot(.x) + 
    #       geom_tile(aes(x=depth2,y=depth1,fill=dist)) +
    #       scale_fill_viridis(option="inferno",name=beta_title_map[measurement]) +
    #       scale_x_discrete(position="top") + 
    #       theme_bw() +
    #       theme(
    #         axis.text.x = element_text(angle=25,vjust=1,hjust=0),
    #         panel.grid = element_blank()
    #       ) +
    #       labs(x="Depth zone",y="Depth zone")
    #   }) %>%
    #   reduce(`+`) 
    wrap_elements(
      .x %>%
        mutate(
          measurement = factor(measurement,levels=c("beta.sor","beta.sim","beta.sne")),
          depth1 = fct_rev(depth1)
        ) %>%
        group_by(measurement) %>%
        group_map(~{
          measurement <- as.character(.y$measurement)
          ggplot(.x) + 
            geom_tile(aes(x=depth2,y=depth1,fill=dist)) +
            scale_fill_viridis(option="inferno",name=beta_title_map[measurement]) +
            scale_x_discrete(position="top") + 
            theme_bw() +
            theme(
              axis.text.x = element_text(angle=25,vjust=1,hjust=0),
              panel.grid = element_blank()
            ) #+
            # labs(x="Depth zone",y="Depth zone")
        }) %>%
        reduce(`+`) 
    ) + plot_annotation(tag_levels = list(beta_label_map[[.y]]))
  }) %>%
  reduce(`/`) #+
  # plot_annotation(tag_levels = "A")
beta_pairs_composite
ggsave(path(figure_dir,"mesophotic_beta_pairs.pdf"),beta_pairs_composite,device=cairo_pdf,width=12,height=9,units="in")

#### Figure: beta diversity heatmaps (animals)
beta_pairs_composite_animals <- beta_pairs_animals %>%
  map2(names(.),~{
    .x %>%
      mutate(
        measurement = factor(measurement,levels=c("beta.sor","beta.sim","beta.sne")),
        depth1 = fct_rev(depth1)
      ) %>%
      group_by(measurement) %>%
      group_map(~{
        measurement <- as.character(.y$measurement)
        ggplot(.x) + 
          geom_tile(aes(x=depth2,y=depth1,fill=dist)) +
          scale_fill_viridis(option="inferno",name=beta_title_map[measurement]) +
          scale_x_discrete(position="top") + 
          theme_bw() +
          theme(
            axis.text.x = element_text(angle=25,vjust=1,hjust=0),
            panel.grid = element_blank()
          ) +
          labs(x="Depth zone",y="Depth zone")
      }) %>%
      reduce(`+`) 
  }) %>%
  reduce(`/`) +
  plot_annotation(tag_levels = "A")
beta_pairs_composite_animals
ggsave(path(figure_dir,"mesophotic_beta_pairs_animals.pdf"),beta_pairs_composite_animals,device=cairo_pdf,width=12,height=9,units="in")


#### Figure: cluster plots
cluster_plotz <- beta_pairs %>%
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
  plot_annotation(tag_levels="A")
cluster_composite
ggsave(path(figure_dir,"mesophotic_cluster.pdf"),cluster_composite,device=cairo_pdf,width=12,height=4,units="in")
# this is the *old* way of doing cluster plots
# merge_factor <- "depth_f"
# cluster_plotz <- distance_methods %>%
#   set_names() %>%
#   map(~{
#     dm <- .x
#     to_plot %>%
#       map2(names(.),~{
#         title <- plot_text2[.y]
#         merged <- merge_samples(.x,merge_factor)
#         comm_merged <- distance(merged,method=dm)
#         labels(comm_merged) <- str_c(labels(comm_merged)," ")
#         ggplot(as.dendrogram(hclust(comm_merged),hang=0.5)) +
#           theme(
#             plot.caption = element_text(hjust=0.5,size=14)
#           ) + 
#           expand_limits(y=-0.25)
#       }) 
#   })


#### Figure: cluster plots (animals)
aminal_cluster_plotz <- beta_pairs_animals %>%
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
  plot_annotation(tag_levels="A")
aminal_cluster_composite
ggsave(path(figure_dir,"mesophotic_cluster_animals.pdf"),aminal_cluster_composite,device=cairo_pdf,width=8,height=4,units="in")

# old way of doing clustering
# merge_factor <- "depth_f"
# aminal_cluster_plotz <- distance_methods %>%
#   set_names() %>%
#   map(~{
#     dm <- .x
#     animals %>%
#       map2(names(.),~{
#         title <- plot_text2[.y]
#         merged <- merge_samples(.x,merge_factor)
#         comm_merged <- distance(merged,method=dm)
#         labels(comm_merged) <- str_c(labels(comm_merged)," ")
#         ggplot(as.dendrogram(hclust(comm_merged),hang=0.5)) +
#           theme(
#             plot.caption = element_text(hjust=0.5,size=14)
#           ) + 
#           expand_limits(y=-0.25)
#       }) 
#   })


#### Figure: shallow vs deep ordinations
# TODO: put permanova results on ordination plots?
zone_groupings <- c("fish" = "depth_zone", "inverts" = "depth_zone45", "metazoans" = "depth_zone45")
ord_zone_plotz <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    # models <- anovas[[.x]]
    c("depth_zone","depth_zone45") %>%
      set_names() %>%
      map(~{
        zone <- .x
        map2(to_plot,names(to_plot),~{
          title <- plot_text2[.y]
          # model <- models[[.y]][[zone]]
          p <- plot_betadisp(.x, group=zone, method=dm, list=TRUE)
          
          # get plot ranges
          # xr <- layer_scales(p$plot)$x$range$range
          # yr <- layer_scales(p$plot)$y$range$range
          
          # permanova <- as_tibble(model$aov.tab,rownames = "factor") %>%
          #   rename(df=2,sum_sq=3,mean_sq=4,f_stat=5,r2=6,pval=7) %>%
          #   filter(factor == zone) %>%
          #   mutate(
          #     across(is.numeric,~round(.x,2)),
          #     str=str_glue("Pseudo-F = {f_stat}, {pval(pval)}\n\n"),
          #     x = max(xr) - diff(xr) * 0.01,
          #     y = max(yr)# - diff(yr) * 0.1
          #   )
          
          p$plot <- p$plot +
            scale_fill_manual(values=pal[c(1,length(pal))],name="Depth Zone") + #,labels=c("Shallow","Deep")) +
            plotz_theme("light") +
            xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
            ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
            # geom_text(data=permanova,mapping=aes(x=x,y=y,label=str),hjust="right",vjust="top") + 
            theme(legend.position = "right")
          p$plot
        }) %>% set_names(names(to_plot))
      })
  })

# main text figure (shallow = 0-30m)
ord_zone_composite <- 
  (ord_zone_plotz$sim$depth_zone$fish + ord_zone_plotz$sim$depth_zone$inverts + ord_zone_plotz$sim$depth_zone$metazoans) +#/  (ord_zone_plotz$sim$depth_zone45$fish + ord_zone_plotz$sim$depth_zone45$inverts + ord_zone_plotz$sim$depth_zone45$metazoans)  +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
ord_zone_composite
ggsave(path(figure_dir,"mesophotic_ord_shallowdeep.pdf"),ord_zone_composite,device=cairo_pdf,width=12,height=4,units="in")

# supplemental figure (shallow = 0-45m)
ord_zone_composite <- 
  (ord_zone_plotz$sim$depth_zone45$fish + ord_zone_plotz$sim$depth_zone45$inverts + ord_zone_plotz$sim$depth_zone45$metazoans) +#/  (ord_zone_plotz$sim$depth_zone4545$fish + ord_zone_plotz$sim$depth_zone4545$inverts + ord_zone_plotz$sim$depth_zone4545$metazoans)  +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
ord_zone_composite
ggsave(path(figure_dir,"mesophotic_ord_shallowdeep45.pdf"),ord_zone_composite,device=cairo_pdf,width=12,height=4,units="in")

#### Figure: depth zone ordinations
ord_plotz <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    map2(to_plot,names(to_plot),~{
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
    }) %>% set_names(names(to_plot))
  })
ord_composite <- 
  (ord_plotz$sim$fish + ord_plotz$sim$inverts + ord_plotz$sim$metazoans) +#/  (ord_plotz$jaccard$fish + ord_plotz$jaccard$inverts + ord_plotz$jaccard$metazoans)  +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
ord_composite

ggsave(path(figure_dir,"mesophotic_ord_samples.pdf"),ord_composite,device=cairo_pdf,width=12,height=4,units="in")

#### Figure: fish depth distributions vs eDNA detections
fish <- communities$fish$raw %>%
  inner_join(communities$fish$tax_data,by="OTU") %>%
  select(domain:species,OTU,matches(sample_pattern)) %>%
  filter(species != "dropped" & !is.na(species))

fish_species <- fish %>%
  distinct(species) %>%
  pull(species)

detected_depth <- fish %>%
  filter(species %in% fish_species) %>%
  select(family,species,matches(sample_pattern)) %>%
  pivot_longer(matches(sample_pattern),names_to="sample",values_to="reads") %>%
  inner_join(sample_data %>% select(id,depth),by=c("sample" = "id")) %>%
  filter(reads >= 10) %>%
  group_by(family,species,depth) %>%
  summarise(n=n(),reads=sum(reads)) %>%
  ungroup() %>%
  left_join(
    fb_tbl("species") %>%
      mutate(sci_name = paste(Genus, Species)) %>%
      # filter(sci_name %in% fish_species) %>%
      select(sci_name, contains('Depth')) ,
    by = c("species" = "sci_name")
  ) %>%
  select(family,species,depth,reads,fb_shallow=DepthRangeShallow,fb_deep=DepthRangeDeep) %>%
  filter(fb_deep < 150)  %>%
  ungroup() %>%
  group_by(species) %>%
  mutate(depth_diff = depth - max(fb_deep)) %>%
  ungroup() %>%
  mutate(
    depth_range = fb_deep-fb_shallow,
    species = fct_reorder(species,depth_range,.desc = TRUE)
    # species = fct_reorder(species,as.numeric(factor(family)))
  ) 

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
    map2(to_plot,names(to_plot),~{
      title <- plot_text2[.y]
      p <- plot_betadisp(.x, group="station_grouping", method=dm, list=TRUE)
      p$plot <- p$plot +
        # ggtitle(title) +
        scale_fill_manual(values=pal[c(1,4,7)],name="Site") +
        plotz_theme("light") + 
        xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
        ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
        theme(legend.position = "right")
      p$plot
    }) %>% set_names(names(to_plot))
  })

ord_site_composite <- 
  (ord_sites$sim$fish + ord_sites$sim$inverts + ord_sites$sim$metazoans) +#/  (ord_sites$jaccard$fish + ord_sites$jaccard$inverts + ord_sites$jaccard$metazoans)  +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
ord_site_composite
ggsave(path(figure_dir,"mesophotic_ord_sites.pdf"),ord_site_composite,device=cairo_pdf,width=12,height=4,units="in")


#### Supplemental Figure: species   accumulations
sup <- function(...) suppressWarnings(suppressMessages(...))
# using normalized datasets
all_models <- map2(comm_ps,names(comm_ps),~{
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
        # # distinct(sites,richness) %>%
        # mutate(depth=.y$depth_f) %>%
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
    # lomolino = try_mod(nls(richness ~ SSlomolino(sites, Asym, xmid, slope),  data=specs)),
    asymp = nls(richness ~ SSasymp(sites, Asym, R0, lrc),  data=specs),
    gompertz = nls(richness ~ SSgompertz(sites, Asym,  xmid, scal), data=specs),
    `michaelis-menten` = nls(richness ~  SSmicmen(sites, Vm, K), data=specs),
    logis = nls(richness ~ SSlogis(sites,  Asym, xmid, scal), data=specs) 
  )
  # AICs <- sort(map_dbl(modls,AIC))
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
  
  # slopes <- diff(predictions$y)/diff(predictions$x)
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
      # ggtitle(tit) +
      xlab("Replicates") + 
      ylab("zOTUs") + 
      theme_bw() #+
      # theme(panel.grid = element_blank())
  })
accum_composite <- 
  accum_plotz$fish + accum_plotz$inverts + accum_plotz$metazoans + 
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
accum_composite
ggsave(path(figure_dir,"mesophotic_accum.pdf"),accum_composite,device=cairo_pdf,width=12,height=4,units="in")



# manuscript tables -------------------------------------------------------
table_dir <- "~/projects/dissertation/manuscript/tables"

write_ms_table <- function(tbl,file,caption="",...) {
  if (caption != "") {
    write_lines(str_c(": ",caption),file)
    write_csv(tbl,file,append=TRUE,col_names = TRUE,...)
  } else {
    write_csv(tbl,file,...)    
  }
}

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
        value = "Raw sequence reads"
      )
    s <- all_samples %>%
      left_join(counts,by="Sample ID") #%>%
      # bind_rows(counts %>% filter(str_detect(`Sample ID`,blank_pattern)))

    rs <- s %>% 
      mutate(
        "Site" = station_grouping,
        "Depth" = depth_f,
        "Assay" = plot_text[.y],
        "Amplification" = !is.na(`Raw sequence reads`),
        "Type" = case_when(
          str_detect(substrate,"control") ~ "field control",
          str_detect(`Sample ID`,sample_pattern) ~ "eDNA Sample",
          str_detect(`Sample ID`,blank_pattern) ~ "extraction blank",
          TRUE ~ ""
        ),
        "Passing" = !is.na(`Raw sequence reads`) & `Raw sequence reads` > min_total
      ) %>%
      select(`Sample ID`,Assay,Site,Depth,`Raw sequence reads`,Type,Passing,Amplification) %>%
      arrange(`Sample ID`)
  }) %>%
  left_join(normed,by=c('Sample ID','Assay'))

rs <- reads_summary %>%
  replace_na(list(`Raw sequence reads`=0,`Normalized sequence reads`=0)) %>%
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
    `Raw sequence reads (16S; 18S; COI)`=str_c(scales::comma(na.omit(`Raw sequence reads`),accuracy = 1),collapse="; "),
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
  # arrange(Site,Depth,`Sample ID`)
  arrange(order,`Sample ID`) %>%
  select(-order)

write_ms_table(rs,path(table_dir,"mesophotic_reads_summary.csv"),"Sequencing results for mesophotic eDNA samples across assay types",na="")

# numbers reported in the text --------------------------------------------

#### eDNA read counts for different categories squished together
cat("sample:\n")
walk2(comm_ps,names(comm_ps),~{
  # f <- merge_samples(.x,"station")
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
