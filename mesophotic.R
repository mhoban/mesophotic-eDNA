library(here)
library(EcolUtils)
library(beyonce)
library(plotly)
library(indicspecies)
library(betapart)
library(Polychrome)
library(patchwork)


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

summary_grouping <- c("station_grouping","station")  # how to group samples for replication report

max_blank <- 5                    # maximum reads in blank to believe it
min_total <- 2e4                   # minimum total reads in sample

insect_classify <- FALSE
include_unassigned <- TRUE        # include zOTUs that didn't blast to anything
filter_unknown <- FALSE

distance_methods <- c("sim","jaccard")
distance_method <- "sim"
plot_theme <- "dark"

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
    depth_zone = cut(depth,c(-Inf,30,Inf),labels=c("shallow","deep")),
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
      column_to_rownames("site") %>%
      rrarefy.perm(n=100)
    ps <- as_ps2(rarefied,.x$tax_data,sample_data)
    pruno <- filter_taxa(ps,function(x) {
      sum(x[negative_controls],na.rm=T) == 0 & sum(x) > 0
    })
    ps <- prune_taxa(pruno,ps)
    ps <- subset_taxa(ps, !family %in% c("Hominidae","Bovidae","Felidae","Salmonidae"))
    
    prune_samples(sample_sums(ps) > 0, ps)
  }) 

# text map for plots
plot_text <- c(
  'fish' = 'Fishes (16S rRNA)',
  'inverts' = '"Eukaryotes" (18S rRNA)',
  'metazoans' = '"Metazoans" (COI)'
)
plot_text2 <- c(
  'fish' = '16S-Fish',
  'inverts' = '18S-Eukaryote',
  'metazoans' = 'COI-Metazoan'
)

# set up depth zone color palette
pal <- rev(beyonce_palette(75))[5:12]
pal <- pal[c(1:3,5,4,6:7)]
pal[4] <- '#827222'


# start here --------------------------------------------------------------



# depth zones -------------------------------------------------------------
# TODO: plot ideas
# a bar plot showing average jaccard dissimilarity for each depth zone pair
# with error bars, I guess

# this plots average jaccard dissimilarity by depth comparison
# is this useful? I don't know!

diff_plotz <- map2(comm_ps,names(comm_ps),~{
  sd <- sample_tibble(.x)
  ddm <- as.matrix(distance(.x,method=distance_method))
  ddmd <- as_tibble(ddm,rownames="left") %>%
    pivot_longer(-left,names_to="right",values_to="jaccard") %>%
    filter(left != right) %>%
    filter(!duplicated(paste0(pmax(left, right), pmin(left, right)))) %>%
    left_join(sd %>% select(sample,depth_left=depth), by=c("left" = "sample")) %>%
    left_join(sd %>% select(sample,depth_right=depth), by=c("right" = "sample")) %>%
    filter(depth_left != depth_right) %>%
    mutate(comparison = str_c(pmin(depth_left,depth_right),"-",pmax(depth_left,depth_right))) %>% 
    select(comparison,jaccard) %>%
    group_by(comparison) %>%
    summarise(sd=sd(jaccard),jaccard=mean(jaccard),min=jaccard-sd,max=jaccard+sd)
  ggplot(ddmd,aes(x=comparison,y=jaccard)) + 
    geom_col(position="dodge") + 
    geom_errorbar(aes(ymin=min,ymax=max),width=0.3) + 
    ggtitle(.y)
})

# PCoA ordinations --------------------------------------------------------

### PCoA ordination by depth zone and station grouping
ord_plotz <- map2(comm_ps,names(comm_ps),~{
  ord <- ordinate(.x,"PCoA",distance=distance_method)
  ord_plot <- plot_ordination(.x,ord,type="samples",color="depth_f",shape="station_grouping") +
    geom_point(size=8,aes(text=station)) + 
    ggtitle(str_glue("PCoA ordination: {plot_text[.y]}")) + 
    # scale_color_manual(values=rev(beyonce_palette(75))[5:12], name="Depth Zone") + 
    scale_color_manual(values=pal, name="Depth Zone") +
    scale_shape(name="Sampling Site") + 
    theme_bw() + 
    theme(
      plot.background = element_rect(fill="#1a1a1a"),
      text = element_text(color="grey80"),
      panel.background = element_rect(fill="#000000"),
      panel.grid = element_blank(),
      legend.background = element_rect(fill="grey27"),
      legend.key = element_rect(fill="grey76", color = "black"),
      legend.key.size = unit(1.4,"line")
    )
  
  ord_plot + 
    guides(shape = guide_legend(override.aes = list(size=5)))
})

library(fs)
walk2(ord_plotz,names(ord_plotz),~{
  dir <- "~/projects/dissertation/admin/committee meetings/february 24, 2022/img"
  ggsave(.x, file=path(dir,str_glue("pcoa_{.y}.svg")), device="svg", width = 16, height = 9, units="in")
  # ggsave(.x, file=str_glue("~/school/presentations/icrs/2021/img/pcoa_{.y}.svg"), device="svg", width = 16, height = 9, units="in")
})


# beta dispersion plots (it's a PCoA with polygons) --------------------------------------------------------
beta_plotz <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    map2(comm_ps,names(comm_ps),~{
      title <- switch(
        dm,
        sim = bquote(.(plot_text2[.y]) ~ "(" * beta[sim] * ")" ),
        jaccard = bquote(.(plot_text2[.y]) ~ "(" * beta[j] * ")" ),
        # jaccard = bquote(beta[j] ~  .(plot_text2[.y]))
      )
      p <- plot_betadisp(.x, group="depth_f", method=dm, list=TRUE)
      p$plot <- p$plot +
        # ggtitle(plot_text[.y]) +
        # ggtitle(expression(beta[j]~tau)) + 
        ggtitle(title) +
        scale_fill_manual(values=pal,name="Depth") +
        plotz_theme("light") +
        xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
        ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
        theme(legend.position = "right")
      p$plot
    }) %>% set_names(names(comm_ps))
  })
# leg <- cowplot::ggdraw(cowplot::get_legend(beta_plotz$sim$fish + theme(legend.position="right")))
beta_composite <- 
  (beta_plotz$sim$fish + beta_plotz$sim$inverts + beta_plotz$sim$metazoans) /  (beta_plotz$jaccard$fish + beta_plotz$jaccard$inverts + beta_plotz$jaccard$metazoans)  +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
beta_composite

ggsave("~/projects/dissertation/manuscript/figures/mesophotic_beta_samples.pdf",beta_composite,device=cairo_pdf,width=12,height=9,units="in")

# beta dispersion plots by site -------------------------------------------

beta_sites <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    map2(comm_ps,names(comm_ps),~{
      title <- switch(
        dm,
        sim = bquote(.(plot_text2[.y]) ~ "(" * beta[sim] * ")" ),
        jaccard = bquote(.(plot_text2[.y]) ~ "(" * beta[j] * ")" ),
        # jaccard = bquote(beta[j] ~  .(plot_text2[.y]))
      )
      p <- plot_betadisp(.x, group="station_grouping", method=dm, list=TRUE)
      p$plot <- p$plot +
        # ggtitle(str_glue("Beta dispersion plot by depth ({plot_text[.y]}")) +
        ggtitle(title) +
        scale_fill_manual(values=pal[c(1,4,7)],name="Site") +
        plotz_theme("light") + 
        xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
        ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
        theme(legend.position = "right")
      p$plot
    }) %>% set_names(names(comm_ps))
  })
# beta_site_plot <- (beta_sites$sim$fish + beta_sites$sim$inverts + beta_sites$sim$metazoans) / (beta_sites$jaccard$fish + beta_sites$jaccard$inverts + beta_sites$jaccard$metazoans)
beta_site_composite <- 
  (beta_sites$sim$fish + beta_sites$sim$inverts + beta_sites$sim$metazoans) /  (beta_sites$jaccard$fish + beta_sites$jaccard$inverts + beta_sites$jaccard$metazoans)  +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
beta_site_composite
ggsave("~/projects/dissertation/manuscript/figures/mesophotic_beta_sites.pdf",beta_site_composite,device=cairo_pdf,width=12,height=9,units="in")


# heatmap plots -----------------------------------------------------------
heatmap_plotz <- map2(comm_ps,names(comm_ps), ~{
  sd <- sample_tibble(.x)
  # td <- taxa_tibble(.x)  
  
  sample_order <- sd %>%
    arrange(depth,station_grouping,sample) %>%
    pull(sample)
  
  # taxa_order <- td %>%
  #   arrange(domain,phylum,class,order,family,genus,species) %>%
  #   pull(zotu)
  
  p <- viridis::viridis_pal(1,0,1)
  pp <- p(10)
  
  hm <- plot_heatmap(
    .x,
    method='PCoA',
    # formula = ~depth,
    # distance='jaccard',
    distance=distance_method,
    sample.label = 'depth',
    taxa.label = 'genus',
    sample.order = sample_order,
    low = pp[1],
    high = pp[10],
    max.label = 50
  ) + 
    theme(
      plot.background = element_rect(fill="#1a1a1a"),
      text = element_text(color="grey80"),
      panel.background = element_rect(fill="#000000"),
      panel.grid = element_blank(),
      legend.background = element_rect(fill="grey27"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) + 
    ggtitle(str_glue("OTU PCoA heatmap by depth ({plot_text[.y]})"))
  # p$scales$scales[[1]]$name <- "My X-Axis"
  hm$scales$scales[[1]]$name <- "Taxa (genus)"
  hm$scales$scales[[2]]$name <- "Depth Zone (m)"
  hm <- hm + 
    labs(fill = "Normalized\nsequence\nreads")
  hm
}) %>% set_names(names(comm_ps))

walk2(heatmap_plotz,names(heatmap_plotz),~{
  dir <- "~/projects/dissertation/admin/committee meetings/february 24, 2022/img"
  ggsave(.x, file=path(dir,str_glue("heatmap_{.y}.svg")), device="svg", width = 16, height = 9, units="in")
  # ggsave(.x, file=str_glue("~/school/presentations/icrs/2021/img/heatmap_{.y}.svg"), device="svg", width = 16, height = 9, units="in")
})


# Bar plots by depth zone -------------------------------------------------

fills <- list(
  fish = "family",
  inverts = "class",
  metazoans = "class"
)

### Bar plots by various taxa and depth zone
bar_plotz <- map2(comm_ps,names(comm_ps),~{
  
  fill <- fills[[.y]]
  
  .x <- merge_samples(.x,"depth_f")
  sd <- sample_data(.x)
  sd$depth_f <- factor(sample_names(.x))
  sd$depth_f <- fct_reorder(sd$depth_f,sd$depth)
  sample_data(.x) <- sd
  otu_table(.x) <- otu_table(wisconsin(otu_table(.x)),taxa_are_rows = FALSE)
  # sample_data(.x)$depth_f <- factor(sample_names(.x),levels=c("Surface","15m","30m","45m","60m","75m","90m"))
  # .x <- subset_taxa(.x,kingdom == "Metazoa" & family != "Hominidae")
  # otu_table(.x) <- otu_table(wisconsin(otu_table(.x)),taxa_are_rows = FALSE)
  
  pdata <- as_tibble(psmelt(.x)) %>%
    # filter(
    #   kingdom == "Metazoa",
    #   family != "Hominidae"
    # ) %>%
    arrange(domain,kingdom,phylum,order,class,family,genus,species) %>%
    mutate(
      across(
        domain:species,
        ~case_when(
          is.na(.x) ~ "Unknown",
          TRUE ~ .x
        )
      ),
      across(
        domain:species,
        ~factor(.x,levels=unique(.x))
      )
    ) %>%
    group_by(depth_f,across(all_of(fill))) %>%
    summarise(Abundance=sum(Abundance)) %>%
    ungroup()
  
  x <- "depth_f"
  y <- "Abundance"
  
  levs <- pdata %>%
    pull(fill) %>%
    unique() %>%
    length()
  
  clr <- unname(Polychrome::createPalette(levs,seedcolors = pal[c(1,4,7)]))
  
  p <- ggplot(pdata, aes_string(x = x, y = y, fill = fill)) +
    geom_col(position = "stack", color = "black") #+ 
    # facet_wrap(~station_grouping)
  p <- p + scale_fill_manual(values=clr, name=fill)
  # p <- p + scale_fill_discrete(name="family") 
  p + 
    ggtitle(str_c(str_glue("Relative abundance by depth zone ({fill}) - "),plot_text[.y])) +
    xlab("Depth Zone") + 
    ylab("Sequence Reads") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = -90, hjust = 0, color="grey65"),
      axis.text.y = element_text(color="grey65"),
      plot.background = element_rect(fill="#1a1a1a"),
      text = element_text(color="grey80"),
      panel.background = element_rect(fill="#000000"),
      panel.grid = element_blank(),
      legend.background = element_rect(fill="grey27"),
      legend.key = element_rect(fill="grey76", color = NA)
    )
}) %>% set_names(names(comm_ps))

walk2(bar_plotz,names(bar_plotz),~{
  dir <- "~/projects/dissertation/admin/committee meetings/february 24, 2022/img"
  ggsave(.x, file=path(dir,str_glue("stackedbar_{.y}.svg")), device="svg", width = 16, height = 9, units="in")
  # ggsave(.x, file=str_glue("~/school/presentations/icrs/2021/img/stackedbar_{.y}.svg"), device="svg", width = 16, height = 9, units="in")
})

# permanova analyses ------------------------------------------------------
anovas <- map2(comm_ps,names(comm_ps),~{
  perm <- 8000
  sd <- sample_tibble(.x)
  dd <-distance(.x,method=distance_method)
  list(
    adonis_shallowdeep = adonis(dd ~ depth_zone, data=sd, permutations=perm),
    adonis_shallowdeep_sites = adonis(dd ~ depth_zone*station_grouping, data=sd, permutations=perm),
    adonis = adonis(dd ~ depth_f, data=sd, permutations = perm), 
    adonis_sites = adonis(dd ~ depth_f*station_grouping, data=sd, permutations = perm), 
    anosim_shallowdeep = anosim(dd,sd$depth_zone, permutations = perm), 
    anosim = anosim(dd, sd$depth_f, permutations = perm),
    # indicators_shallowdeep = multipatt(community, sd$depth_zone, control = how(nperm=perm)),
    # indicators = multipatt(community, sd$depth_f, control = how(nperm=perm)),
    # simper = simper(community, sd$depth_f, permutations = perm),
    # pairwise = adonis.pair(dd, sd$depth_f, nper = perm, corr.method = "bonferroni")
    pairwise = pairwise_adonis(dd, sd$depth_f, permutations = perm,correction = 'BH')
  )
}) %>% set_names(names(comm_ps))


write_rds(simpson_anovas,here::here("output","simpson_anovas.rds"))
write_rds(jaccard_anovas,here::here("output","jaccard_anovas.rds"))


# cluster plots -----------------------------------------------------------
merge_factor <- "depth_f"
cluster_plotz <- comm_ps %>%
  map2(names(.),~{
    comm_dist <- distance(.x,method=distance_method)
    merged <- merge_samples(.x,merge_factor)
    comm_merged <- distance(merged,method=distance_method)
    
    list(
      samples = hclust(comm_dist),
      merged = hclust(comm_merged)
    )
  })

# beta diversity calculations ---------------------------------------------

beta_diversity <- map(comm_ps,~{
  comm <- as(otu_table(.x),"matrix")
  # do the abundance-based beta partitioning
  bc <- betapart.core.abund(comm)
  comm[comm > 0] <- 1
  # now the presence-absence one
  bp <- betapart.core(comm)
  
   
  bcm <- beta.multi.abund(bc,"bray")
  bpm <- beta.multi(bp)
  
  list(abund = bcm, regular = bpm) 
}) %>% set_names(names(comm_ps))
