library(here)
library(EcolUtils)
library(beyonce)
library(plotly)
library(indicspecies)
library(betapart)


# Community setup and initial analysis ------------------------------------

# shorten calls to suppressMessages (readr outputs all sorts of nonsense)
sm <- suppressMessages            

# set the project, directory, and markers for this analysis
project <- "mesophotic"
# project <- "marianas"
# project <- "seagrant"
project_dir <- here::here("working","flow-output",project)
sample_pattern <- "^BMSP[0-9]+$"  # pattern to match sample IDs
# sample_pattern <- "^SGP[0-9]+$"  # pattern to match sample IDs
# sample_pattern <- "^NMI[0-9]+$"  # pattern to match sample IDs
# blank_pattern <- "^Blank-[0-9]+$"      # pattern to match extraction blank IDs
blank_pattern <- "^B[0-9]+$"      # pattern to match extraction blank IDs
abundance_threshold <- 0.001     # minimum threshold for relative abundance

max_blank <- 5                    # maximum reads in blank to believe it
min_total <- 12                   # minimum total reads in sample

# load our OTU tables and sample data
source("analyze.R")

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
    rarefied <- .x$relative %>%
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
  }) %>%
  set_names(c("fish","inverts","metazoans"))

descriptions <- c(
  'fish' = 'Fishes (16S rRNA)',
  'inverts' = '"Eukaryotes" (18S rRNA)',
  'metazoans' = '"Metazoans" (COI)'
)

# set up depth zone color palette
pal <- rev(beyonce_palette(75))[5:12]
pal <- pal[c(1:3,5,4,6:7)]
pal[4] <- '#827222'

# PCoA ordinations --------------------------------------------------------

### PCoA ordination by depth zone and station grouping
ord_plotz <- map2(comm_ps,names(comm_ps),~{
  # .x <- subset_taxa(.x,domain == 'Animals')
  # # animals <- tax_glom(animals, taxrank = "genus")
  # .x <- subset_taxa(.x, !family %in% c("Hominidae","Bovidae","Felidae","Salmonidae"))
  # .x <- prune_samples(sample_sums(.x) > 0, .x)
  ord <- ordinate(.x,"PCoA","jaccard")
  ord_plot <- plot_ordination(.x,ord,type="samples",color="depth_f",shape="station_grouping") +
    geom_point(size=8,aes(text=station)) + 
    ggtitle(str_glue("PCoA ordination: {descriptions[.y]}")) + 
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
}) %>% set_names(names(comm_ps))

walk2(ord_plotz,names(ord_plotz),~{
  ggsave(.x, file=str_glue("~/school/presentations/icrs/2021/img/pcoa_{.y}.svg"), device="svg", width = 16, height = 9, units="in")
})


# Bar plots by depth zone -------------------------------------------------

### Bar plots by family and depth zone
bar_plotz <- map2(comm_ps,names(comm_ps),~{
  
  pdata <- as_tibble(psmelt(.x)) %>%
    filter(
      domain == "Animals",
      family != "Hominidae"
    ) %>%
    arrange(domain,phylum,order,class,family,genus,species) %>%
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
    group_by(depth_f,family) %>%
    summarise(Abundance=sum(Abundance))
  
  x <- "depth_f"
  fill <- "family"
  y <- "Abundance"
  
  levs <- pdata %>%
    pull(fill) %>%
    levels() %>%
    length()
  
  clr <- unname(Polychrome::createPalette(levs,seedcolors = pal[c(1,4,7)]))
  
  p <- ggplot(pdata, aes_string(x = x, y = y, fill = fill)) +
    geom_col(position = "fill", color = "black")
  p <- p + scale_fill_manual(values=clr, name="family")
  # p <- p + scale_fill_discrete(name="family") 
  p + 
    ggtitle(str_c("Relative abundance by depth zone (families) - ",descriptions[.y])) +
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
  ggsave(.x, file=str_glue("~/school/presentations/icrs/2021/img/stackedbar_{.y}.svg"), device="svg", width = 16, height = 9, units="in")
})


# heatmap plots -----------------------------------------------------------

heatmap_plotz <- map2(comm_ps,names(comm_ps), ~{
  animals <- .x  
  # animals <- subset_taxa(animals,domain == 'Animals')
  # animals <- tax_glom(animals, taxrank = "genus")
  animals <- subset_taxa(animals, !family %in% c("Hominidae","Bovidae","Felidae","Salmonidae"))
  animals <- prune_samples(sample_sums(animals) > 0, animals)
  
  sd <- sample_tibble(animals)
  td <- taxa_tibble(animals)  
  
  sample_order <- sd %>%
    arrange(depth,station_grouping,sample) %>%
    pull(sample)
  
  taxa_order <- td %>%
    arrange(domain,phylum,class,order,family,genus,species) %>%
    pull(zotu)
  
  p <- viridis::viridis_pal(1,0,1)
  pp <- p(10)
  
  hm <- plot_heatmap(
    animals,
    method='PCoA',
    # formula = ~depth,
    distance='jaccard',
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
    ggtitle(str_glue("OTU PCoA heatmap by depth ({descriptions[.y]})"))
  # p$scales$scales[[1]]$name <- "My X-Axis"
  hm$scales$scales[[1]]$name <- "Taxa (genus)"
  hm$scales$scales[[2]]$name <- "Depth Zone (m)"
  hm <- hm + 
    labs(fill = "Normalized\nsequence\nreads")
  hm
}) %>% set_names(names(comm_ps))

walk2(heatmap_plotz,names(heatmap_plotz),~{
  ggsave(.x, file=str_glue("~/school/presentations/icrs/2021/img/heatmap_{.y}.svg"), device="svg", width = 16, height = 9, units="in")
})


# permanova analyses ------------------------------------------------------

animals_only <- FALSE

anovas <- map2(comm_ps,names(comm_ps),~{
  perm <- 8000
  if (animals_only) {
    .x <- subset_taxa(.x,domain == 'Animals')
    .x <- prune_samples(sample_sums(.x) > 0, .x)
  }
  community <- otu_table(.x)
  sd <- sample_tibble(.x)
  dd <- vegdist(community, method = "jaccard")
  list(
    adonis_shallowdeep = adonis(dd ~ depth_zone, data=sd, permutations=perm),
    adonis = adonis(dd ~ depth_f, data=sd, permutations = perm), 
    anosim_shallowdeep = anosim(dd,sd$depth_zone, permutations = perm),
    anosim = anosim(dd, sd$depth_f, permutations = perm),
    # indicators_shallowdeep = multipatt(community, sd$depth_zone, control = how(nperm=perm)),
    # indicators = multipatt(community, sd$depth_f, control = how(nperm=perm)),
    # simper = simper(community, sd$depth_f, permutations = perm),
    # pairwise = adonis.pair(dd, sd$depth_f, nper = perm, corr.method = "bonferroni")
    pairwise = pairwise_adonis(dd, sd$depth_f, permutations = perm,correction = 'BH')
  )
}) %>% set_names(names(comm_ps))

write_rds(anovas,here::here("data","anovas.rds"))


# betadisper plots --------------------------------------------------------

beta_plotz <- map2(comm_ps,names(comm_ps),~{
  p <- plot_betadisp(.x, 'depth_f', 'jaccard', list=TRUE)
  p$plot <- p$plot +
    ggtitle(str_glue("Beta dispersion plot by depth ({descriptions[.y]}")) +
    scale_fill_manual(values=pal,name="Depth") +
    theme_bw() + 
    theme(
      plot.background = element_rect(fill="#1a1a1a"),
      text = element_text(color="grey80"),
      panel.background = element_rect(fill="#000000"),
      panel.grid = element_blank(),
      legend.background = element_rect(fill="grey27"),
      legend.key = element_rect(fill="grey76", color = "black"),
      legend.key.size = unit(1.4,"line")
    ) +
    xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
    ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})"))
  p
}) %>% set_names(names(comm_ps))

walk2(beta_plotz,names(beta_plotz),~{
  ggsave(.x$plot, file=str_glue("~/school/presentations/icrs/2021/img/beta_{.y}.svg"), device="svg", width = 16, height = 9, units="in")
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
  bpm <- beta.multi(bp, "jaccard")
  
  list(abund = bcm, regular = bpm) 
}) %>% set_names(names(comm_ps))
