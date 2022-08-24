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

# create the metazoan subset of our three communities
minimals <- 1000                                      # minimum reads to retain a sample
animals <- comm_ps %>%
  map(~{
    f <- subset_taxa(.x,kingdom == "Metazoa")         # filter by metazoans
    f <- prune_samples(sample_sums(f) >= minimals,f)  # filter by reads >= minimum
    new_otu <- otu_table(f) %>%                       # re-rarefy new otu table to equal depth
      as("matrix") %>%
      rrarefy.perm(n=100)
    otu_table(f) <- otu_table(new_otu,taxa_are_rows = FALSE) # reassign new rarefied otu table
    return(f)
  })
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
  'metazoans' = 'COI-Metazoan'
)

# set up depth zone color palette
pal <- rev(beyonce_palette(75))[5:12]
pal <- pal[c(1:3,5,4,6:7)]
pal[4] <- '#827222'


# start here --------------------------------------------------------------


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
        adonis_shallowdeep = adonis(dd ~ depth_zone, data=sd, permutations=perm),
        adonis_shallowdeep45 = adonis(dd ~ depth_zone45, data=sd, permutations=perm),
        # adonis_sites = adonis(dd ~ depth_f*station_grouping, data=sd, permutations = perm), 
        # adonis_shallowdeep_sites = adonis(dd ~ depth_zone*station_grouping, data=sd, permutations=perm),
        adonis_depth = adonis(dd ~ depth_f, data=sd, permutations = perm), 
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



# beta diversity calculations ---------------------------------------------
beta_diversity <- map(comm_ps,~{
  # .x <- merge_samples(.x,"depth_f")
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
beta_diversity_animals <- map(animals,~{
  # .x <- merge_samples(.x,"depth_f")
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


# manuscript figure plots -------------------------------------------------
#### Figure: cluster plots
merge_factor <- "depth_f"
cluster_plotz <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    comm_ps %>%
      map2(names(.),~{
        # title <- switch(
        #   dm,
        #   sim = bquote(.(plot_text2[.y]) ~ "(" * beta[sim] * ")" ),
        #   jaccard = bquote(.(plot_text2[.y]) ~ "(" * beta[j] * ")" )
        # )
        title <- plot_text2[.y]
        merged <- merge_samples(.x,merge_factor)
        comm_merged <- distance(merged,method=dm)
        # labels(comm_merged) <- labels(comm_merged) %>%
        #   map_chr(~str_c(.x," "))
        labels(comm_merged) <- str_c(labels(comm_merged)," ")
        ggplot(as.dendrogram(hclust(comm_merged),hang=0.5)) +
          theme(
            plot.caption = element_text(hjust=0.5,size=14)
          ) + 
          # labs(caption=title) +
          expand_limits(y=-0.25)
          # ggtitle(title)
      }) 
  })

cluster_composite <- (cluster_plotz$sim$fish + cluster_plotz$sim$inverts + cluster_plotz$sim$metazoans) +
  plot_annotation(tag_levels="A")
cluster_composite
ggsave("~/projects/dissertation/manuscript/figures/mesophotic_cluster.pdf",cluster_composite,device=cairo_pdf,width=12,height=4,units="in")

#### Figure: cluster plots (animals)
merge_factor <- "depth_f"
aminal_cluster_plotz <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    animals %>%
      map2(names(.),~{
        # title <- switch(
        #   dm,
        #   sim = bquote(.(plot_text2[.y]) ~ "(" * beta[sim] * ")" ),
        #   jaccard = bquote(.(plot_text2[.y]) ~ "(" * beta[j] * ")" )
        # )
        title <- plot_text2[.y]
        merged <- merge_samples(.x,merge_factor)
        comm_merged <- distance(merged,method=dm)
        # labels(comm_merged) <- labels(comm_merged) %>%
        #   map_chr(~str_c(.x," "))
        labels(comm_merged) <- str_c(labels(comm_merged)," ")
        ggplot(as.dendrogram(hclust(comm_merged),hang=0.5)) +
          theme(
            plot.caption = element_text(hjust=0.5,size=14)
          ) + 
          # labs(caption=title) +
          expand_limits(y=-0.25)
          # ggtitle(title)
      }) 
  })

aminal_cluster_composite <- (aminal_cluster_plotz$sim$inverts + aminal_cluster_plotz$sim$metazoans) +
  plot_annotation(tag_levels="A")
aminal_cluster_composite
ggsave("~/projects/dissertation/manuscript/figures/mesophotic_cluster_animals.pdf",aminal_cluster_composite,device=cairo_pdf,width=8,height=4,units="in")

#### Figure: shallow vs deep ordinations
zone_groupings <- c("fish" = "depth_zone", "inverts" = "depth_zone45", "metazoans" = "depth_zone45")
beta_zone_plotz <- c("depth_zone","depth_zone45") %>%
  set_names() %>%
  map(~{
    dm <- "sim"
    zone <- .x
    map2(to_plot,names(to_plot),~{
      # title <- switch(
      #   dm,
      #   sim = bquote(.(plot_text2[.y]) ~ "(" * beta[sim] * ")" ),
      #   jaccard = bquote(.(plot_text2[.y]) ~ "(" * beta[j] * ")" ),
      #   # jaccard = bquote(beta[j] ~  .(plot_text2[.y]))
      # )
      title <- plot_text2[.y]
      # p <- plot_betadisp(.x, group="depth_zone", method=dm, list=TRUE)
      p <- plot_betadisp(.x, group=zone, method=dm, list=TRUE)
      p$plot <- p$plot +
        # ggtitle(title) +
        scale_fill_manual(values=pal[c(1,length(pal))],name="Depth Zone") + #,labels=c("Shallow","Deep")) +
        plotz_theme("light") +
        xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
        ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
        theme(legend.position = "right")
      p$plot
    }) %>% set_names(names(to_plot))
  })

# main text figure (shallow = 0-30m)
beta_zone_composite <- 
  (beta_zone_plotz$depth_zone$fish + beta_zone_plotz$depth_zone$inverts + beta_zone_plotz$depth_zone$metazoans) +#/  (beta_zone_plotz$depth_zone45$fish + beta_zone_plotz$depth_zone45$inverts + beta_zone_plotz$depth_zone45$metazoans)  +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
beta_zone_composite
ggsave("~/projects/dissertation/manuscript/figures/mesophotic_beta_shallowdeep.pdf",beta_zone_composite,device=cairo_pdf,width=12,height=4,units="in")

# supplemental figure (shallow = 0-45m)
beta_zone_composite <- 
  (beta_zone_plotz$depth_zone45$fish + beta_zone_plotz$depth_zone45$inverts + beta_zone_plotz$depth_zone45$metazoans) +#/  (beta_zone_plotz$depth_zone4545$fish + beta_zone_plotz$depth_zone4545$inverts + beta_zone_plotz$depth_zone4545$metazoans)  +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
beta_zone_composite
ggsave("~/projects/dissertation/manuscript/figures/mesophotic_beta_shallowdeep45.pdf",beta_zone_composite,device=cairo_pdf,width=12,height=4,units="in")

#### Figure: depth zone ordinations
beta_plotz <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    map2(to_plot,names(to_plot),~{
      # title <- switch(
      #   dm,
      #   sim = plot_text2[.y]bquote(.(plot_text2[.y]) ~ "(" * beta[sim] * ")" ),
      #   jaccard = bquote(.(plot_text2[.y]) ~ "(" * beta[j] * ")" )
      # )
      title <- plot_text2[.y]
      p <- plot_betadisp(.x, group="depth_f", method=dm, list=TRUE, expand=TRUE)
      p$plot <- p$plot +
        # ggtitle(title) +
        # geom_point(aes(x=x,y=y,color=depth_f,shape=station_grouping),size=4,alpha=0.7) +
        scale_fill_manual(values=pal,name="Depth") +
        # scale_color_manual(values=alpha(pal,0.5),name="Depth") + 
        # scale_shape(name="Sampling Site") + 
        plotz_theme("light") +
        xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
        ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
        theme(legend.position = "right") +
        guides(color="none")
      p$plot
    }) %>% set_names(names(to_plot))
  })
# leg <- cowplot::ggdraw(cowplot::get_legend(beta_plotz$sim$fish + theme(legend.position="right")))
beta_composite <- 
  (beta_plotz$sim$fish + beta_plotz$sim$inverts + beta_plotz$sim$metazoans) +#/  (beta_plotz$jaccard$fish + beta_plotz$jaccard$inverts + beta_plotz$jaccard$metazoans)  +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
beta_composite

ggsave("~/projects/dissertation/manuscript/figures/mesophotic_beta_samples.pdf",beta_composite,device=cairo_pdf,width=12,height=4,units="in")

#### Figure: fish depth distributions
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
  mutate(
    depth_range = fb_deep-fb_shallow,
    species = fct_reorder(species,depth_range,.desc = TRUE)
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
  # geom_errorbar(aes(y=species,xmin=fb_shallow,xmax=fb_deep),width=0.5,color="dodgerblue4") +
  geom_errorbar(aes(x=species,ymin=fb_shallow,ymax=fb_deep),width=0.5,color="dodgerblue4") +
  geom_point(aes(x=species,y=depth,size=reads,color=recode_factor(factor(str_c(depth,'m')), "0m" = "Surface")),alpha=0.7) +
  scale_color_manual(values=pal,name="Depth zone") +
  scale_size(name="Rarefied read count",labels=scales::comma,range=c(2,12),breaks=c(1e2,1e3,1e4,5e4,1e5)) +
  # scale_y_discrete(limits=rev) +
  # scale_y_continuous(breaks=seq(0,200,by=50),minor_breaks=seq(0,200,by=10)) +
  scale_y_reverse(breaks=seq(0,300,by=30),minor_breaks=seq(0,300,by=15)) +
  scale_x_discrete(position="top") + 
  plotz_theme("light") +
  theme(
    axis.text.x = element_text(face="italic",angle=30,vjust=1,hjust=0),
    legend.position = "right",
    legend.key = element_rect(color=NA),
    panel.grid.major = element_line(color="grey95"),
    panel.grid.minor = element_line(color="grey97"),
    panel.grid.major.x = element_blank()
    # legend.background = element_rect(fill="grey95")
    # panel.grid =  
    # legend.key.size = unit(3,"lines") 
  ) +
  guides(color = guide_legend(override.aes = list(size=7))) +
  ylab("Detection depth (m)") + 
  xlab("Species") 
depth_plotz
ggsave("~/projects/dissertation/manuscript/figures/mesophotic_fish_depth.pdf",depth_plotz,device=cairo_pdf,width=12,height=10,units="in")

#### Supplemental Figure: ordinations by site
beta_sites <- distance_methods %>%
  set_names() %>%
  map(~{
    dm <- .x
    map2(to_plot,names(to_plot),~{
      # title <- switch(
      #   dm,
      #   sim = bquote(.(plot_text2[.y]) ~ "(" * beta[sim] * ")" ),
      #   jaccard = bquote(.(plot_text2[.y]) ~ "(" * beta[j] * ")" ),
      #   # jaccard = bquote(beta[j] ~  .(plot_text2[.y]))
      # )
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
# beta_site_plot <- (beta_sites$sim$fish + beta_sites$sim$inverts + beta_sites$sim$metazoans) / (beta_sites$jaccard$fish + beta_sites$jaccard$inverts + beta_sites$jaccard$metazoans)
beta_site_composite <- 
  (beta_sites$sim$fish + beta_sites$sim$inverts + beta_sites$sim$metazoans) +#/  (beta_sites$jaccard$fish + beta_sites$jaccard$inverts + beta_sites$jaccard$metazoans)  +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
beta_site_composite
ggsave("~/projects/dissertation/manuscript/figures/mesophotic_beta_sites.pdf",beta_site_composite,device=cairo_pdf,width=12,height=4,units="in")


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
      geom_errorbar(data=d,aes(x=sites,ymin=richness-sd,ymax=richness+sd),width=0.2) + 
      geom_line(data=p,aes(x=x,y=y),color="blue",size=1.1) + 
      geom_vline(xintercept=.x$asymptote,color="red") + 
      scale_color_manual(values=pal,name="Depth Zone") + 
      xlim(1,25) +
      # ggtitle(tit) +
      xlab("Replicates") + 
      ylab("zOTUs") + 
      theme_bw()
  })
accum_composite <- 
  accum_plotz$fish + accum_plotz$inverts + accum_plotz$metazoans + 
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
accum_composite
ggsave("~/projects/dissertation/manuscript/figures/mesophotic_accum.pdf",accum_composite,device=cairo_pdf,width=12,height=4,units="in")

