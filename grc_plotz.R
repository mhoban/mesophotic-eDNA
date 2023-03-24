# cluster plots -----------------------------------------------------------
figure_dir <- "/Users/deadbilly/Documents/school/presentations/grc mesophotic/GRC/img/figures"
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

beta_pairs_composite <- beta_pairs %>%
  imap(~{
    dataset_name <- .y 
    label_map = beta_label_map[[dataset_name]]
    .x %>% imap(~{
      dm <- .y 
      .x %>%
        imap(~{
          marker <- .y
          # wrap_elements(
            .x %>%
              mutate(
                measurement = factor(measurement,levels=beta_map$factor_levels[[dm]]),
                depth1 = fct_rev(depth1)
              ) %>%
              group_by(measurement) %>%
              group_map(~{
                measurement <- as.character(.y$measurement)
                ggplot(.x) + 
                  geom_tile(aes(x=depth2,y=depth1,fill=dist),color="grey80") +
                  scale_fill_viridis(option="magma",name=beta_title_map[[dm]][measurement]) +
                  scale_x_discrete(position="top") + 
                  # theme_bw() +
                  theme(
                    axis.text.x = element_text(angle=25,vjust=1,hjust=0),
                    panel.grid = element_blank()
                  ) +
                  labs(x="Depth zone",y="Depth zone") + 
                  plotz_theme("transparent")
              }) %>%
               reduce(`+`) & plotz_theme("transparent") & theme(legend.position = "right")
          # ) 
        }) 
    })
  })

beta_pairs_composite$all$sim %>%
  iwalk(~{
    plotz <- .y
    ggsave(str_glue("{figure_dir}/beta_{plotz}.svg"),.x,device=svg,width=12,height=3.5,units="in")
  })


cluster_plotz <- beta_pairs %>%
  imap(~{
    .x %>%
      imap(~{
        dm <- .y 
        .x %>%
          imap(~{
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
            dend <- hclust(dd) %>%
              as.dendrogram() %>%
              color_branches(k=1,col="grey80") %>%
              color_labels(k=1,col="grey80") 
            labels_cex(dend) <- 2
            ggplot(dend %>% set("branches_lwd",3),hang=0.5,size=20) +
              theme(
                plot.caption = element_text(hjust=0.5,size=14),
              ) + 
              expand_limits(y=-0.35) + 
              ggfun::theme_transparent() +
              theme(
                # plot.background = element_rect(fill="#1a1a1a"),
                text = element_text(color="grey80"),
                # panel.background = element_rect(fill="#000000"),
                # legend.background = element_rect(fill="grey27"),
                legend.key = element_rect(fill="grey76", color = "black"),
              ) 
          }) #%>%
          # reduce(`+`) + 
          # plot_annotation(tag_levels="A") &
          # theme(plot.tag = element_text(face="bold"))
      })
  })

cluster_plotz$all$sim$metazoans
cluster_plotz$all$sim %>%
  iwalk(~{
    plotz <- .y
    ggsave(str_glue("{figure_dir}/{plotz}_cluster.svg"),.x,device=svg,width=12,height=8,units="in")
  })

ggsave(str_glue("{figure_dir}/fish_cluster_bray.svg"),cluster_plotz$all$bray$fish,device=svg,width=12,height=8,units="in")

animal_clusters <- 
  (cluster_plotz$animals$sim$inverts + theme(plot.margin = unit(c(0,1,0,0),"cm"))) +
  (cluster_plotz$animals$sim$metazoans + theme(plot.margin = unit(c(0,0,0,1),"cm"))) &
  plotz_theme("transparent") & 
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank() 
  )
ggsave(str_glue("{figure_dir}/animal_clusters.svg"),animal_clusters,device=svg,width=12,height=6,units="in")

benthic_clusters <- 
  (cluster_plotz$benthic$sim$inverts + theme(plot.margin = unit(c(0,1.5,0,0),"cm"))) +
  (cluster_plotz$benthic$sim$metazoans + theme(plot.margin = unit(c(0,0,0,1.5),"cm"))) &
  plotz_theme("transparent") & 
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank() 
  )
benthic_clusters
ggsave(str_glue("{figure_dir}/benthic_clusters.svg"),benthic_clusters,device=svg,width=12,height=7,units="in")



# accumulation ------------------------------------------------------------

ggsave(str_glue("{figure_dir}/mesophotic_accum.svg"),accum_composite,device=svg,width=12,height=4,units="in")



# depth zone ordination ---------------------------------------------------

depth_zones <- c("depth_zone","depth_zone45")
ord_zone_composite <- datasets %>%
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
                p <- plot_betadisp(.x, group=zone, method=dm, list=TRUE,binary=dm != "bray")
                
                p$plot <- p$plot +
                  scale_fill_manual(values=pal[c(1,length(pal))],name="Depth Zone",guide="none") + 
                  plotz_theme("transparent") +
                  # ggfun::theme_transparent() +
                  xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
                  ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
                  theme(legend.position = "right")
                p$plot
              }) #%>%
              # reduce(`+`) +
              # plot_layout(guides="collect") +
              # plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold"))
          })
      })
  })
ord_zone_composite$all$sim$depth_zone45$fish
fish_ord <- ord_zone_composite$all$sim$depth_zone$fish +
  ord_zone_composite$all$sim$depth_zone45$fish +
  plot_layout(guides="collect") & 
  plotz_theme("transparent")
ggsave(str_glue("{figure_dir}/fish_zone.svg"),fish_ord,device=svg,width=10,height=6,units="in")
fish_ord2 <- ord_zone_composite$all$bray$depth_zone$fish +
  ord_zone_composite$all$bray$depth_zone45$fish +
  plot_layout(guides="collect") & 
  plotz_theme("transparent")
ggsave(str_glue("{figure_dir}/fish_zone_bray.svg"),ord_zone_composite$all$bray$depth_zone45$fish,device=svg,width=8,height=6,units="in")

ord_zone_composite$all$sim$depth_zone45 %>%
  iwalk(~{
    plotz <- .y
    ggsave(str_glue("{figure_dir}/{plotz}_ord_all.svg"),.x,device=svg,width=12,height=8,units="in")
  })

# sampling depth ordination -----------------------------------------------

ord_composite <- datasets %>%
  map(~{
    dataset <- .x
    distance_methods %>%
      set_names() %>%
      map(~{
        dm <- .x
        dataset %>%
          imap(~{
            title <- plot_text2[.y]
            p <- plot_betadisp(.x, group="depth_f", method=dm, list=TRUE, expand=TRUE,binary=dm != "bray")
            p$plot <- p$plot +
              scale_fill_manual(values=pal,name="Depth",guide="none") +
              plotz_theme("transparent") +
              xlab(str_glue("Principle Coordinate 1 ({scales::percent(p$x_var,accuracy=0.1)})")) +
              ylab(str_glue("Principle Coordinate 2 ({scales::percent(p$y_var,accuracy=0.1)})")) +
              theme(
                legend.position = "right",
                axis.title = element_text(size=24)
              ) +
              guides(color="none")
            p$plot
          }) #%>%
          # reduce(`+`) +
          # plot_layout(guides="collect") +
          # plot_annotation(tag_levels = "A") &
          # theme(plot.tag = element_text(face="bold"))
      })
  })
ord_composite$all$sim$inverts
# ggsave("~/tmp/text.svg",ord_composite$all$sim$inverts,device=svg,width=12,height=8,units="in")

ord_composite$all$sim %>%
  iwalk(~{
    plotz <- .y
    ggsave(str_glue("{figure_dir}/{plotz}_ord_depth.svg"),.x,device=svg,width=12,height=8,units="in")
  })

ord_composite$all$bray$fish
ggsave(str_glue("{figure_dir}/fish_ord.svg"),ord_composite$all$sim$fish,device=svg,width=12,height=10,units="in")
ggsave(str_glue("{figure_dir}/fish_ord_bray.svg"),ord_composite$all$bray$fish,device=svg,width=12,height=10,units="in")


# depth distributions -----------------------------------------------------

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
## according to Muff & pals, a significant difference might indicate a niche shift

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
    # species = fct_reorder(species,as.numeric(factor(family)))
  ) #%>%
  # select(all_of(names(dd)))
  
# what is the most number of observations for a species?
# we need this to plot the depth ranges elegantly
# since it'll plot the line range as many time as there
# are observerations, we just duplicate the observations
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

depth_plotz <- ggplot(depth_data) + 
  geom_errorbar(aes(x=species,ymin=fb_shallow,ymax=fb_deep),width=0.5,color="grey60") +
  geom_point(aes(x=species,y=depth,size=reads,color=recode_factor(factor(str_c(depth,'m')), "0m" = "Surface")),alpha=0.7) +
  stat_summary(aes(x=species,y=depth),geom="point",fun="mean",col="grey60",size=10,shape="-") + 
  scale_color_manual(values=pal,name="Depth zone") +
  scale_size(name="Read count",labels=scales::comma,range=c(2,12),breaks=c(1e2,1e3,1e4,5e4,1e5), guide=guide_legend(override.aes = list(color="grey80"))) +
  scale_y_reverse(breaks=seq(0,300,by=30),minor_breaks=seq(0,300,by=15)) +
  scale_x_discrete(position="top") + 
  plotz_theme("transparent") +
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
ggsave(path("{figure_dir}","fish_depth.svg"),depth_plotz,device=svg,width=12,height=10,units="in")



# marker pie charts -------------------------------------------------------

marker_pal <- c("COI" = "#648FFF", "12S" = "#DC267F", "16S" = "#FFB000")
hitlist <- read_csv(here::here("data","hitlist.csv")) 
totes <- nrow(hitlist)

hl_summary <- hitlist %>% 
  mutate( across(starts_with("marker_") | starts_with("geotagged"),~.x > 0) ) %>%
  select(starts_with("marker_") | starts_with("geotagged")) %>%
  colSums() %>%
  enframe(name = "measurement",value="count") %>%
  filter(
    str_detect(measurement,"^marker_total_") | str_detect(measurement,"^geotagged_specimen_") | 
      str_detect(measurement,"^marker_specimen_in_real")
  ) %>% 
  extract(measurement,"^(.+)_(coi|16s|12s)$",into=c("measurement","marker")) %>%
  mutate(
    percent = count/totes,
    marker = str_to_upper(marker),
    marker = factor(marker,levels=c("COI","16S","12S"))
  ) %>%
  # pivot_wider(-percent,names_from = "measurement",values_from = "count") %>%
  pivot_wider(-count,names_from = "measurement",values_from = "percent") %>%
  rename(percent_total = marker_total, percent_geo_specimen = geotagged_specimen, percent_hawaii_specimen = marker_specimen_in_real) %>%
  mutate(x=4)

rot <- pi
p <- ggplot(hl_summary) + 
  geom_col(aes(x=x,y=percent_total,fill=marker)) + 
  scale_fill_manual(values=marker_pal,name="Marker",guide="none") +
  ylim(c(0,1)) +
  facet_wrap(~marker) + 
  geom_vline(xintercept = c(3.5,4.5),color="grey80") + 
  coord_polar(theta="y",start=rot,direction=1) +
  geom_text(aes(label=scales::percent(percent_total)),x=0,y=0,size=8,color="grey80") +
  # coord_polar(theta="y") + 
  xlim(c(0.2,4.5))  + 
  theme_minimal() + 
  plotz_theme("transparent") + 
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "right",
    strip.text = element_text(color="grey80",size=14)
  )
p
ggsave("{figure_dir}/marker_percent.svg",p,device=svg,width=9,height=4,units="in")

p <- ggplot(hl_summary) + 
  geom_col(aes(x=x,y=percent_geo_specimen,fill=marker)) + 
  scale_fill_manual(values=marker_pal,name="Marker",guide="none") +
  ylim(c(0,1)) + 
  facet_wrap(~marker) + 
  geom_vline(xintercept = c(3.5,4.5),color="grey80") + 
  coord_polar(theta="y",start=rot,direction=1) +
  geom_text(aes(label=scales::percent(percent_geo_specimen)),x=0,y=0,size=8,color="grey80") +
  # coord_polar(theta="y") + 
  xlim(c(0.2,4.5))  + 
  theme_minimal() + 
  plotz_theme("transparent") + 
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "right",
    strip.text = element_text(color="grey80",size=14)
  )
p
ggsave("{figure_dir}/marker_percent_geo.svg",p,device=svg,width=9,height=4,units="in")


p <- ggplot(hl_summary) + 
  geom_col(aes(x=x,y=percent_hawaii_specimen,fill=marker)) + 
  scale_fill_manual(values=marker_pal,name="Marker",guide="none") +
  ylim(c(0,1)) + 
  facet_wrap(~marker) + 
  geom_vline(xintercept = c(3.5,4.5),color="grey80") + 
  coord_polar(theta="y",start=rot,direction=1) +
  geom_text(aes(label=scales::percent(percent_hawaii_specimen,accuracy = 0.1)),x=0,y=0,size=8,color="grey80") +
  # coord_polar(theta="y") + 
  xlim(c(0.2,4.5))  + 
  theme_minimal() + 
  plotz_theme("transparent") + 
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "right",
    strip.text = element_text(color="grey80",size=14)
  )
p
ggsave("{figure_dir}/marker_percent_hawaii.svg",p,device=svg,width=9,height=4,units="in")


fishes <- read_csv("data/mundy_fish_list.csv") %>%
  mutate(
    occurrence = case_when(
      str_detect(occurrence,"endemic") ~ "endemic",
      str_detect(occurrence,"native") ~ "native",
      str_detect(occurrence,"introduced") ~ "native",
      TRUE ~ occurrence
    )
  )

fishes %>%
  count(occurrence)

ccas <- datasets$complete$all %>%
  keep_names(~.x != "fish") %>%
  map(~{
    cc <- .x
    # cc <- subset_samples(cc,station_grouping != "Honaunau North")
    ord <- ordinate(
      cc,
      method = "CCA",
      distance = "jaccard",
      binary = TRUE,
      formula = ~ depth + temp + station_grouping
      # formula = ~ depth + temp
    )
    p <- plot_ordination(cc,ord,type="samples",color="depth_f",shape="station_grouping") +
      geom_point(size=8) + 
      # ggtitle(str_glue("PCoA ordination: {plot_text[.y]}")) + 
      # scale_color_manual(values=rev(beyonce_palette(75))[5:12], name="Depth Zone") + 
      scale_color_manual(values=pal, name="Depth Zone") +
      scale_shape(name="Sampling Site") + 
      plotz_theme("light") +
      guides(shape=guide_legend(override.aes = list(color="#cccccc")))
    arrowmat <- vegan::scores(ord, display = "bp")
    # Add labels, make a data.frame
    arrowdf <- arrowmat %>%
      as_tibble(rownames = "labels") %>%
      mutate(labels = str_replace(labels,"station_grouping",""))
    # Define the arrow aesthetic mapping
    arrow_map <- aes(xend = CCA1, 
                     yend = CCA2, 
                     x = 0, 
                     y = 0, 
                     shape = NULL, 
                     color = NULL, 
                     label = labels)
    label_map <- aes(x = 1.3 * CCA1, 
                     y = 1.3 * CCA2, 
                     shape = NULL, 
                     color = NULL, 
                     label = labels)
    arrowhead = arrow(length = unit(0.02, "npc"))
    pp <- p + 
      geom_segment(
        mapping = arrow_map, 
        size = .5, 
        data = arrowdf, 
        color = "gray", 
        arrow = arrowhead
      ) + 
      geom_text(
        mapping = label_map, 
        size = 10,  
        data = arrowdf, 
        show.legend = FALSE,
        color='#cccccc'
      ) + 
      plotz_theme("transparent") +
      theme(legend.position="right", axis.title = element_text(size=24))
    pp
  })
ccas %>%
  iwalk(~{
    plotz <- .y
    ggsave(str_glue("{figure_dir}/cca_{plotz}.svg"),.x,device=svg,width=12,height=8,units="in")
  })

caps <- datasets$complete$all %>%
  # keep_names(~.x != "fish") %>%
  map(~{
    cc <- .x
    # cc <- subset_samples(cc,station_grouping != "Honaunau North")
    cc <- ps_standardize(cc,"relative")
    ord <- ordinate(
      cc,
      method = "CAP",
      distance = "bray",
      binary = TRUE,
      formula = ~ depth + temp + station_grouping
      # formula = ~ depth + temp
    )
    p <- plot_ordination(cc,ord,type="samples",color="depth_f",shape="station_grouping") +
      geom_point(size=8) + 
      # ggtitle(str_glue("PCoA ordination: {plot_text[.y]}")) + 
      # scale_color_manual(values=rev(beyonce_palette(75))[5:12], name="Depth Zone") + 
      scale_color_manual(values=pal, name="Depth Zone") +
      scale_shape(name="Sampling Site") + 
      plotz_theme("light") +
      guides(shape=guide_legend(override.aes = list(color="#cccccc")))
    arrowmat <- vegan::scores(ord, display = "bp")
    # Add labels, make a data.frame
    arrowdf <- arrowmat %>%
      as_tibble(rownames = "labels") %>%
      mutate(labels = str_replace(labels,"station_grouping",""))
    # Define the arrow aesthetic mapping
    arrow_map <- aes(xend = CAP1, 
                     yend = CAP2, 
                     x = 0, 
                     y = 0, 
                     shape = NULL, 
                     color = NULL, 
                     label = labels)
    label_map <- aes(x = 1.3 * CAP1, 
                     y = 1.3 * CAP2, 
                     shape = NULL, 
                     color = NULL, 
                     label = labels)
    arrowhead = arrow(length = unit(0.02, "npc"))
    pp <- p + 
      geom_segment(
        mapping = arrow_map, 
        size = .5, 
        data = arrowdf, 
        color = "gray", 
        arrow = arrowhead
      ) + 
      geom_text_repel(
        mapping = label_map, 
        size = 9,  
        data = arrowdf, 
        show.legend = FALSE,
        color='#cccccc'
      ) + 
      plotz_theme("transparent") +
      theme(legend.position="right", axis.title = element_text(size=24))
    pp +
      stat_ellipse(aes(group=depth_zone45))
  })
caps$fish
caps %>% reduce(`/`)
caps$inverts 
caps$metazoans 

caps %>%
  iwalk(~{
    plotz <- .y
    ggsave(str_glue("{figure_dir}/cap_{plotz}.pdf"),.x,device=cairo_pdf,width=12,height=8,units="in")
  })
  
plotz <- comm_ps %>%
  keep_names(~.x != "fish") %>%
  map(~{
    pdata <- psmelt(.x) %>%
      as_tibble() %>%
      mutate(
        thing = case_when(
          kingdom == "Metazoa" & is.na(phylum) ~ "Unidentified metazoan",
          kingdom == "Metazoa" & !is.na(phylum) ~ "Known animal phylum",
          !is.na(kingdom) ~ "Non-animal phylum",
          is.na(kingdom) ~ "Unidentified eukaryote"
        )
      ) %>%
      group_by(depth_zone45,thing) %>%
      arrange(thing) %>%
      mutate(
        thing = factor(thing),
        thing = fct_relevel(thing,"Unidentified eukaryote","Unidentified metazoan","Known animal phylum")
      ) %>%
      summarise(reads=sum(Abundance)) %>%
      ungroup() %>%
      group_by(depth_zone45) %>%
      summarise(reads=reads/sum(reads),thing=thing) %>%
      mutate(x=4)
      # replace_na(list(="Unidentified Eukaryote"))
    marker_pal <- c("Unidentified eukaryote" = "#648FFF", "Unidentified metzoan" = "#DC267F", "Known animal phylum" = "#FFB000","Non-animal phylum" = "#65E9E1")
    
    rot <- pi
    p <- ggplot(pdata) + 
      geom_col(aes(x=x,y=reads,fill=thing)) + 
      scale_fill_manual(values=marker_pal,name="Taxa") +
      ylim(c(0,1)) +
      facet_wrap(~depth_zone45) +
      geom_vline(xintercept = c(3.5,4.5),color="grey80") + 
      coord_polar(theta="y",start=rot,direction=1) +
      # geom_text(aes(label=scales::percent(percent_total)),x=0,y=0,size=8,color="grey80") +
      # coord_polar(theta="y") + 
      xlim(c(0.2,4.5))  + 
      theme_minimal() + 
      plotz_theme("transparent") + 
      theme(
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=24),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right",
        strip.text = element_text(color="grey80",size=14)
      )
    p
  })

plotz %>%
  imap(~{
    ggsave(str_glue("{figure_dir}/{.y}_unidentified.pdf"),.x,device=cairo_pdf,width=9,height=4,units="in")
  })
