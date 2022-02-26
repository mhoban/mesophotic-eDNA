library(vegan)
library(sars)

sup <- function(...) suppressWarnings(suppressMessages(...))

all_models <- map2(communities,names(communities),~{
  raw_data <- .x$raw
  specs <- sample_data %>%
    group_by(depth_f) %>%
    group_map(~{
      s <- raw_data %>% 
        select(OTU,any_of(.x$id)) %>%
        pivot_longer(matches(sample_pattern),names_to="site",values_to="reads") %>%
        pivot_wider(names_from="OTU",values_from="reads") %>%
        column_to_rownames("site") %>%
        decostand("pa",margin=NULL)
      sup(curve <- specaccum(s,method="random",permutations=1000))
      # tibble(sites=curve$sites,richness=curve$richness,depth=.y$depth_f)
      as_tibble(curve$perm) %>%
        mutate(sites=row_number()) %>%
        select(sites,everything()) %>%
        pivot_longer(-sites,names_to = "permutation", values_to = "richness") %>%
        distinct(sites,richness) %>%
        mutate(depth=.y$depth_f) %>%
        select(sites,richness,depth)
    }) %>% 
    set_names(levels(sample_data$depth_f)) %>%
    bind_rows()
  
  sup(models <- sar_average(data=specs %>% select(sites,richness)))
  # plot(models)
  # summary(models)
  
  best_mod <- names(which.min(models$details$ics))
  sup(predictor <- do.call(str_c("sar_",best_mod),list(data=specs %>% select(sites,richness), grid_start = "exhaustive", grid_n = 1000)))
  
  # summary(predictor)
  # plot(predictor)
  
  # predictor <- models
  predictions <- sar_pred(predictor,c(1:1000))
  slopes <- diff(predictions$Prediction)/diff(predictions$Area)
  slopes <- c(1,slopes)
  
  # plot(Prediction ~ Area, data=predictions)
  ggplot(predictions) + 
    geom_line(aes(x=Area,y=Prediction))
  
  asymp <- predictions$Area[slopes < 1][1]
  asymp_taxa <- predictions$Prediction[slopes < 1][1]
  taxa_6 <- predictions$Prediction[6]
  percent_taxa <- taxa_6/asymp_taxa
  cat("Community:",.y,"\n")
  cat("Best model:",best_mod,"\n")
  cat("Asymptote:",asymp,"\n")
  cat("Taxa at asymptote:",asymp_taxa,"\n")
  cat("Taxa at 6 replicates:",taxa_6,"\n")
  cat("Percent of taxa represented by 6 replicates:",percent_taxa,"\n")
  cat("\n\n")
  list(models = models, best_model_name = best_mod, best_model = predictor, asymptote=asymp, asymptote_taxa=asymp_taxa, our_taxa = taxa_6, accum_data=specs)
  # plot(Prediction ~ Area, data=predictions)
})


all_models <- map2(comm_ps,names(comm_ps),~{
  raw_data <- otu_table(.x) %>%
    as("matrix") %>%
    as_tibble(rownames="sample")
  # plotted <<- FALSE
  specs <- sample_data %>%
    group_by(depth_f) %>%
    group_map(~{
      s <- raw_data %>% 
        filter(sample %in% .x$id) %>%
        column_to_rownames("sample") %>%
        decostand("pa",margin=NULL)
      sup(curve <- specaccum(s,method="random",permutations=1000))
      # plot(curve, add = plotted, col = "deepskyblue3", ylab="Number of taxa", ci.lty = 0, lty = 1, xlab = "Number of replicates", xlim=c(0,6), ylim=c(0,600))
      # if (!plotted) {
      #   plotted <<- TRUE
      # }
      # list(spec=curve,depth=.y$depth_f,data=s)
      tibble(sites=curve$sites,richness=curve$richness,depth=.y$depth_f)
    }) %>% 
    set_names(levels(sample_data$depth_f)) %>%
    bind_rows()
  
  sup(models <- sar_average(data=specs %>% select(sites,richness)))
  # plot(models)
  # summary(models)
  
  best_mod <- names(which.min(models$details$ics))
  sup(predictor <- do.call(str_c("sar_",best_mod),list(data=specs %>% select(sites,richness), grid_start = "exhaustive", grid_n = 1000)))
  
  # summary(predictor)
  # plot(predictor)
  
  # predictor <- models
  predictions <- sar_pred(predictor,c(1:1000))
  slopes <- diff(predictions$Prediction)/diff(predictions$Area)
  slopes <- c(1,slopes)
  
  # plot(Prediction ~ Area, data=predictions)
  
  asymp <- predictions$Area[slopes < 1][1]
  asymp_taxa <- predictions$Prediction[slopes < 1][1]
  taxa_6 <- predictions$Prediction[6]
  percent_taxa <- taxa_6/asymp_taxa
  cat("Community:",.y,"\n")
  cat("Best model:",best_mod,"\n")
  cat("Asymptote:",asymp,"\n")
  cat("Taxa at asymptote:",asymp_taxa,"\n")
  cat("Taxa at 6 replicates:",taxa_6,"\n")
  cat("Percent of taxa represented by 6 replicates:",percent_taxa,"\n")
  cat("\n\n")
  list(models = models, best_model_name = best_mod, best_model = predictor, asymptote=asymp, asymptote_taxa=asymp_taxa, our_taxa = taxa_6, accum_data=specs)
  # plot(Prediction ~ Area, data=predictions)
})


accum_data <- comm_ps %>%
  map(~{
    raw_data <- otu_table(.x) %>%
      as("matrix") %>%
      as_tibble(rownames="sample")
    # plotted <<- FALSE
    sample_data %>%
      group_by(depth_f) %>%
      group_map(~{
        s <- raw_data %>% 
          filter(sample %in% .x$id) %>%
          column_to_rownames("sample") %>%
          decostand("pa",margin=NULL)
        sup(curve <- specaccum(s,method="random",permutations=1000))
        # plot(curve, add = plotted, col = "deepskyblue3", ylab="Number of taxa", ci.lty = 0, lty = 1, xlab = "Number of replicates", xlim=c(0,6), ylim=c(0,600))
        # if (!plotted) {
        #   plotted <<- TRUE
        # }
        # list(spec=curve,depth=.y$depth_f,data=s)
        tibble(sites=curve$sites,richness=curve$richness,depth=.y$depth_f)
      }) %>% 
      set_names(levels(sample_data$depth_f)) %>%
      bind_rows()
  })

accum_data <- map2(communities,names(communities),~{
  raw_data <- .x$relative
  sample_data %>%
    group_by(depth_f) %>%
    group_map(~{
      s <- raw_data %>% 
        select(OTU,any_of(.x$id)) %>%
        pivot_longer(matches(sample_pattern),names_to="site",values_to="reads") %>%
        pivot_wider(names_from="OTU",values_from="reads") %>%
        column_to_rownames("site") %>%
        decostand("pa",margin=NULL)
      sup(curve <- specaccum(s,method="random",permutations=1000))
      # plot(curve, add = plotted, col = "deepskyblue3", ylab="Number of taxa", ci.lty = 0, lty = 1, xlab = "Number of replicates", xlim=c(0,6), ylim=c(0,350))
      if (!plotted) {
        plotted <<- TRUE
      }
      # list(spec=curve,depth=.y$depth_f,data=s)
      tibble(sites=curve$sites,richness=curve$richness,depth=.y$depth_f)
    }) %>% 
    set_names(levels(sample_data$depth_f)) %>%
    bind_rows()
})
 
names_map <- c("m16S" = "fishes", "m18S" = "inverts", "mCOI" = "metazoans")
plotz <- all_models %>%
  map2(names(.),~{
    ggplot(.x$accum_data) + 
      geom_line(aes(x=sites,y=richness,color=depth)) + 
      ggtitle(.y)
  })
relative_plotz



# using nls ---------------------------------------------------------------

# using raw data
all_models <- map2(communities,names(communities),~{
  raw_data <- .x$raw
  specs <- sample_data %>%
    group_by(depth_f) %>%
    group_map(~{
      s <- raw_data %>% 
        select(OTU,any_of(.x$id)) %>%
        pivot_longer(matches(sample_pattern),names_to="site",values_to="reads") %>%
        pivot_wider(names_from="OTU",values_from="reads") %>%
        column_to_rownames("site") %>%
        decostand("pa",margin=NULL)
      sup(curve <- specaccum(s,method="random",permutations=1000))
      # tibble(sites=curve$sites,richness=curve$richness,depth=.y$depth_f)
      as_tibble(curve$perm) %>%
        mutate(sites=row_number()) %>%
        select(sites,everything()) %>%
        pivot_longer(-sites,names_to = "permutation", values_to = "richness") %>%
        distinct(sites,richness) %>%
        mutate(depth=.y$depth_f) %>%
        select(sites,richness,depth)
    }) %>% 
    set_names(levels(sample_data$depth_f)) %>%
    bind_rows()
  
  try_mod <- function(e) {
    tryCatch(
      e,
      error = function(e) NULL
    )
  }
  
  modls <- list(
    lomolino = try_mod(nls(richness ~ SSlomolino(sites, Asym, xmid, slope),  data=specs)),
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
  
  slopes <- diff(predictions$y)/diff(predictions$x)
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

# using normalized datasets
all_models <- map2(comm_ps,names(comm_ps),~{
  raw_data <- otu_table(.x) %>%
    as("matrix") %>%
    as_tibble(rownames="sample")
  specs <- sample_data %>%
    group_by(depth_f) %>%
    group_map(~{
      s <- raw_data %>% 
        filter(sample %in% .x$id) %>%
        column_to_rownames("sample") %>%
        decostand("pa",margin=NULL)
      sup(curve <- specaccum(s,method="random",permutations=1000))
      as_tibble(curve$perm) %>%
        mutate(sites=row_number()) %>%
        select(sites,everything()) %>%
        pivot_longer(-sites,names_to = "permutation", values_to = "richness") %>%
        distinct(sites,richness) %>%
        mutate(depth=.y$depth_f) %>%
        select(sites,richness,depth)
    }) %>% 
    set_names(levels(sample_data$depth_f)) %>%
    bind_rows()
  
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
  
  slopes <- diff(predictions$y)/diff(predictions$x)
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
