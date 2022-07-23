comm_ps %>%
  walk2(names(.),~{
    ntax <- ntaxa(.x)
    assay <- plot_text[.y]
    cat(str_glue("{assay}:"),"\n")
    cat(str_glue("{ntax} zOTUs"),"\n")
    taxa <- taxa_tibble(.x) %>%
      arrange(domain,kingdom,phylum,class,order,family,genus,species) %>%
      mutate(
        order = case_when(
          # some manual tweaks
          family == "Pomacentridae" ~ "Ovalentaria",
          TRUE ~ order
        ),
        family = case_when(
          genus %in% c("Calotomus","Chlorurus","Scarus") ~ "Scarinae",
          TRUE ~ family
        )
      ) 
    read_summary <- psmelt(.x) %>%
      mutate(
        order = case_when(
          # some manual tweaks
          family == "Pomacentridae" ~ "Ovalentaria",
          TRUE ~ order
        ),
        family = case_when(
          genus %in% c("Calotomus","Chlorurus","Scarus") ~ "Scarinae",
          TRUE ~ family
        )
      ) %>%
      group_by(family) %>%
      summarise(zotus = length(unique(OTU)), readcount = sum(Abundance)) %>%
      arrange(desc(readcount)) %>%
      mutate(readcount=scales::number(readcount,big.mark = ","))
    print(read_summary)
    cat("number of zotus identified to:\n")
    print(taxa %>%
            mutate(across(domain:species,~!is.na(.x))) %>%
            select(domain:species) %>%
            colSums())
    cat("\nidentification methods:\n")
    print(table(taxa$id_method))
    cat("\nnumber of identified taxonomic categories:\n")
    print(taxa %>%
            mutate(across(domain:species,~!is.na(.x) & !duplicated(.x))) %>%
            select(domain:species) %>%
            colSums())
    cat("------------------------------------------------------\n\n")
  })


.x <- comm_ps$metazoans
f <- psmelt(.x)
ff <- f %>% 
  pivot_longer(domain:species,names_to="level",values_to="taxon") %>%
  group_by(level,taxon) %>%
  summarise(reads = sum(Abundance))
ff
f %>%
  filter(kingdom == "Metazoa") %>%
  group_by(domain,kingdom,phylum,class,order,family,genus,species) %>%
  summarise(reads = sum(Abundance,na.rm=TRUE),zotus=length(unique(OTU))) %>%
  arrange(desc(reads)) %>%
  View()

View(ff)
  
read_summary <- psmelt(.x) %>%
  mutate(
    order = case_when(
      # some manual tweaks
      family == "Pomacentridae" ~ "Ovalentaria",
      TRUE ~ order
    ),
    family = case_when(
      genus %in% c("Calotomus","Chlorurus","Scarus") ~ "Scarinae",
      TRUE ~ family
    )
  ) %>%
  group_by(family) %>%
  summarise(zotus = length(unique(OTU)), readcount = sum(Abundance)) %>%
  arrange(desc(readcount)) %>%
  mutate(readcount=scales::number(readcount,big.mark = ","))