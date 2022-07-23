write_ms_table <- function(tbl,file,caption="",...) {
  if (caption != "") {
    write_lines(str_c(": ",caption),file)
    write_csv(tbl,file,append=TRUE,col_names = TRUE,...)
  } else {
    write_csv(tbl,file,...)    
  }
}
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
    Site = fct_relevel(Site,"none",after=Inf)
  ) %>%
  arrange(Site,Depth,`Sample ID`)

write_ms_table(rs,"output/reads_summary.csv","Sequencing results for mesophotic eDNA samples across assay types",na="")

