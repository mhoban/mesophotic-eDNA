library(tidyverse)
library(dbplyr)

coacross <- function(...) {
  coalesce(!!!across(...))
}


marker <- "coi"


worms <- read_delim("~/projects/barcoding/fishes/data/worms_dump.tsv",delim="\t",quote='"',escape_double=TRUE,col_types = str_c(rep("c",32),collapse=""))
worms <- worms %>%
  select(scientificName, taxonomicStatus, taxonRank, parentNameUsage, kingdom:specificEpithet)

blast <- read_tsv(str_glue("working/flow-output/mesophotic/{marker}/07_blast_seq/seq_blast_Result.tab"),
                  col_names = c("qseqid","sseqid","staxids","sscinames","scomnames","sskingdoms",
                                "pident","length","qlen","slen","mismatch","gapopen","gaps","qstart",
                                "qend","sstart","send","stitle","evalue","bitscore","qcovs","qcovhsp")) %>%
  select(zotu=qseqid, seq_id=sseqid, taxid=staxids, match=pident, length, query_cov=qcovs)  %>%
  separate_rows(taxid,sep=";") %>%
  mutate(taxid=as.numeric(taxid))

ranked_lineages <- read_delim(
  "data/taxdump/rankedlineage.dmp",delim="|",
  col_names=c("taxid","name","species","genus","family","order","class","phylum","kingdom","superkingdom"), 
  col_types="nccccccccc")

otus <- read_tsv(str_glue("working/flow-output/mesophotic/{marker}_collapsed_taxa.tab")) %>%
  select(domain:species,OTU) %>%
  rowwise() %>%
  mutate(dropped=sum(c_across(domain:species) == "dropped"))


ot <- otus %>%
  select(OTU,dropped) %>%
  filter(dropped >= 5) %>%
  left_join(blast,by=c("OTU" = "zotu")) %>%
  inner_join(ranked_lineages, by="taxid") %>%
  rowwise() %>%
  mutate(taxa = sum(!is.na(c_across(species:kingdom)))) %>% 
  filter(taxa >= 3) %>%
  select(zotu=OTU,dropped,taxa,match,taxid,name:superkingdom) %>%
  left_join(worms %>% filter(taxonRank == "Genus", taxonomicStatus == "accepted") , by=c("genus" = "scientificName"), suffix=c("","_family")) %>%
  left_join(worms %>% filter(taxonRank == "Family", taxonomicStatus == "accepted"), by=c("family" = "scientificName"), suffix=c("","_order")) %>%
  left_join(worms %>% filter(taxonRank == "Order", taxonomicStatus == "accepted") , by=c("order" = "scientificName"), suffix=c("","_class")) %>%
  left_join(worms %>% filter(taxonRank == "Class", taxonomicStatus == "accepted") , by=c("class" = "scientificName"), suffix=c("","_phylum")) %>%
  left_join(worms %>% filter(taxonRank == "Phylum" | taxonRank == "Phylum (Division)", taxonomicStatus == "accepted") %>% select(scientificName,kingdom), by=c("phylum" = "scientificName"), suffix=c("","_kingdom")) %>%
  mutate(
    family = coacross(starts_with("family")),
    order = coacross(starts_with("order")),
    class = coacross(starts_with("class")),
    phylum = coacross(starts_with("phylum")),
    kingdom = coacross(starts_with("kingdom"))
  ) %>%
  select(1:superkingdom)

View(ot)
View(otus %>% filter(OTU %in% ot$zotu))

most_frequent <- function(x, min=1, mult = FALSE) {
  if (!is.factor(x)) x <- factor(x)
  A <- tabulate(x)
  if (isTRUE(mult)) {
    levels(x)[A == max(A)]
  } 
  idx <- which.max(A)
  # mf <- levels(x)[which.max(A)]
  if (min > 1) {
    if (A[idx] >= min) {
      return(levels(x)[idx])
    }  else {
      return(NA)
    }
  } else {
    return(levels(x)[idx])
  }
}


ot_merged <- ot %>% 
  # head(1000) %>%
  group_by(zotu) %>%
  summarise(
    across(
      species:superkingdom,
      most_frequent
    ),
    matches = n()
  ) 

otu_table <- read_tsv(str_glue("working/flow-output/mesophotic/{marker}_collapsed_taxa.tab"))

otutus <- otu_table %>%
  mutate(
    across(domain:species, ~ifelse(.x == "dropped",NA,.x))
  ) %>%
  left_join(ot_merged, by=c("OTU" = "zotu"), suffix=c("","_merge")) %>%
  mutate(
    species = coacross(starts_with("species")),
    genus = coacross(starts_with("genus")),
    family = coacross(starts_with("family")),
    order = coacross(starts_with("order")),
    class = coacross(starts_with("class")),
    phylum = coacross(starts_with("phylum")),
    kingdom = coacross(starts_with("kingdom")),
    domain = coacross(starts_with("domain")),
  ) %>%
  select(names(otu_table)) %>%
  mutate(
    across(domain:species, ~ifelse(is.na(.x),"dropped",.x))
  )
write_tsv(otutus,str_glue("working/flow-output/mesophotic/{marker}_merged.tab"))

ff <- tribble(
  ~one,~two,~three,~four,
  "won",NA,NA,"won",
  "tow",NA,NA,NA,
  NA,"thar",NA,NA,
  "foor","foor","foor","foor",
  NA,NA,NA,NA
)

fff <- ff %>%
  rowwise() %>% 
  mutate(real_thing = list(na.omit(c_across(one:four))))


df <- tibble(a = c(1, 2), b = c(4, 3), c = c(5, 7))
map_dfc(.x = c(combn(rev(names(df)), 2, simplify = FALSE),
               combn(names(df), 2, simplify = FALSE)),
        ~ df %>%
          rowwise() %>%
          transmute(!!paste(.x, collapse = "_") := reduce(c_across(all_of(.x)), `-`)) %>%
          ungroup())


ff %>%
  # rowwise() %>% 
  mutate(eal_thing = coalesce(!!!vars_select(one:four)))

f <- list(apple="one",button="TWO",sockface="three")
ff <- map(f,str_to_upper)
