library(vegan)
library(tidyverse)
library(reshape2)
library(scales)
library(psych)
otu_s <- read_tsv("../tax_count.S.norm") %>%
  select(-matches(".*Rein00([6-9]|10)|Sal.*")) %>%
  column_to_rownames("Taxonomy") %>% t()
group <- read_tsv("../metadata.txt") %>% 
  filter(.,!grepl(".*Rein00([6-9]|10)|Sal.*", SampleID)) %>% 
  # mutate(sample=SampleID) %>%
  column_to_rownames("SampleID")
adonis2(otu_s ~ Group,
        group,
        permutations = 999,
        distance = 'bray')
adonis2(otu_s ~ Replicate,
        group,
        permutations = 999,
        distance = 'bray')
