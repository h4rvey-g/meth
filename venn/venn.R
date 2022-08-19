library(tidyverse)
library(eulerr)
otu_s <- read_tsv("../tax_count.S.norm") %>%  
    select(-Taxonomy)
otu_s[otu_s>0] <- 1
fit <- euler(otu_s)
