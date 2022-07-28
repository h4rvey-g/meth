library(vegan)
library(tidyverse)
otu <- read_tsv("./tax.S.txt")
group <- read_tsv("./metadata.txt")
