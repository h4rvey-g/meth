library(clusterProfiler)
library(tidyverse)
library(readxl)
library(DOSE)
library(ggh4x)
dat <- read_xlsx("./KEGG enrichment.xlsx") %>%
  mutate(
    GeneRatio = parse_ratio(GeneRatio)
  ) %>%
  arrange(level)
p <- ggplot(dat, aes(x = GeneRatio, y = interaction(Level_3, Level_2, Level_1, sep = "&"))) +
  geom_point(aes(size = Count, color = -log10(pvalue))) +
  guides(
    y = guide_axis_nested(delim = "&")
  ) +
  theme_classic()
ggsave("./KEGG enrichment.pdf", p, width = 25, height = 10)
