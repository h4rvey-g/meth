---
title: "Results"
output: github_document
---



```r
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
  column_to_rownames("SampleID")
```
自变量为时期

```r
adonis2(otu_s ~ Group,
        group,
        permutations = 999,
        distance = 'bray')
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = otu_s ~ Group, data = group, permutations = 999, distance = "bray")
##          Df SumOfSqs      R2      F Pr(>F)    
## Group     3  0.76748 0.45339 4.4238  0.001 ***
## Residual 16  0.92528 0.54661                  
## Total    19  1.69276 1.00000                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
自变量为每个个体

```r
adonis2(otu_s ~ Replicate,
        group,
        permutations = 999,
        distance = 'bray')
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = otu_s ~ Replicate, data = group, permutations = 999, distance = "bray")
##           Df SumOfSqs      R2      F Pr(>F)
## Replicate  1  0.12999 0.07679 1.4972  0.193
## Residual  18  1.56278 0.92321              
## Total     19  1.69276 1.00000
```


