Results
================

``` r
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
  mutate(sample=SampleID) %>% 
  column_to_rownames("SampleID")
```

自变量为时期

``` r
adonis2(otu_s ~ Group,
        group,
        permutations = 999,
        distance = 'bray')
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

自变量为每个个体，结果显示No residual component，不太清楚原因

``` r
adonis2(otu_s ~ sample,
        group,
        permutations = 999,
        distance = 'bray')
```

    ## No residual component
    ## 
    ## adonis2(formula = otu_s ~ sample, data = group, permutations = 999, distance = "bray")
    ##          Df SumOfSqs R2 F Pr(>F)
    ## Model    19   1.6928  1         
    ## Residual  0   0.0000  0         
    ## Total    19   1.6928  1
