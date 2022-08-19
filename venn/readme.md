Venn图
================

``` r
library(tidyverse)
library(eulerr)
library(scales)
library(RColorBrewer)
otu_s <- read_tsv("../tax_count.S.norm") %>%
  select(-matches(".*Rein00([6-9]|10)|Sal.*|Taxonomy"))
otu_collapse <- transmute(
  rowwise(otu_s),
  Meth_Pre = max(c_across(contains("Pre"))),
  Meth_Acq = max(c_across(contains("Acq"))),
  Meth_Ext = max(c_across(contains("Ext"))),
  Meth_Rein = max(c_across(contains("Rein")))
)
otu_collapse[otu_collapse>0] <- 1
```

不设置权重时

``` r
fit_orig <- euler(otu_collapse)
fill_color <- brewer.pal(4,"Set1") 
plot(fit_orig, quantities = T,
     fills = list(fill = fill_color, alpha = 0.9),
     labels = list(col = "white", font = 2))
```

![](venn_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

都包含的太多。所以将四个组都有的设置权重0.01

``` r
weights <-
  otu_collapse %>% rowwise() %>% transmute(weights = ifelse(sum(c_across(everything(
  ))) == 4, 0.01, 1))
fit <-
  euler(otu_collapse, shape = "ellipse", weights = weights)
quantity <- tibble(
  count = fit_orig[["original.values"]],
  percent = percent(count / sum(count), accuracy = 0.01),
  all = paste(count,"(",percent,")")
)
plot(fit, quantities = list(labels=quantity$all,col = "white",cex=0.8),
     fills = list(fill = fill_color, alpha = 0.9),
     labels = list(col = "white", font = 2))
```

![](venn_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->