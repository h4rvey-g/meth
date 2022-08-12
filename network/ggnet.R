library(igraph)
library(psych)
library(WGCNA)
library(OTUtable)
library(microbiome)
library(ggClusterNet) # remotes::install_github("taowenmicro/ggClusterNet")
library(tidyverse)
library(ggnewscale)
library(ggrepel)
# setwd(dir = "./network")
#导入数据
otu_s <- read_tsv("../tax_count.S.norm")
otu_p <- read_csv("../tax_count.P.norm.csv")
otu_g <- read_csv("../tax_count.G.norm.csv")
group <- read_tsv("../metadata.txt") %>%
  column_to_rownames("SampleID") %>%
  sample_data()
# r cutoff可选择MENA推荐的值,也可都设为0.6;出图的区别不大
# r_cut <- c(0.51, 0.51, 0.43, 0.33, 0.46, 0.68, 0.41, 0.33)
r_cut <- rep(0.6, 8)
group_list <- c(
  "Meth_Acq", "Meth_Ext", "Meth_Pre", "Meth_Rein",
  "Sal_Acq", "Sal_Ext", "Sal_Pre", "Sal_Rein"
)

#------------------------------------------------------#
# 主要的函数
#------------------------------------------------------#
get_network <- function(otu, group, i, r_cut,fname) {
  otu_relative <- transform(otu %>% column_to_rownames("Taxonomy"), transform = "compositional")
  otu <- otu %>%
    filter(Taxonomy %in%
             (OTUtable::filter_taxa(otu_relative, abundance = 0.01, persistence = 3) %>% rownames())) %>%
    column_to_rownames("Taxonomy") %>%
    otu_table(taxa_are_rows = TRUE)
  physeq <- phyloseq(otu, group) # 转为phyloseq格式
  physeq <- prune_samples(x = physeq, !grepl(".*Rein00[5-9]", sample_names(physeq))) # 去除Rein005-Rein009的样本
  topsp <- names(sort(taxa_sums(physeq), TRUE)[1:500])
  physeq <- prune_taxa(topsp, physeq)
  group_list <- c(
    "Meth_Acq", "Meth_Ext", "Meth_Pre", "Meth_Rein",
    "Sal_Acq", "Sal_Ext", "Sal_Pre", "Sal_Rein"
  )
  i <- group_list[i] #获取第i个组
  oldDF <- as(sample_data(physeq), "data.frame")
  newDF <- subset(oldDF, Group == i)
  sample_data(physeq) <- sample_data(newDF)
  # physeq <- subset_samples(physeq,Group == i) 这个函数有问题,弃之
  # 计算相关性
  result <- cor_Big_micro(
    ps = physeq,
    N = ntaxa(physeq),
    r.threshold = r_cut,
    p.threshold = 0.05,
    method = "spearman"
  )
  cor <- result[[1]]
  table(cor[cor > 0])
  
  #--提取相关矩阵
  model_igraph.2 <-
    function(cor = cor,
             method = "cluster_fast_greedy",
             seed = 12,
             Top_M = 20) {
      igraph <- graph.adjacency(cor, weighted = TRUE, mode = "undirected")
      igraph <- simplify(igraph)
      bad.vs <- V(igraph)[degree(igraph) == 0]
      igraph <- delete.vertices(igraph, bad.vs)
      
      col_g <- "#C1C1C1"
      cols <-
        colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(Top_M)
      E(igraph)$correlation <- E(igraph)$weight
      E(igraph)$weight <- abs(E(igraph)$weight)
      set.seed(seed)
      if (method == "cluster_walktrap") {
        fc <- cluster_walktrap(igraph, weights = abs(E(igraph)$weight))
      }
      if (method == "cluster_edge_betweenness") {
        fc <-
          cluster_edge_betweenness(igraph, weights = abs(E(igraph)$weight))
      }
      if (method == "cluster_fast_greedy") {
        fc <- cluster_fast_greedy(igraph, weights = abs(E(igraph)$weight))
      }
      if (method == "cluster_spinglass") {
        fc <- cluster_spinglass(igraph, weights = abs(E(igraph)$weight))
      }
      V(igraph)$modularity <- membership(fc)
      V(igraph)$label <- V(igraph)$name
      V(igraph)$label <- NA
      modu_sort <- V(igraph)$modularity %>%
        table() %>%
        sort(decreasing = T)
      top_num <- Top_M
      modu_name <- names(modu_sort[1:Top_M])
      modu_cols <- cols[1:length(modu_name)]
      names(modu_cols) <- modu_name
      V(igraph)$color <- V(igraph)$modularity
      V(igraph)$color[!(V(igraph)$color %in% modu_name)] <- col_g
      V(igraph)$color[(V(igraph)$color %in% modu_name)] <-
        modu_cols[match(V(igraph)$color[(V(igraph)$color %in%
                                           modu_name)], modu_name)]
      V(igraph)$frame.color <- V(igraph)$color
      E(igraph)$color <- col_g
      for (i in modu_name) {
        col_edge <- cols[which(modu_name == i)]
        otu_same_modu <-
          V(igraph)$name[which(V(igraph)$modularity ==
                                 i)]
        E(igraph)$color[(data.frame(as_edgelist(igraph))$X1 %in%
                           otu_same_modu) &
                          (data.frame(as_edgelist(igraph))$X2 %in%
                             otu_same_modu)] <- col_edge
      }
      sub_net_layout <-
        layout_with_fr(igraph, niter = 999, grid = "nogrid")
      data <- as.data.frame(sub_net_layout)
      data$OTU <- igraph::get.vertex.attribute(igraph)$name
      colnames(data) <- c("X1", "X2", "elements")
      tem <- V(igraph)$modularity
      tem[!tem %in% modu_name] <- "mini_model"
      tem[tem %in% modu_name] <-
        paste("model_", tem[tem %in% modu_name],
              sep = "")
      row.names(data) <- data$elements
      dat <- data.frame(
        orig_model = V(igraph)$modularity,
        model = tem,
        color = V(igraph)$color,
        OTU = igraph::get.vertex.attribute(igraph)$name,
        X1 = data$X1,
        X2 = data$X2
      )
      return(list(data, dat, igraph))
    }
  result2 <- model_igraph.2(cor = cor,
                            method = "cluster_fast_greedy",
                            seed = 12)
  # 节点的信息
  node <- result2[[1]]
  head(node)
  
  # browser()
  # 添加domain, phylum, cpp交集
  dat <- result2[[2]]
  if (fname == "s") {
    dat <- dat %>% mutate(
      domain = str_extract(OTU, "(?<=d_)[[:alpha:]]*(?=\\|+)"),
      phylum = str_extract(OTU, "(?<=p_)[[:alpha:]]*(?=\\|+)")
    )
  } else if (fname == "g") {
    dat <- dat %>% mutate(
      domain = str_extract(OTU, "(?<=d_)[[:alpha:]]*(?=\\|+)"),
      phylum = str_extract(OTU, "(?<=p_)[[:alpha:]]*(?=\\|+)")
    )
  } else {
    dat <- dat %>% mutate(
      domain = str_extract(OTU, "(?<=d_)[[:alpha:]]*(?=\\|+)"),
      phylum = str_extract(OTU, "(?<=p_)[[:alpha:]]*")
    )
  }

  
  head(dat)
  tem <- data.frame(mod = dat$model, col = dat$color) %>%
    dplyr::distinct(mod, .keep_all = TRUE)
  col <- tem$col
  names(col) <- tem$mod
  
  # ---node节点注释#-----------
  # otu_table <- as.data.frame(t(vegan_otu(physeq)))
  # tax_table <- as.data.frame(vegan_tax(physeq))
  # nodes <- nodeadd(plotcord = node, otu_table = otu_table, tax_table = tax_table)
  # head(nodes)
  #-----计算边#--------
  edge <- edgeBuild(cor = cor, node = node)
  colnames(edge)[8] <- "cor"
  head(edge)
  # 合并节点和边的信息
  tem2 <- dat %>%
    dplyr::select(OTU, model, color) %>%
    dplyr::right_join(edge, by = c("OTU" = "OTU_1")) %>%
    dplyr::rename(OTU_1 = OTU,
                  model1 = model,
                  color1 = color)
  head(tem2)
  
  tem3 <- dat %>%
    dplyr::select(OTU, model, color) %>%
    dplyr::right_join(edge, by = c("OTU" = "OTU_2")) %>%
    dplyr::rename(OTU_2 = OTU,
                  model2 = model,
                  color2 = color)
  head(tem3)
  
  tem4 <- tem2 %>% inner_join(tem3)
  head(tem4)
  #根据相关性设置边的颜色
  edge2 <- tem4 %>% mutate(color = ifelse(cor == "+", "red", "blue"))
  # 画图,geom_segment为边, point为点, geom_text为文字
  p <- ggplot() +
    geom_segment(aes(
      x = X1,
      y = Y1,
      xend = X2,
      yend = Y2,
      color = color
    ),
    data = edge2,
    size = .5) +
    scale_color_discrete(name = "Correlation",
                         labels = c("Positive", "Negative")) +
    new_scale_color() +
    geom_point(aes(X1, X2, shape = domain, color = phylum),
               data = dat,
               size = 2) +
    # scale_colour_manual(values = col) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    scale_shape_discrete(name = "Domain") +
    scale_color_discrete(name = "Phylum") +
    theme(panel.background = element_blank()) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(legend.background = element_rect(colour = NA)) +
    theme(panel.background = element_rect(fill = "white", colour = NA)) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    labs(title = paste0(i, "_",fname))
  p
  ggsave(paste0("./result/network_plot/",i, "_",fname, ".pdf"), p, width = 13, height = 12)
  #
  #
  
  
  # result4 = nodeEdge(cor = cor)
  # #提取变文件
  # edge = result4[[1]]
  # #--提取节点文件
  # node = result4[[2]]
  # browser()
  # igraph  = igraph::graph_from_data_frame(edge, directed = FALSE, vertices = node)
  # res = ZiPiPlot(igraph = igraph, method = "cluster_fast_greedy")
  # p <- res[[1]]
}
# get_network(otu = otu, group = group, physeq = ps, i = 1, r_cut = r_cut[1])
# temp <- get_dissimilarity(otu = otu, group = group, physeq = ps, i = 1, r_cut = r_cut[1])

# ------------------------------------------------------------
# 进行网络的计算
# ------------------------------------------------------------

wrap_levels <- function(otu,group,r_cut,fname) {
  total_list <- list(rep(list(otu), 8), rep(list(group), 8), 1:8,
                     r_cut, rep(list(fname), 8))
  # get_network(otu = otu, group = group, physeq = ps, i = 1, r_cut = r_cut[1], cpp = cpp)
  pwalk(total_list, get_network)
}
wrap_levels(otu = otu_s, group = group,r_cut = r_cut,fname="s")
wrap_levels(otu = otu_p, group = group,r_cut = r_cut,fname="p")
wrap_levels(otu = otu_g, group = group,r_cut = r_cut,fname="g")
