
# 加载包
library(igraph)
library(psych)
library(WGCNA)
library(OTUtable)
library(stringr)
library(dplyr)
library(microbiome)
# tax.s.rel <- read.table("sample.txt",header = T, sep="\t", comment.char="", quote = "",row.names = 1)
otu = read.table("./tax.S.full.txt",header =T,row.names=1,sep = "\t",comment.char="", quote = "")
tax.s.rel <- transform(otu,transform = "compositional")
# otu$taxonomy <- paste0("otu",1:7188)
tax.s.rel <- dplyr::filter(otu,rownames(otu)%in%rownames(OTUtable::filter_taxa(tax.s.rel,abundance=.01,persistence=3)))
group <- c("Meth_Acq","Meth_Ext","Meth_Pre","Meth_Rein",
           "Sal_Acq","Sal_Ext","Sal_Pre","Sal_Rein")
net_data <- data.frame(Edges_number="",
                       Nodes_number="",
                       "connectance"="",
                       "Average degree"="",
                       "Average path length"="",
                       "Diameter"="",
                       "edge connectivity"="",
                       "Clustering coefficient"="",
                       "Betweenness centralization"="",
                       "Degree centralization"=""
)
otu_lap <- list() #交集
for (i in 1:5) {
        otu_select = select(tax.s.rel, contains(paste0(group[i])))
        occor = corAndPvalue(t(otu_select), use = "all.obs", method =
                                     "spearman")
        occor.r = occor$cor # 取相关性矩阵R值
        occor.p = occor$p # 取相关性矩阵p值
        occor.p <- p.adjust(occor.p, method = "bonferroni")
        
        
        r_cut <- switch (i,
                         0.51, 0.51, 0.43, 0.33, 0.46, 0.68, 0.41, 0.33)
        occor.r[occor.p > 0.5 |abs(occor.r) < r_cut |abs(occor.r) == 1] = 0
        # 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
        
        
        # 构建igraph对象
        # occor.r[occor.r != 0] = 1
        igraph = graph_from_adjacency_matrix(occor.r,
                                             mode = "undirected",
                                             weighted = T,
                                             diag = FALSE)
        
        bad.vs = V(igraph)[degree(igraph) == 0]
        igraph = delete.vertices(igraph, bad.vs)
        
        
        
        # 将igraph weight属性赋值到igraph.weight
        igraph.weight = E(igraph)$weight
        
        # 做图前去掉igraph的weight权重，因为做图时某些layout会受到其影响
        E(igraph)$weight = NA
        set.seed(123)
        E.color = igraph.weight
        E.color = ifelse(E.color > 0, "red", ifelse(E.color < 0, "blue", "grey"))
        E(igraph)$color = as.character(E.color)
        taxnomy <-
                data.frame(
                        domain = str_extract(V(igraph)$name, "d_[[:alpha:]]*(?=\\|+)"),
                        otu = str_extract(V(igraph)$name, "s_.*")
                )
        library(FELLA)
        add.vertex.shape("triangle",
                         clip = shapes("circle")$clip,
                         plot = FELLA:::mytriangle)
        taxnomy$domain %>% str_replace(., "d_Bacteria", "circle") %>% str_replace(., "d_Viruses", "triangle") %>% str_replace(., "d_Archaea", "square") ->
                V(igraph)$shape
        fc = cluster_fast_greedy(igraph, weights = NULL)# cluster_walktrap cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
        modularity = modularity(igraph, membership(fc))
        # 按照模块为节点配色
        comps = membership(fc)
        colbar = rainbow(max(comps))
        V(igraph)$color = colbar[comps]
        
        # l <- layout.fruchterman.reingold(igraph, niter=5000, area=vcount(igraph)^4*10)
        # windowsFonts(TimesNewRoman  = windowsFont("Times New Roman"))
        pdf(paste0(group[i], ".pdf"),family="serif")
        # par(pin=c(10,3))
        plot(
                igraph,
                main = paste0(group[i]),
                vertex.frame.color = NA,
                vertex.label = NA,
                edge.lty = 1,
                edge.width = .5,
                edge.curved = TRUE,
                edge.color = E(igraph)$color,
                margin = c(0, 0, 0, 1),
                layout = layout.fruchterman.reingold,
                vertex.shape = V(igraph)$shape,
                vertex.size = 1
        )
        
        legend(
                "topright",
                title = as.expression(bquote(bold("Legend"))),
                c(
                        as.expression(bquote(bold(
                                "Correlation"
                        ))),
                        "Positive correlation",
                        "Negative correlation",
                        as.expression(bquote(bold("Domain"))),
                        "Archaea",
                        "Bacteria",
                        "Viruses"
                ),
                lty = c(NA, 1, 1, rep(NA, 4)),
                col = c(NA, "red", "blue", NA, rep("black", 4)),
                pch = c(rep(NA, 4), 15, 16, 17)
        )
        
        dev.off()
        
        #输出数据
        E(igraph)$weight = abs(igraph.weight)
        net_data[nrow(net_data)+1,] <-
                c(
                        length(E(igraph)),
                        length(V(igraph)),
                        edge_density(igraph, loops = FALSE),
                        mean(igraph::degree(igraph)),
                        average.path.length(igraph),
                        diameter(
                                igraph,
                                directed = FALSE,
                                unconnected = TRUE,
                                weights = NULL
                        ),
                        edge_connectivity(igraph),
                        transitivity(igraph),
                        centralization.betweenness(igraph)$centralization,
                        centralization.degree(igraph)$centralization
                )
        otu_lap <- append(otu_lap,list(V(igraph)$name))
        
        
}
net_data <- net_data[-1,]
rownames(net_data) <- group
write.table(net_data,file = "./network/network data.txt",sep = "\t",quote = F,row.names=T)
#找overlap
names(otu_lap) <- group
to_species <- function(otu_ele){
        otu_ele <- str_extract(otu_ele,pattern = "s_.*")
        return(otu_ele)
}
otu_lap2 <- rapply(otu_lap,to_species,how="replace")
otu_lap2 <- append(otu_lap2,list(cpp=names$otu))
intersections <- list()
for (i in 1:8) {
        tmp1 <- intersect(otu_lap2[[i]],otu_lap2[[9]])
        intersections <- append(intersections,list(tmp1))
}
names(intersections) <- group
# lapply(intersection, write, "./stamp/overlap/cpp_and_network_overlap.txt", append=TRUE, ncolumns=1000)
# intersection <- Reduce(intersect,otu_lap2)
# write.table(intersection,file = "./stamp/overlap/cpp_and_network_overlap.txt",sep = "\t",quote = F,row.names=T)
sink("./stamp/overlap/cpp_and_network_overlap.txt")
print(intersections)
sink()

#cpp与pathway的交集
#1 找出显著pathway和otu
sig_acq <-
        rownames(
                read.table(
                        "./stamp/overlap/associate_cpm_meth.txt",
                        header = T,
                        row.names = 1,
                        sep = "\t",
                        comment.char = "",
                        quote = ""
                )
        )
read.table(
        "./stamp/overlap/pathabundance_cpm_meth.pcl",
        header = T,
        row.names = 1,
        sep = "\t",
        comment.char = "",
        quote = ""
) %>% rownames() %>% str_split_fixed(., "\\|", 2) %>% as.data.frame() -> all_acq
all_acq <- filter(all_acq, grepl("g__", V2))
str_split_fixed(all_acq$V2, "\\.(?=s_)", 2)[, 2] %>% str_replace(., "__", "_") -> all_acq$V2
all_acq <- filter(all_acq, V1 %in% sig_acq)
str_replace_all(names$otu, "[[:punct:]]|[[:blank:]]", "_") %>% str_replace(., "__", "_") ->
        names$otu_trans
# 2 取交集
tmp1 <- inner_join(all_acq, names, by = c("V2" = "otu_trans"))
tmp1 <- tmp1[, c("V1", "otu")]
colnames(tmp1)[1] <- "pathway"
write.table(
        tmp1,
        file = "./stamp/overlap/cpp_and_pathway_METH_overlap.txt",
        sep = "\t",
        quote = F,
        row.names = F
)
