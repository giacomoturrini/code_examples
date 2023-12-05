library(Seurat)

options(future.globals.maxSize = 10000000000)
options(Seurat.object.assay.version = "v5")

# Load R image with astros ----
all_astro_merge <- readRDS("all_astrocytes_integrated.Rds")

## Insert the information from schist clustering ----
level_2 <- read.csv("schist_level_2.csv", header=FALSE)
schist_vector <- setNames(level_2$V2, level_2$V1)
all_astro_merge[["schist_clusters"]] <- schist_vector[row.names(all_astro_merge@meta.data)]
Idents(all_astro_merge) <- 'schist_clusters'

## Search for DEG with schist clusters ----
DefaultAssay(all_astro_merge) <- "RNA"
all_astro_merge <- JoinLayers(all_astro_merge)
markers_schist <- FindAllMarkers(all_astro_merge, only.pos = TRUE, logfc.threshold = 0.1)
markers_schist <-markers_schist[order(markers_schist$avg_log2FC, decreasing=TRUE),]
markers_schist <- split(markers_schist, f = markers_schist$cluster)

# Write xlsx files for the markers
library(writexl)
write_xlsx(markers_schist[[1]],"/mnt/remo/prj/transcriptome_cerebellum/dataset/3fold_integ/schistmarkers0.xlsx")
write_xlsx(markers_schist[[2]],"/mnt/remo/prj/transcriptome_cerebellum/dataset/3fold_integ/schistmarkers1.xlsx")
write_xlsx(markers_schist[[3]],"/mnt/remo/prj/transcriptome_cerebellum/dataset/3fold_integ/schistmarkers2.xlsx")
write_xlsx(markers_schist[[4]],"/mnt/remo/prj/transcriptome_cerebellum/dataset/3fold_integ/schistmarkers3.xlsx")
write_xlsx(markers_schist[[5]],"/mnt/remo/prj/transcriptome_cerebellum/dataset/3fold_integ/schistmarkers4.xlsx")

## Plot cluster dataset composition ----
# Name clusters as in thesis
all_astro_merge@meta.data[which(all_astro_merge@meta.data$datasets == "Goetz"),]$datasets <- "dLMU"
all_astro_merge@meta.data[which(all_astro_merge@meta.data$datasets == "Telley_FC"),]$datasets <- "dLau"
all_astro_merge@meta.data[which(all_astro_merge@meta.data$datasets == "Macosko"),]$datasets <- "Kozareva"

library(dplyr)
library(ggplot2)

#Count cells by Identity 
table1<- all_astro_merge@meta.data %>% group_by(cca_clusters) %>% summarise(n.cells=n())

knitr::kable(table1,caption = "Number of Cluster Members")


#Count cells by other cluster level
table<- all_astro_merge@meta.data %>% group_by(cca_clusters,datasets) %>% summarise(n.cells=n()) %>%
  group_by(datasets) %>%
  mutate(Clust.Perc = round((n.cells / sum(n.cells)*100),0))

#Count cells by other cluster level

a<-ggplot(table1, aes(y=reorder(as.factor(cca_clusters),n.cells), x=n.cells)) + geom_col(width = 0.9,position="dodge") + 
  ggtitle("Number of cells by Type") + 
  theme(axis.text.x=element_text(angle = 90),
        axis.text.y=element_text(size=6)) + ylab("")+ NoLegend()

b<-ggplot(table, aes(reorder(as.factor(cca_clusters),n.cells), x=n.cells, fill=datasets)) + geom_col(width = 0.9,position="fill" ) + 
  ggtitle("Datasets percentage in each Cluster") + theme(axis.text.x=element_text(angle = 90),
                                                   axis.text.y=element_text(size=8),
                                                   legend.position="right", 
                                                   legend.key.size = unit(0.4,"cm"),
                                                   legend.text = element_text(size=8)
  )+ylab("Clusters") + xlab("Percentage")
library(ggpubr)
ggarrange(a,b)
## Clustering Analysis with Seurat ----
all_astro_merge <- FindNeighbors(all_astro_merge, reduction = "integrated.cca", dims = 1:30)
all_astro_merge <- FindClusters(all_astro_merge, resolution = 0.5, cluster.name = "cca_clusters")
Idents(all_astro_merge) <- 'cca_clusters'
markers_cca <- FindAllMarkers(all_astro_merge, only.pos = TRUE, logfc.threshold = 0.1)
markers_cca <-markers_cca[order(markers_cca$avg_log2FC, decreasing=TRUE),]
markers_cca <- split(markers_cca, f = markers_cca$cluster)

all_astro_merge <- RunUMAP(all_astro_merge, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

save.image("markers_all_astrocytes_integrated.RDS")

p1 <- DimPlot(
  all_astro_merge,
  reduction = "umap.cca",
  group.by = c("datasets"),
  combine = FALSE)

p2 <- DimPlot(
  all_astro_merge,
  reduction = "umap.cca",
  group.by = c("sex"),
  combine = FALSE)


pdf(file="pdf/integrated_cca_datasets.pdf")
p1
p2
dev.off()

# Write xlsx files for the markers
library(writexl)
write_xlsx(markers_cca[[1]],"/mnt/remo/prj/transcriptome_cerebellum/dataset/3fold_integ/ccamarkers0.xlsx")
write_xlsx(markers_cca[[2]],"/mnt/remo/prj/transcriptome_cerebellum/dataset/3fold_integ/ccamarkers1.xlsx")
write_xlsx(markers_cca[[3]],"/mnt/remo/prj/transcriptome_cerebellum/dataset/3fold_integ/ccamarkers2.xlsx")
write_xlsx(markers_cca[[4]],"/mnt/remo/prj/transcriptome_cerebellum/dataset/3fold_integ/ccamarkers3.xlsx")

## Compare clusters markers with markers identified in the developmental dataset ----
Markers_all <- read.csv("Astromarkers_AstroAll.txt")
CNA <- Markers_all[Markers_all$cluster==9,]
CNA <- CNA[CNA$avg_log2FC > 0.5,]
CNA <- CNA[CNA$pct.1 > 0.5,]
CNA <- CNA$gene

BG <- Markers_all[Markers_all$cluster==7,]
BG <- BG[BG$avg_log2FC > 0.5,]
BG <- BG[BG$pct.1 > 0.5,]
BG <- BG$gene

GLA <- Markers_all[Markers_all$cluster==6,]
GLA <- GLA[GLA$avg_log2FC > 0.5,]
GLA <- GLA[GLA$pct.1 > 0.5,]
GLA <- GLA$gene

WMA <- Markers_all[Markers_all$cluster==5,]
WMA <- WMA[WMA$avg_log2FC > 0.5,]
WMA <- WMA[WMA$pct.1 > 0.5,]
WMA <- WMA$gene

library(readxl)

Schist_markers0 <- read_excel("schistmarkers0.xlsx")
Schist_markers0 <-Schist_markers0[Schist_markers0$avg_log2FC > 0.5,]
Schist_markers0 <- Schist_markers0[Schist_markers0$pct.1 > 0.25,]
Schist_markers0 <- Schist_markers0$gene

Schist_markers1 <- read_excel("schistmarkers1.xlsx")
Schist_markers1 <- Schist_markers1[Schist_markers1$avg_log2FC > 0.5,]
Schist_markers1 <- Schist_markers1[Schist_markers1$pct.1 > 0.25,]
Schist_markers1 <- Schist_markers1$gene

Schist_markers2 <- read_excel("schistmarkers2.xlsx")
Schist_markers2 <- Schist_markers2[Schist_markers2$avg_log2FC > 0.5,]
Schist_markers2 <- Schist_markers2[Schist_markers2$pct.1 > 0.25,]
Schist_markers2 <- Schist_markers2$gene

Schist_markers3 <- read_excel("schistmarkers3.xlsx")
Schist_markers3 <- Schist_markers3[Schist_markers3$avg_log2FC > 0.5,]
Schist_markers3 <- Schist_markers3[Schist_markers3$pct.1 > 0.25,]
Schist_markers3 <- Schist_markers3$gene

Schist_markers4 <- read_excel("schistmarkers4.xlsx")
Schist_markers4 <- Schist_markers4[Schist_markers4$avg_log2FC > 0.5,]
Schist_markers4 <- Schist_markers4[Schist_markers4$pct.1 > 0.25,]
Schist_markers4 <- Schist_markers4$gene

Cca_markers0 <- read_excel("ccamarkers0.xlsx")
Cca_markers0 <- Cca_markers0[Cca_markers0$avg_log2FC > 0.5,]
Cca_markers0 <- Cca_markers0[Cca_markers0$pct.1 > 0.5,]
Cca_markers0 <- Cca_markers0$gene

Cca_markers1 <- read_excel("ccamarkers1.xlsx")
Cca_markers1 <- Cca_markers1[Cca_markers1$avg_log2FC > 0.5,]
Cca_markers1 <- Cca_markers1[Cca_markers1$pct.1 > 0.5,]
Cca_markers1 <- Cca_markers1$gene

Cca_markers2 <- read_excel("ccamarkers2.xlsx")
Cca_markers2 <- Cca_markers2[Cca_markers2$avg_log2FC > 0.5,]
Cca_markers2 <- Cca_markers2[Cca_markers2$pct.1 > 0.5,]
Cca_markers2 <- Cca_markers2$gene

Cca_markers3 <- read_excel("ccamarkers3.xlsx")
Cca_markers3 <- Cca_markers3[Cca_markers3$avg_log2FC > 0.5,]
Cca_markers3 <- Cca_markers3[Cca_markers3$pct.1 > 0.5,]
Cca_markers3 <- Cca_markers3$gene

library(ComplexHeatmap)
list <- list(Schist_markers0, Schist_markers1, Schist_markers2, Schist_markers3, Schist_markers4, BG, GLA, WMA, CNA)
names(list) <- c("Cluster0", "Cluster1", "Cluster2", "Cluster3", "Cluster4" ,"BG", "GLA", "WMA", "CNA")
comb_matr_intersect <- make_comb_mat(list, mode = "intersect")
comb_matr_distinct <- make_comb_mat(list, mode = "distinct")

UpSet(comb_matr_intersect[c(29:32,34:43)], set_order = c("BG", "GLA", "WMA", "CNA", "Cluster0", "Cluster1", "Cluster2", "Cluster3", "Cluster4"), comb_order = order(comb_size(comb_matr_intersect[c(29:32,34:43)]), decreasing = TRUE))
UpSet(comb_matr_distinct[c(43:46, 49:50)], set_order = c("BG", "GLA", "WMA", "CNA", "Cluster0", "Cluster1", "Cluster2", "Cluster3", "Cluster4"), comb_order = order(comb_size(comb_matr_distinct[c(43:46, 49:50)]), decreasing = TRUE))

list <- list(Cca_markers0, Cca_markers1, Cca_markers2, Cca_markers3, BG, GLA, WMA, CNA)
names(list) <- c("Cluster0", "Cluster1", "Cluster2", "Cluster3" ,"BG", "GLA", "WMA", "CNA")
comb_matr_intersect_cca <- make_comb_mat(list, mode = "intersect")
comb_matr_distinct_cca <- make_comb_mat(list, mode = "distinct")

UpSet(comb_matr_intersect_cca[c(29:32,34:43)], set_order = c("BG", "GLA", "WMA", "CNA", "Cluster0", "Cluster1", "Cluster2", "Cluster3"), comb_order = order(comb_size(comb_matr_intersect[c(29:32,34:43)]), decreasing = TRUE))
UpSet(comb_matr_distinct_cca[c(27:29, 31:32, 34:39)], set_order = c("BG", "GLA", "WMA", "CNA", "Cluster0", "Cluster1", "Cluster2", "Cluster3"), comb_order = order(comb_size(comb_matr_distinct[c(27:29, 31:32, 34:39)]), decreasing = TRUE))

# Check the expression of known markers to classify the clusters ----
Idents(all_astro_merge) <- all_astro_merge$datasets

library(ggplot2)
library(patchwork)
library(ggpubr)
##GAT3 ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Slc6a11", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Slc6a11 - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Slc6a11", group.by = "schist_clusters"
) + NoLegend() + ggtitle("Slc6a11 - schist Clusters")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_schist_GAT3.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p3 <- DimPlot(
  all_astro_merge,
  reduction = "umap.cca",
  group.by = c("cca_clusters"),
  combine = FALSE, label = T)

p4 <- DimPlot(
  all_astro_merge,
  reduction = "umap.cca",
  group.by = c("clusters_for_GO"),
  combine = FALSE, label = T)

Dimplots <- wrap_plots(c(p3, p4), ncol = 2)

p5 <- FeaturePlot(all_astro_merge, features = "Slc6a11", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")

Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_schist_GAT3.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Slc6a11", split.by = "cca_clusters", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitCcaClusters_GAT3.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

p9 <- FeaturePlot(all_astro_merge, features = "Slc6a11", split.by = "datasets", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitStudy_GAT3.pdf", plot =p9, device =cairo_pdf, width = 297, height = 210, units = "mm")

##AQP4 ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Aqp4", group.by = "cca_clusters"
) + NoLegend() + ggtitle("AQP4 - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Aqp4", group.by = "schist_clusters"
) + NoLegend() + ggtitle("AQP4 - schist Clusters")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_schist_AQP4.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Aqp4", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")

Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_schist_Aqp4.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Aqp4", split.by = "cca_clusters", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitCcaClusters_Aqp4.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")


##Gdf10 ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Gdf10", group.by = "cca_clusters"
) + NoLegend() + ggtitle("GDF10 - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Gdf10", group.by = "schist_clusters"
) + NoLegend() + ggtitle("GDF10 - schist Clusters")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_schist_GDF10.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Gdf10", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")

Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_schist_Gdf10.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Gdf10", split.by = "cca_clusters", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitCcaClusters_Gdf10.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

##Fam107a ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Fam107a", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Fam107A - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Fam107a", group.by = "schist_clusters"
) + NoLegend() + ggtitle("Fam107a - schist Clusters")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_schist_Fam107a.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Fam107a", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")

Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_schist_Fam107a.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Fam107a", split.by = "cca_clusters", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitCcaClusters_Fam107a.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")


##Slc1a3 ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Slc1a3", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Slc1a3 - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Slc1a3", group.by = "schist_clusters"
) + NoLegend() + ggtitle("Slc1a3 - schist Clusters")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_schist_Glast.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Slc1a3", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")

Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_schist_Slc1a3.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Slc1a3", split.by = "cca_clusters", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitCcaClusters_Slc1a3.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

##Gria1 ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Gria1", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Gria1 - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Gria1", group.by = "schist_clusters"
) + NoLegend() + ggtitle("Gria1 - schist Clusters")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_schist_Gria1.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Gria1", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")

Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_seurat_Gria1.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Gria1", split.by = "cca_clusters", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitCcaClusters_Gria1.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

##Gria4 ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Gria4", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Gria4 - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Gria4", group.by = "schist_clusters"
) + NoLegend() + ggtitle("Gria4 - schist Clusters")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_schist_Gria4.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Gria4", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")

Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_schist_Gria4_.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Gria4", split.by = "cca_clusters", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitCcaClusters_Gria4.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

##Myoc ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Myoc", group.by = "cca_clusters"
) + NoLegend() + ggtitle("GLAST - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Myoc", group.by = "mnn_clusters"
) + NoLegend() + ggtitle("Myoc - mnn Clusters")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_mnn_Myoc.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Myoc", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")
p6 <- FeaturePlot(all_astro_merge, features = "Myoc", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.mnn")

Featureplots <- p5 | p6

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_mnn_Myoc.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Myoc", split.by = "cca_clusters", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitCcaClusters_Myoc.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

## Slc1a2 ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Slc1a2", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Slc1a2 - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Slc1a2", group.by = "schist_clusters"
) + NoLegend() + ggtitle("Myoc - schist Clusters")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_schist_Slc1a2.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Slc1a2", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")
p6 <- FeaturePlot(all_astro_merge, features = "Slc1a2", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.mnn")

Featureplots <- p5 | p6

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_schist_Slc1a2.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Slc1a2", split.by = "cca_clusters", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitCcaClusters_Slc1a2.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

##Tekt5 ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Tekt5", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Tekt5 - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Tekt5", group.by = "schist_clusters"
) + NoLegend() + ggtitle("Tekt5 - schist Clusters")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_schist_Tekt5.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Tekt5", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")
p6 <- FeaturePlot(all_astro_merge, features = "Tekt5", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.mnn")

Featureplots <- p5 | p6

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_schist_Tekt5.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Tekt5", split.by = "cca_clusters", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitCcaClusters_Tekt5.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

##Gabra6 ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Gabra6", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Gabra6 - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Gabra6", group.by = "clusters_for_GO"
) + NoLegend() + ggtitle("Gabra6 - manual clusters ")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_manual_Gabra6.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Gabra6", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")
Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_manual_Gabra6.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Gabra6", split.by = "clusters_for_GO", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitmanualClusters_Gabra6.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

##Mctp1 ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Mctp1", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Mctp1 - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Mctp1", group.by = "clusters_for_GO"
) + NoLegend() + ggtitle("Mctp1 - manual clusters ")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_manual_Mctp1.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Mctp1", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")
Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_manual_Mctp1.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Mctp1", split.by = "clusters_for_GO", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitmanualClusters_Mctp1.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

##Olfm3 ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Olfm3", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Olfm3 - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Olfm3", group.by = "clusters_for_GO"
) + NoLegend() + ggtitle("Olfm3 - manual clusters ")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_manual_Olfm3.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Olfm3", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")
Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_manual_Olfm3.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Olfm3", split.by = "clusters_for_GO", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitmanualClusters_Olfm3.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

##Tenm1 ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Tenm1", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Tenm1 - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Tenm1", group.by = "clusters_for_GO"
) + NoLegend() + ggtitle("Tenm1 - manual clusters ")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_manual_Tenm1.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Tenm1", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")
Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_manual_Tenm1.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Tenm1", split.by = "clusters_for_GO", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitmanualClusters_Tenm1.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

##Pde10a ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Pde10a", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Pde10a - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Pde10a", group.by = "clusters_for_GO"
) + NoLegend() + ggtitle("Pde10a - manual clusters ")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_manual_Pde10a.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Pde10a", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")
Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_manual_Pde10a.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Pde10a", split.by = "clusters_for_GO", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitmanualClusters_Pde10a.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

##Ablim1 ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Ablim1", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Ablim1 - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Ablim1", group.by = "clusters_for_GO"
) + NoLegend() + ggtitle("Ablim1 - manual clusters ")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_manual_Ablim1.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Ablim1", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")
Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_manual_Ablim1.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Ablim1", split.by = "clusters_for_GO", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitmanualClusters_Ablim1.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

##Fabp7 ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Fabp7", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Fabp7 - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Fabp7", group.by = "clusters_for_GO"
) + NoLegend() + ggtitle("Fabp7 - manual clusters ")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_manual_Fabp7.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Fabp7", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")
Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_manual_Fabp7.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Fabp7", split.by = "clusters_for_GO", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitmanualClusters_Fabp7.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

##Hopx ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Hopx", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Hopx - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Hopx", group.by = "clusters_for_GO"
) + NoLegend() + ggtitle("Hopx - manual clusters ")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_manual_Hopx.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Hopx", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")
Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_manual_Hopx.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Hopx", split.by = "clusters_for_GO", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitmanualClusters_Hopx.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

##Vim ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Vim", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Vim - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Vim", group.by = "clusters_for_GO"
) + NoLegend() + ggtitle("Vim - manual clusters ")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_manual_Vim.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Vim", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")
Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_manual_Vim.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Vim", split.by = "clusters_for_GO", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitmanualClusters_Vim.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")

##Snap25 ----
p1 <- VlnPlot(
  all_astro_merge,
  features = "Snap25", group.by = "cca_clusters"
) + NoLegend() + ggtitle("Snap25 - cca Clusters")

p2 <- VlnPlot(
  all_astro_merge,
  features = "Snap25", group.by = "clusters_for_GO"
) + NoLegend() + ggtitle("Snap25 - manual clusters ")

p3 <- p1 | p2
ggsave(filename = "pdf/VlnPlot_cca_vs_manual_Snap25.pdf", plot =p3, device =cairo_pdf, width = 297, height = 210, units = "mm")

p5 <- FeaturePlot(all_astro_merge, features = "Snap25", pt.size = 1,max.cutoff="q99",min.cutoff="q1", reduction = "umap.cca")
Featureplots <- p5

p7 <- ggarrange(Dimplots, Featureplots, ncol = 1, nrow = 2)
ggsave(filename = "pdf/Featureplot_cca_vs_manual_Snap25.pdf", plot =p7, device =cairo_pdf, width = 297, height = 210, units = "mm")

p8 <- FeaturePlot(all_astro_merge, features = "Snap25", split.by = "clusters_for_GO", reduction = "umap.cca",pt.size = 1,max.cutoff="q99",min.cutoff="q1")
ggsave(filename = "pdf/Featureplot_splitmanualClusters_Snap25.pdf", plot =p8, device =cairo_pdf, width = 297, height = 210, units = "mm")


# Plot top20 markers for each cluster on Visium dataset to help classifying the clusters ----
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

all_astro_merge@active.ident<-all_astro_merge$cca_clusters
integrated_cca <- all_astro_merge
saveRDS(integrated_cca, file = "all_astrocytes_integrated_cca.Rds")

##Adult Visium dataset
InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "posterior1")
brain2 <- LoadData("stxBrain", type = "posterior2")
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)
brain2 <- RunPCA(object = brain2,verbose = FALSE)
brain2 <- RunUMAP(object = brain2, dims = 1:30)
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
brain <- RunPCA(object=brain, verbose = FALSE)
brain <- RunUMAP(object=brain, dims = 1:30)

DefaultAssay(brain2) <- "SCT"
DefaultAssay(brain) <- "SCT"

crb <- subset(brain, posterior1_imagerow < 370 & posterior1_imagecol > 310, invert = FALSE)

## Select the markers
#cca_mark_select0 <- markers_cca[[1]][markers_cca[[1]]$avg_log2FC>0.5 & markers_cca[[1]]$pct.1 - markers_cca[[1]]$pct.2 > 0.3 ,]$gene
#cca_mark_select1 <- markers_cca[[2]][markers_cca[[2]]$avg_log2FC>0.5 & markers_cca[[2]]$pct.1 - markers_cca[[2]]$pct.2 > 0.3 ,]$gene
#cca_mark_select2 <- markers_cca[[3]][markers_cca[[3]]$avg_log2FC>0.5 & markers_cca[[3]]$pct.1 - markers_cca[[3]]$pct.2 > 0.3 ,]$gene
#cca_mark_select3 <- markers_cca[[4]][markers_cca[[4]]$avg_log2FC>0.5 & markers_cca[[4]]$pct.1 - markers_cca[[4]]$pct.2 > 0.3 ,]$gene
cca_mark_select0 <- markers_cca[[1]][1:50,]$gene
cca_mark_select1 <- markers_cca[[2]][1:50,]$gene
cca_mark_select2 <- markers_cca[[3]][1:50,]$gene
cca_mark_select3 <- markers_cca[[4]][1:50,]$gene

#schist_mark_select0 <- markers_schist[[1]][markers_schist[[1]]$avg_log2FC>0.5 & markers_schist[[1]]$pct.1 - markers_schist[[1]]$pct.2 > 0.3 ,]$gene
#schist_mark_select1 <- markers_schist[[2]][markers_schist[[2]]$avg_log2FC>0.5 & markers_schist[[2]]$pct.1 - markers_schist[[2]]$pct.2 > 0.3 ,]$gene
#schist_mark_select2 <- markers_schist[[3]][markers_schist[[3]]$avg_log2FC>0.5 & markers_schist[[3]]$pct.1 - markers_schist[[3]]$pct.2 > 0.3 ,]$gene
#schist_mark_select3 <- markers_schist[[4]][markers_schist[[4]]$avg_log2FC>0.5 & markers_schist[[4]]$pct.1 - markers_schist[[4]]$pct.2 > 0.3 ,]$gene
#schist_mark_select4 <- markers_schist[[5]][markers_schist[[5]]$avg_log2FC>0.5 & markers_schist[[5]]$pct.1 - markers_schist[[5]]$pct.2 > 0.3 ,]$gene
schist_mark_select0 <- markers_schist[[1]][1:20,]$gene
schist_mark_select1 <- markers_schist[[2]][1:20,]$gene
schist_mark_select2 <- markers_schist[[3]][1:20,]$gene
schist_mark_select3 <- markers_schist[[4]][1:20,]$gene
schist_mark_select4 <- markers_schist[[5]][1:20,]$gene

genes.list_cca <- list(cca_mark_select0, cca_mark_select1, cca_mark_select2, cca_mark_select3)
enrich.name <- "50cca_clusters_in_visium"
crb <- AddModuleScore(crb,
                                 features = genes.list_cca,
                                 pool = NULL,
                                 n.bin = 5,
                                 seed = 1,
                                 ctrl = length(genes.list_cca),
                                 k = FALSE,
                                 name = enrich.name,
                                 random.seed = 1)

pdf("pdf/50cca_markers_in_Visiumbrain.pdf")
SpatialFeaturePlot(crb, features = "50cca_clusters_in_visium1", pt.size.factor = 1)
SpatialFeaturePlot(crb, features = "50cca_clusters_in_visium2", pt.size.factor = 1)
SpatialFeaturePlot(crb, features = "50cca_clusters_in_visium3", pt.size.factor = 1)
SpatialFeaturePlot(crb, features = "50cca_clusters_in_visium4", pt.size.factor = 1)
dev.off()

genes.list_schist <- list(schist_mark_select0, schist_mark_select1, schist_mark_select2, schist_mark_select3, schist_mark_select4)
enrich.name <- "20schist_clusters_in_visium"
brain <- AddModuleScore(brain,
                        features = genes.list_schist,
                        pool = NULL,
                        n.bin = 5,
                        seed = 1,
                        ctrl = length(genes.list_schist),
                        k = FALSE,
                        name = enrich.name,
                        random.seed = 1)

pdf("pdf/20schist_markers_in_Visiumbrain.pdf")
SpatialFeaturePlot(brain, features = "20schist_clusters_in_visium1", pt.size.factor = 1)
SpatialFeaturePlot(brain, features = "20schist_clusters_in_visium2", pt.size.factor = 1)
SpatialFeaturePlot(brain, features = "20schist_clusters_in_visium3", pt.size.factor = 1)
SpatialFeaturePlot(brain, features = "20schist_clusters_in_visium4", pt.size.factor = 1)
SpatialFeaturePlot(brain, features = "20schist_clusters_in_visium5", pt.size.factor = 1)
dev.off()

## Integrate the datasets ----
all_astro_merge <- readRDS("all_astrocytes_integrated.Rds")

all_astro_merge@active.ident<-all_astro_merge$cca_clusters

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "posterior1")
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
brain <- RunPCA(object=brain, verbose = FALSE)
brain <- RunUMAP(object=brain, dims = 1:30)
DefaultAssay(brain) <- "SCT"

crb <- subset(brain, posterior1_imagerow < 370 & posterior1_imagecol > 310, invert = FALSE)
crb <- SCTransform(crb, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

Idents(all_astro_merge) <- all_astro_merge$cca_clusters
DefaultAssay(all_astro_merge) <- "SCT"
anchors <- FindTransferAnchors(reference = all_astro_merge, query = crb, normalization.method = "SCT", 
                               recompute.residuals = F)

predictions.assay <- TransferData(anchorset = anchors, refdata = all_merge$cca_clusters, prediction.assay = TRUE,
                                  weight.reduction = crb[["pca"]], dims = 1:30)
crb[["predictions"]] <- predictions.assay

DefaultAssay(crb) <- "predictions"

pdf("E:/Giacomo/Cerebellum_adult_mouse_single_cell/Plots/cca_clusters_in_visium.pdf")
SpatialFeaturePlot(crb, features = "0", pt.size.factor = 1, ncol = 1, crop = TRUE)
SpatialFeaturePlot(crb, features = "1", pt.size.factor = 1, ncol = 1, crop = TRUE)
SpatialFeaturePlot(crb, features = "2", pt.size.factor = 1, ncol = 1, crop = TRUE)
SpatialFeaturePlot(crb, features = "3", pt.size.factor = 1, ncol = 1, crop = TRUE)
dev.off()


# Plot cca_clusters divided for dataset
library(ggplot2)
library(dplyr)

clusters <- Idents(all_astro_merge)
studies <- all_astro_merge$datasets
dataframe <- data.frame(clusters, studies)
cluster_freq <- dataframe %>% 
  group_by(studies, clusters) %>% 
  summarise(frequency= n()) %>% 
  ungroup() %>% 
  mutate(total_cells= ave(frequency, clusters, FUN = sum)) %>%
  mutate(relative_frequency = frequency/total_cells)

plot <- ggplot(cluster_freq, aes(x = clusters, y = relative_frequency, fill = studies)) +
  geom_bar(stat = "identity") +
  #scale_fill_manual(values = c("Group 1" = "blue", "Group 2" = "red", "Group 3" = "green")) +  # Customize fill colors as desired
  labs(x = "Cluster", y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "right") +
  #scale_x_discrete(limits = cluster_order) +
  scale_y_continuous(labels = scales::comma)

ggsave(filename = "pdf/Frequency_datasets_across_clusters.pdf", plot =plot, device =cairo_pdf, width = 297, height = 210, units = "mm")



# DEG analysis ----

load("Adult_annotated_Markers.Rds") #load Astro_adult and other objects from the processing in "integration_development_adult.R"

## obtain vectors of the markers of adult clusters ----
pGLA <- Markers_presumed[[1]]
pGLA <- pGLA[pGLA$avg_log2FC > 0.5,]
pGLA <- pGLA[pGLA$pct.1 > 0.5,]
pGLA <- pGLA$gene

pCNA <- Markers_presumed[[2]]
pCNA <- pCNA[pCNA$avg_log2FC > 0.5,]
pCNA <- pCNA[pCNA$pct.1 > 0.5,]
pCNA <- pCNA$gene

pGLA2 <- Markers_presumed[[3]]
pGLA2 <- pGLA2[pGLA2$avg_log2FC > 0.5,]
pGLA2 <- pGLA2[pGLA2$pct.1 > 0.5,]
pGLA2 <- pGLA2$gene

pBG <- Markers_presumed[[4]]
pBG <- pBG[pBG$avg_log2FC > 0.5,]
pBG <- pBG[pBG$pct.1 > 0.5,]
pBG <- pBG$gene

pImmature_like <- Markers_presumed[[5]]
pImmature_like <- pImmature_like[pImmature_like$avg_log2FC > 0.5,]
pImmature_like <- pImmature_like[pImmature_like$pct.1 > 0.5,]
pImmature_like <- pImmature_like$gene

## Look at the intersection between developmental markers and adult ones ----
library(ComplexHeatmap)
list <- list(pGLA, pCNA, pGLA2, pBG, pImmature_like, BG, GLA, WMA, CNA)
names(list) <- c("presumed GLA", "presumed CNA", "Presumed GLA2", "Presumed BG", "Presumed Immature Like" ,"BG", "GLA", "WMA", "CNA")
comb_matr_intersect <- make_comb_mat(list, mode = "intersect")
comb_matr_distinct <- make_comb_mat(list, mode = "distinct")

pdf("pdf/markers_compared_adult_devel.pdf")
UpSet(comb_matr_distinct[c(35:38, 41:45)], set_order = c("BG", "GLA", "WMA", "CNA", "presumed GLA", "presumed CNA", "Presumed WMA", "Presumed BG", "Presumed Immature"), comb_order = order(comb_size(comb_matr_distinct[c(35:38, 41:45)]), decreasing = TRUE))
UpSet(comb_matr_intersect[c(4:5,20:29,33:36,41:42,44,55:77)], set_order = c("BG", "GLA", "WMA", "CNA", "presumed GLA", "presumed CNA", "Presumed WMA", "Presumed BG", "Presumed Immature"), comb_order = order(comb_size(comb_matr_intersect[c(4:5,20:29,33:36,41:42,44,55:77)]), decreasing = TRUE))
dev.off()

## Select genes expressed by all the adult types

list <- list(pGLA, pCNA, pWMA, pBG, pImmature)
names(list) <- c("presumed GLA", "presumed CNA", "Presumed WMA", "Presumed BG", "Presumed Immature")
comb_matr_intersect <- make_comb_mat(list, mode = "intersect")

#GO analysis ----
library(topGO)
library(ViSEAGO)
library(clusterProfiler, verbose = FALSE)
library(org.Mm.eg.db, verbose = FALSE)

BG.entrez1 <- clusterProfiler::bitr(pBG,fromType = "SYMBOL",toType = "ENTREZID", OrgDb = org.Mm.eg.db)
BG.entrez1 <- BG.entrez1$ENTREZID

GLA.entrez1 <- clusterProfiler::bitr(pGLA,fromType = "SYMBOL",toType = "ENTREZID", OrgDb = org.Mm.eg.db)
GLA.entrez1 <- GLA.entrez1$ENTREZID

WMA.entrez1 <- clusterProfiler::bitr(pWMA,fromType = "SYMBOL",toType = "ENTREZID", OrgDb = org.Mm.eg.db)
WMA.entrez1 <- WMA.entrez1$ENTREZID

CNA.entrez1 <- clusterProfiler::bitr(pCNA,fromType = "SYMBOL",toType = "ENTREZID", OrgDb = org.Mm.eg.db)
CNA.entrez1 <- CNA.entrez1$ENTREZID

Immature.entrez1 <- clusterProfiler::bitr(pImmature,fromType = "SYMBOL",toType = "ENTREZID", OrgDb = org.Mm.eg.db)
Immature.entrez1 <- Immature.entrez1$ENTREZID

background <- rownames(Astro_adult)
bkgd.genes.entrez <- clusterProfiler::bitr(background,fromType = "SYMBOL",toType = "ENTREZID", OrgDb = org.Mm.eg.db)
bkgd.genes.entrez <- bkgd.genes.entrez$ENTREZID

## connect to Bioconductor
Bioconductor<-ViSEAGO::Bioconductor2GO()

## load GO annotations from Bioconductor
myGENE2GO <- ViSEAGO::annotate(
  "org.Mm.eg.db",
  Bioconductor
)

## create topGOdata for BP for each list of DE genes
BP_BG <-ViSEAGO::create_topGOdata(
  geneSel= BG.entrez1,
  allGenes= bkgd.genes.entrez,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

BP_GLA <-ViSEAGO::create_topGOdata(
  geneSel= GLA.entrez1,
  allGenes= bkgd.genes.entrez,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

BP_CNA <-ViSEAGO::create_topGOdata(
  geneSel= CNA.entrez1,
  allGenes= bkgd.genes.entrez,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

BP_WMA <-ViSEAGO::create_topGOdata(
  geneSel= WMA.entrez1,
  allGenes= bkgd.genes.entrez,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

BP_Immature <-ViSEAGO::create_topGOdata(
  geneSel= Immature.entrez1,
  allGenes= bkgd.genes.entrez,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

## perform topGO tests
classic_BP_BG<-topGO::runTest(
  BP_BG,
  algorithm ="classic",
  statistic = "fisher",
  cutOff=0.05
)

classic_BP_GLA<-topGO::runTest(
  BP_GLA,
  algorithm ="classic",
  statistic = "fisher",
  cutOff=0.05
)

classic_BP_CNA <-topGO::runTest(
  BP_CNA,
  algorithm ="classic",
  statistic = "fisher",
  cutOff=0.05
)

classic_BP_WMA<-topGO::runTest(
  BP_WMA,
  algorithm ="classic",
  statistic = "fisher",
  cutOff=0.05
)

classic_BP_Immature<-topGO::runTest(
  BP_Immature,
  algorithm ="classic",
  statistic = "fisher",
  cutOff=0.05
)

# merge topGO results, each cell type separately
BP_sResults_BG <-ViSEAGO::merge_enrich_terms(
  cutoff = 0.05,
  Input=list(
    condition=c("BP_BG",
                "classic_BP_BG")
  )
)

BP_sResults_GLA <-ViSEAGO::merge_enrich_terms(
  cutoff = 0.05,
  Input=list(
    condition=c("BP_GLA",
                "classic_BP_GLA")
  )
)

BP_sResults_WMA <-ViSEAGO::merge_enrich_terms(
  cutoff = 0.05,
  Input=list(
    condition=c("BP_WMA",
                "classic_BP_WMA")
  )
)

BP_sResults_CNA <-ViSEAGO::merge_enrich_terms(
  cutoff = 0.05,
  Input=list(
    condition=c("BP_CNA",
                "classic_BP_CNA")
  )
)

BP_sResults_Immature <-ViSEAGO::merge_enrich_terms(
  cutoff = 0.05,
  Input=list(
    condition=c("BP_Immature",
                "classic_BP_Immature")
  )
)

# merge topGO results, all types together
BP_sResults<-ViSEAGO::merge_enrich_terms(
  cutoff=0.05,
  Input=list(
    BG=c(
      "BP_BG",
      "classic_BP_BG"
    ),
    GLA=c(
      "BP_GLA",
      "classic_BP_GLA"
    ),
    CNA=c(
      "BP_CNA",
      "classic_BP_CNA"
    ),
    WMA=c(
      "BP_WMA",
      "classic_BP_WMA"
    ),
    Immature=c(
      "BP_Immature",
      "classic_BP_Immature"
    )
  )
)

# display a summary
BP_sResults_BG
BP_sResults_GLA
BP_sResults_WMA
BP_sResults_CNA
BP_sResults_Immature
BP_sResults

#save outputs
ViSEAGO::show_table(
  BP_sResults_BG,
  "out/ViseaGO/BP_sResults_BG_sole.xls"
)

ViSEAGO::show_table(
  BP_sResults_GLA,
  "out/ViseaGO/BP_sResults_GLA_sole.xls"
)

ViSEAGO::show_table(
  BP_sResults_WMA,
  "out/ViseaGO/BP_sResults_WMA_sole.xls"
)

ViSEAGO::show_table(
  BP_sResults_CNA,
  "out/ViseaGO/BP_sResults_CNA_sole.xls"
)


ViSEAGO::show_table(
  BP_sResults_Immature,
  "out/ViseaGO/BP_sResults_Immature_sole.xls"
)

ViSEAGO::show_table(
  BP_sResults,
  "out/ViseaGO/BP_sResults.xls"
)

# barchart of significant (or not) GO terms by comparison
ViSEAGO::GOcount(BP_sResults)

# display intersections
ViSEAGO::Upset(
  BP_sResults,
  file="out/ViseaGO/upset.csv"
)

# Initialize myGOs

myGOs_BG<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_sResults_BG
)

myGOs_CNA<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_sResults_CNA
)

myGOs_GLA<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_sResults_GLA
)

myGOs_WMA<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_sResults_WMA
)

myGOs_Immature<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_sResults_Immature
)

myGOs<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_sResults
)

# Compute Semantic Similarity

myGOs_BG<-ViSEAGO::compute_SS_distances(
  myGOs_BG,
  distance="Wang"
)

myGOs_CNA<-ViSEAGO::compute_SS_distances(
  myGOs_CNA,
  distance="Wang"
)

myGOs_GLA<-ViSEAGO::compute_SS_distances(
  myGOs_GLA,
  distance="Wang"
)

myGOs_WMA<-ViSEAGO::compute_SS_distances(
  myGOs_WMA,
  distance="Wang"
)

myGOs_Immature<-ViSEAGO::compute_SS_distances(
  myGOs_Immature,
  distance="Wang"
)

myGOs<-ViSEAGO::compute_SS_distances(
  myGOs,
  distance="Wang"
)
# Create GOterms heatmap
Wang_clusters_wardD2_BG<-ViSEAGO::GOterms_heatmap(
  myGOs_BG,
  showIC=FALSE,
  showGOlabels = FALSE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=2,
        minClusterSize =5
      )
    )
  ),
  samples.tree=NULL
)

Wang_clusters_wardD2_GLA<-ViSEAGO::GOterms_heatmap(
  myGOs_GLA,
  showIC=FALSE,
  showGOlabels = FALSE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=2,
        minClusterSize =5
      )
    )
  ),
  samples.tree=NULL
)

Wang_clusters_wardD2_WMA<-ViSEAGO::GOterms_heatmap(
  myGOs_WMA,
  showIC=FALSE,
  showGOlabels = FALSE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=2,
        minClusterSize =5
      )
    )
  ),
  samples.tree=NULL
)

Wang_clusters_wardD2_CNA<-ViSEAGO::GOterms_heatmap(
  myGOs_CNA,
  showIC=FALSE,
  showGOlabels = FALSE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=2,
        minClusterSize =5
      )
    )
  ),
  samples.tree=NULL
)

Wang_clusters_wardD2_Immature<-ViSEAGO::GOterms_heatmap(
  myGOs_Immature,
  showIC=FALSE,
  showGOlabels =FALSE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=2,
        minClusterSize =5
      )
    )
  ),
  samples.tree=NULL
)


Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
  myGOs,
  showIC=FALSE,
  showGOlabels =FALSE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=2,
        minClusterSize =5
      )
    )
  ),
  samples.tree=NULL
)

save.image("ViSEAGO_heatmaps_generated.Rds")

# Print the heatmaps

ViSEAGO::show_heatmap(
  Wang_clusters_wardD2_Immature,
  "GOterms"
)

ViSEAGO::show_heatmap(
  Wang_clusters_wardD2_BG,
  "GOterms"
)

ViSEAGO::show_heatmap(
  Wang_clusters_wardD2_WMA,
  "GOterms"
)

ViSEAGO::show_heatmap(
  Wang_clusters_wardD2_GLA,
  "GOterms"
)

ViSEAGO::show_heatmap(
  Wang_clusters_wardD2_CNA,
  "GOterms"
)

ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms"
)

ViSEAGO::MDSplot(
  Wang_clusters_wardD2_GLA,
  "GOterms"
)

ViSEAGO::MDSplot(
  Wang_clusters_wardD2_WMA,
  "GOterms"
)

ViSEAGO::MDSplot(
  Wang_clusters_wardD2_CNA,
  "GOterms"
)

ViSEAGO::MDSplot(
  Wang_clusters_wardD2_BG,
  "GOterms"
)

ViSEAGO::MDSplot(
  Wang_clusters_wardD2,
  "GOterms"
)

ViSEAGO::show_table(
  Wang_clusters_wardD2_CNA,
  "out/ViseaGO/BP_sResults_clusters_CNA_sole.csv"
)

ViSEAGO::show_table(
  Wang_clusters_wardD2_BG,
  "out/ViseaGO/BP_sResults_clusters_BG_sole.csv"
)

ViSEAGO::show_table(
  Wang_clusters_wardD2_GLA,
  "out/ViseaGO/BP_sResults_clusters_GLA_sole.csv"
)

ViSEAGO::show_table(
  Wang_clusters_wardD2_WMA,
  "out/ViseaGO/BP_sResults_clusters_WMA_sole.csv"
)

ViSEAGO::show_table(
  Wang_clusters_wardD2_Immature,
  "out/ViseaGO/BP_sResults_clusters_Immature_sole.csv"
)

ViSEAGO::show_table(
  Wang_clusters_wardD2,
  "out/ViseaGO/BP_sResults_clusters_all_sole.csv"
)

##Reduce and visualize lists of Gene Ontology terms by identifying redundance based on semantic similarity----
library(rrvgo)
library(readr)
GO_analysis_CNA  <- read_delim("out/ViseaGO/BP_sResults_clusters_CNA_sole.csv", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

simMatrix_CNA <- calculateSimMatrix(GO_analysis_CNA$GO.ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")
scores_CNA <- setNames(GO_analysis_CNA$`condition.-log10_pvalue`, GO_analysis_CNA$GO.ID)
reducedTerms_CNA <- reduceSimMatrix(simMatrix_CNA,
                                scores_CNA,
                                threshold=0.7,
                                orgdb="org.Mm.eg.db")
write.csv(reducedTerms_CNA, file = "out/ViseaGo/reducedTerms_CNA.csv")

heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)
scatterPlot(simMatrix, reducedTerms)
treemapPlot(reducedTerms_GLA)
wordcloudPlot(reducedTerms, min.freq=1, colors="black")
rrvgo::shiny_rrvgo()

GO_analysis_GLA  <- read_delim("out/ViseaGO/BP_sResults_clusters_GLA_sole.csv", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)

simMatrix_GLA <- calculateSimMatrix(GO_analysis_GLA$GO.ID,
                                    orgdb="org.Mm.eg.db",
                                    ont="BP",
                                    method="Rel")
scores_GLA <- setNames(GO_analysis_GLA$`condition.-log10_pvalue`, GO_analysis_GLA$GO.ID)
reducedTerms_GLA <- reduceSimMatrix(simMatrix_GLA,
                                    scores_GLA,
                                    threshold=0.7,
                                    orgdb="org.Mm.eg.db")
write.csv(reducedTerms_CNA, file = "out/ViseaGo/reducedTerms_GLA.csv")

GO_analysis_BG  <- read_delim("out/ViseaGO/BP_sResults_clusters_BG_sole.csv", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)

simMatrix_BG <- calculateSimMatrix(GO_analysis_BG$GO.ID,
                                    orgdb="org.Mm.eg.db",
                                    ont="BP",
                                    method="Rel")
scores_BG <- setNames(GO_analysis_BG$`condition.-log10_pvalue`, GO_analysis_BG$GO.ID)
reducedTerms_BG <- reduceSimMatrix(simMatrix_BG,
                                    scores_BG,
                                    threshold=0.7,
                                    orgdb="org.Mm.eg.db")
write.csv(reducedTerms_BG, file = "out/ViseaGo/reducedTerms_BG.csv")

GO_analysis_WMA  <- read_delim("out/ViseaGO/BP_sResults_clusters_WMA_sole.csv", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

simMatrix_WMA <- calculateSimMatrix(GO_analysis_WMA$GO.ID,
                                   orgdb="org.Mm.eg.db",
                                   ont="BP",
                                   method="Rel")
scores_WMA <- setNames(GO_analysis_WMA$`condition.-log10_pvalue`, GO_analysis_WMA$GO.ID)
reducedTerms_WMA <- reduceSimMatrix(simMatrix_WMA,
                                   scores_WMA,
                                   threshold=0.7,
                                   orgdb="org.Mm.eg.db")
write.csv(reducedTerms_WMA, file = "out/ViseaGo/reducedTerms_WMA.csv")

GO_analysis_Immature  <- read_delim("out/ViseaGO/BP_sResults_clusters_Immature_sole.csv", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)

simMatrix_Immature <- calculateSimMatrix(GO_analysis_Immature$GO.ID,
                                    orgdb="org.Mm.eg.db",
                                    ont="BP",
                                    method="Rel")
scores_Immature <- setNames(GO_analysis_Immature$`condition.-log10_pvalue`, GO_analysis_Immature$GO.ID)
reducedTerms_Immature <- reduceSimMatrix(simMatrix_Immature,
                                    scores_Immature,
                                    threshold=0.7,
                                    orgdb="org.Mm.eg.db")
write.csv(reducedTerms_Immature, file = "out/ViseaGo/reducedTerms_Immature.csv")

## Isolate terms shared between clusters
upset <- read_delim("out/ViseaGO/upset.csv", 
                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)

GOs_common <- as.character(upset[upset$combs=="BG-GLA-CNA-WMA-Immature",3])
common <- unlist(strsplit(GOs_common, ";"))
reducedTerms_common <- reducedTerms_BG[which(reducedTerms_BG$go %in% common),]

write.csv(reducedTerms_common, "GO_common.csv", row.names = T)
## Isolate terms specific for single clusters

GOs_BG <- as.character(upset[upset$combs=="BG",3])
BG_exclusive <- unlist(strsplit(GOs_BG, ";"))
reducedTerms_BG_excl <- reducedTerms_BG[which(reducedTerms_BG$go %in% BG_exclusive),]
write.csv(reducedTerms_BG_excl, "GO_BG.csv", row.names = T)


GOs_GLA <- as.character(upset[upset$combs=="GLA",3])
GLA_exclusive <- unlist(strsplit(GOs_GLA, ";"))
reducedTerms_GLA_excl <- reducedTerms_GLA[which(reducedTerms_GLA$go %in% GLA_exclusive),]
write.csv(reducedTerms_GLA_excl, "GO_GLA.csv", row.names = T)


GOs_CNA <- as.character(upset[upset$combs=="CNA",3])
CNA_exclusive <- unlist(strsplit(GOs_CNA, ";"))
reducedTerms_CNA_excl <- reducedTerms_CNA[which(reducedTerms_CNA$go %in% CNA_exclusive),]
write.csv(reducedTerms_CNA_excl, "GO_CNA.csv", row.names = T)


GOs_WMA <- as.character(upset[upset$combs=="WMA",3])
WMA_exclusive <- unlist(strsplit(GOs_WMA, ";"))
reducedTerms_WMA_excl <- reducedTerms_WMA[which(reducedTerms_WMA$go %in% WMA_exclusive),]
write.csv(reducedTerms_WMA_excl, "GO_WMA.csv", row.names = T)

GOs_Immature <- as.character(upset[upset$combs=="Immature",3])
Immature_exclusive <- unlist(strsplit(GOs_Immature, ";"))
reducedTerms_Immature_excl <- reducedTerms_Immature[which(reducedTerms_Immature$go %in% Immature_exclusive),]
write.csv(reducedTerms_Immature_excl, "GO_Immature.csv", row.names = T)

###After manual curation, add info of genes and Pvalue to selected BPs----
library(readxl)
reducedTerms_BG <- read_excel("All_terms_collapsed.xlsx", 
                                  sheet = "Sum Up BG")

BP_sResults_BG <- read.delim("out/ViseaGO/BP_sResults_clusters_BG_sole.csv", 
                                 sep = "\t")

reducedTerms_BG$`-log10Pvalue` <- NA
for (i in 1:nrow(reducedTerms_BG)) {
  for (r in 1:nrow(BP_sResults_BG)) {
    if (is.na(reducedTerms_BG$term[i]) || is.na(BP_sResults_BG$term[r])) {
      next  # Skip to the next iteration if either value is missing
    }
    if (reducedTerms_BG$term[i] == BP_sResults_BG$term[r]) {
      reducedTerms_BG$`-log10Pvalue`[i] <- BP_sResults_BG$condition..log10_pvalue[r]
      reducedTerms_BG$Genes[i] <- BP_sResults_BG$condition.Significant_genes_symbol[r]
    }
  }
}

library(ggplot2)

reducedTerms_BG$group <- c(rep("1", 5), rep("2",3), rep("3",2), rep("4",2), rep("5", 2), rep("6",1))

reducedTerms_BG <- reducedTerms_BG[order(reducedTerms_BG$`-log10Pvalue`), ]

# Create the barplot
ggplot(reducedTerms_BG, aes(x = reorder(term, `-log10Pvalue`), y = `-log10Pvalue`, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +  
  labs(x = "Term", y = "-log10Pvalue") +  
  theme_minimal() +  
  theme(axis.text.y = element_text(size = 8))  

write.csv(reducedTerms_BG, file = "reducedTerms_BG_annotated.csv")

reducedTerms_GLA <- read_excel("All_terms_collapsed.xlsx", 
                              sheet = "Sum Up GLA")

BP_sResults_GLA <- read_delim("out/ViseaGO/BP_sResults_clusters_GLA_sole.csv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

for (i in 1:nrow(reducedTerms_GLA)) {
  for (r in 1:nrow(BP_sResults_GLA)) {
    if (is.na(reducedTerms_GLA$term[i]) || is.na(BP_sResults_GLA$term[r])) {
      next  # Skip to the next iteration if either value is missing
    }
    if (reducedTerms_GLA$term[i] == BP_sResults_GLA$term[r]) {
      reducedTerms_GLA$`-log10Pvalue`[i] <- BP_sResults_GLA$`condition.-log10_pvalue`[r]
      reducedTerms_GLA$Genes[i] <- BP_sResults_GLA$condition.Significant_genes_symbol[r]
    }
  }
}



reducedTerms_GLA$group <- c(rep("1",5),rep("2",3),rep("3",4),rep("4",4))

reducedTerms_GLA <- reducedTerms_GLA[order(reducedTerms_GLA$`-log10Pvalue`), ]


# Create the barplot
ggplot(reducedTerms_GLA, aes(x = reorder(term, `-log10Pvalue`), y = `-log10Pvalue`, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +  
  labs(x = "Term", y = "-log10Pvalue") +  
  theme_minimal() +  
  theme(axis.text.y = element_text(size = 8)) 

write.csv(reducedTerms_GLA, file = "reducedTerms_GLA_annotated.csv")

reducedTerms_CNA <- read_excel("All_terms_collapsed.xlsx", 
                               sheet = "Sum Up CNA")

BP_sResults_CNA <- read_delim("out/ViseaGO/BP_sResults_clusters_CNA_sole.csv", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

for (i in 1:nrow(reducedTerms_CNA)) {
  for (r in 1:nrow(BP_sResults_CNA)) {
    if (is.na(reducedTerms_CNA$term[i]) || is.na(BP_sResults_CNA$term[r])) {
      next  # Skip to the next iteration if either value is missing
    }
    if (reducedTerms_CNA$term[i] == BP_sResults_CNA$term[r]) {
      reducedTerms_CNA$`-log10Pvalue`[i] <- BP_sResults_CNA$`condition.-log10_pvalue`[r]
      reducedTerms_CNA$Genes[i] <- BP_sResults_CNA$condition.Significant_genes_symbol[r]
    }
  }
}

reducedTerms_CNA$group <- c(rep("1",3),rep("2",2), rep("3",3),rep("4",3),rep("5",2)) 
reducedTerms_CNA <- reducedTerms_CNA[order(reducedTerms_CNA$`-log10Pvalue`), ]


# Create the barplot
ggplot(reducedTerms_CNA, aes(x = reorder(term, `-log10Pvalue`), y = `-log10Pvalue`, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +  
  labs(x = "Term", y = "-log10Pvalue") +  
  theme_minimal() +  
  theme(axis.text.y = element_text(size = 8)) 

write.csv(reducedTerms_CNA, file = "reducedTerms_CNA_annotated.csv")

reducedTerms_Immature <- read_excel("All_terms_collapsed.xlsx", 
                               sheet = "Sum Up Immature")

BP_sResults_Immature <- read_delim("out/ViseaGO/BP_sResults_clusters_Immature_sole.csv", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

for (i in 1:nrow(reducedTerms_Immature)) {
  for (r in 1:nrow(BP_sResults_Immature)) {
    if (is.na(reducedTerms_Immature$term[i]) || is.na(BP_sResults_Immature$term[r])) {
      next  # Skip to the next iteration if either value is missing
    }
    if (reducedTerms_Immature$term[i] == BP_sResults_Immature$term[r]) {
      reducedTerms_Immature$`-log10Pvalue`[i] <- BP_sResults_Immature$`condition.-log10_pvalue`[r]
      reducedTerms_Immature$Genes[i] <- BP_sResults_Immature$condition.Significant_genes_symbol[r]
    }
  }
}

reducedTerms_Immature$group <- c(rep("1",11), rep("2",1), rep("3",4), rep("4",3))
reducedTerms_Immature <- reducedTerms_Immature[order(reducedTerms_Immature$`-log10Pvalue`), ]


# Create the barplot
ggplot(reducedTerms_Immature, aes(x = reorder(term, `-log10Pvalue`), y = `-log10Pvalue`, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +  
  labs(x = "Term", y = "-log10Pvalue") +  
  theme_minimal() +  
  theme(axis.text.y = element_text(size = 8)) 

write.csv(reducedTerms_Immature, file = "reducedTerms_Immature_annotated.csv")


# Load the ggplot2 library
library(ggplot2)

# Assuming you have a data.frame named reducedTerms_BG
# with columns 'term' and '-log10Pvalue'

# Create the barplot
ggplot(reducedTerms_BG, aes(x = term, y = `-log10Pvalue`)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +  # To create a longitudinal barplot
  labs(x = "Term", y = "-log10Pvalue") +  # Label the axes
  theme_minimal() +  # Choose a minimal theme
  theme(axis.text.y = element_text(size = 8))  # Adjust the text size on the y-axis