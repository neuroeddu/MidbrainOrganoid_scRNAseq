#!/usr/bin/env Rscript

# import all the libraries
library("ggplot2")
library("Seurat")
library("cowplot")
library("clustree")
library(patchwork)
library(dplyr)
library("Matrix")

# testing a script for using in base R so I can run remotely on grumio or compute canada
print("I'm starting the script")
# read in data
# file and pathway
filepath <- "/home/rhalena/scRNAseq/BrainCommRevision/CombinedIntegratedAST23CON_21-07-2021.rds"
MBO.combined <- readRDS(filepath)

# set the file pathway
# setwd("~/scRNAseq/BrainComm/outs")

# get cluster markers for different resolutions 
# DefaultAssay(MBO.combined) <- "integrated"
# i = 0.1
# 
# MBO.combined <- FindClusters(MBO.combined, resolution = i)
# DefaultAssay(MBO.combined) <- "RNA"
# # I want to get the cluster markers from the RNA not just the values used for integration
# ClusterMarkers <- FindAllMarkers(MBO.combined, logfc.threshold = 0.1, return.thresh = 0.05, only.pos = TRUE)
# top20 <-ClusterMarkers %>% group_by(cluster) %>% top_n(n=20, wt = avg_log2FC)
# top5 <- ClusterMarkers %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC)
# # save all the positive markers
# write.csv(ClusterMarkers, "SClusterMarkersPos21072021-res0.1.csv")
# # save the top 20 positive markers
# write.csv(top20, "SClusterMarkersTOP20_Pos21072021-res0.1.csv")
# 
# MBO.combined <- ScaleData(MBO.combined, verbose = FALSE)
# png("HeatMaptop5Markersres0.1.png")
# DoHeatmap(MBO.combined, features = top5$gene, size=3, angle =90, group.bar.height = 0.02)
# dev.off()
# 
# # make heat map of selected markers
# feature_list = c("DLX2","SOX2","PAX6","SOX9","HES1","ITGA6","MAP2","NCAM1","CD24","GRIA2","GRIN2B","GABBR1","TH","CORIN","CALB1","KCNJ6","CXCR4","ITGA6","NES","SLC1A3","CD44","AQP4","S100B","GFAP","ALDH1L1","EAAT1","EAAT2", "PDGFRA","CLDN11","VIM","VCAM1")
# ft2 = c("ACE2","ACE1","SOX9","HES1","ITGA6","MAP2","NCAM1","CD24","GRIA2","GRIN2B","GABBR1","TH","CORIN","CALB1","KCNJ6","CXCR4","ITGA6","NES","SLC1A3","CD44","AQP4","S100B", "PDGFRA","CLDN11","VIM","VCAM1")
# 
# png("HeatMapMarkerlistres0.1.png")
# DoHeatmap(MBO.combined, group.by = "seurat_clusters", features = feature_list)
# dev.off()



# DefaultAssay(MBO.combined) <- "integrated"
# i = 0.3
# 
# MBO.combined <- FindClusters(MBO.combined, resolution = i)
# DefaultAssay(MBO.combined) <- "RNA"
# # I want to get the cluster markers from the RNA not just the values used for integration
# ClusterMarkers <- FindAllMarkers(MBO.combined, logfc.threshold = 0.1, return.thresh = 0.05, only.pos = TRUE)
# top20 <-ClusterMarkers %>% group_by(cluster) %>% top_n(n=20, wt = avg_log2FC)
# top5 <- ClusterMarkers %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC)
# # save all the positive markers
# write.csv(ClusterMarkers, "SClusterMarkersPos21072021-res0.3.csv")
# # save the top 20 positive markers
# write.csv(top20, "SClusterMarkersTOP20_Pos21072021-res0.3.csv")
# 
# MBO.combined <- ScaleData(MBO.combined, verbose = FALSE)
# png("HeatMaptop5Markersres0.3.png")
#   DoHeatmap(MBO.combined, features = top5$gene, size=3, angle =90, group.bar.height = 0.02)
# dev.off()
# 
# # make heat map of selected markers
# feature_list = c("DLX2","SOX2","PAX6","SOX9","HES1","ITGA6","MAP2","NCAM1","CD24","GRIA2","GRIN2B","GABBR1","TH","CORIN","CALB1","KCNJ6","CXCR4","ITGA6","NES","SLC1A3","CD44","AQP4","S100B","GFAP","ALDH1L1","EAAT1","EAAT2", "PDGFRA","CLDN11","VIM","VCAM1")
# ft2 = c("ACE2","ACE1","SOX9","HES1","ITGA6","MAP2","NCAM1","CD24","GRIA2","GRIN2B","GABBR1","TH","CORIN","CALB1","KCNJ6","CXCR4","ITGA6","NES","SLC1A3","CD44","AQP4","S100B", "PDGFRA","CLDN11","VIM","VCAM1")
# 
# png("HeatMapMarkerlistres0.3.png")
#   DoHeatmap(MBO.combined, group.by = "seurat_clusters", features = feature_list)
# dev.off()

# resolution 0.6

DefaultAssay(MBO.combined) <- "integrated"
i = 0.6

MBO.combined <- FindClusters(MBO.combined, resolution = i)
DefaultAssay(MBO.combined) <- "RNA"
# I want to get the cluster markers from the RNA not just the values used for integration
ClusterMarkers <- FindAllMarkers(MBO.combined, logfc.threshold = 0.1, return.thresh = 0.05, only.pos = TRUE)
top20 <-ClusterMarkers %>% group_by(cluster) %>% top_n(n=20, wt = avg_log2FC)
top5 <- ClusterMarkers %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC)
# save all the positive markers
write.csv(ClusterMarkers, "SClusterMarkersPos21072021-res0.6.csv")
# save the top 20 positive markers
write.csv(top20, "SClusterMarkersTOP20_Pos21072021-res0.6.csv")

MBO.combined <- ScaleData(MBO.combined, verbose = FALSE)
png("HeatMaptop5Markersres0.6.png")
DoHeatmap(MBO.combined, features = top5$gene, size=3, angle =90, group.bar.height = 0.02)
dev.off()

# make heat map of selected markers
feature_list = c("DLX2","SOX2","PAX6","SOX9","HES1","ITGA6","MAP2","NCAM1","CD24","GRIA2","GRIN2B","GABBR1","TH","CORIN","CALB1","KCNJ6","CXCR4","ITGA6","NES","SLC1A3","CD44","AQP4","S100B","GFAP","ALDH1L1","EAAT1","EAAT2", "PDGFRA","CLDN11","VIM","VCAM1")
ft2 = c("ACE2","ACE1","SOX9","HES1","ITGA6","MAP2","NCAM1","CD24","GRIA2","GRIN2B","GABBR1","TH","CORIN","CALB1","KCNJ6","CXCR4","ITGA6","NES","SLC1A3","CD44","AQP4","S100B", "PDGFRA","CLDN11","VIM","VCAM1")

png("HeatMapMarkerlistres0.6.png")
DoHeatmap(MBO.combined, group.by = "seurat_clusters", features = feature_list)
dev.off()


# resolution 0.9

DefaultAssay(MBO.combined) <- "integrated"
i = 0.9

MBO.combined <- FindClusters(MBO.combined, resolution = i)
DefaultAssay(MBO.combined) <- "RNA"
# I want to get the cluster markers from the RNA not just the values used for integration
ClusterMarkers <- FindAllMarkers(MBO.combined, logfc.threshold = 0.1, return.thresh = 0.05, only.pos = TRUE)
top20 <-ClusterMarkers %>% group_by(cluster) %>% top_n(n=20, wt = avg_log2FC)
top5 <- ClusterMarkers %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC)
# save all the positive markers
write.csv(ClusterMarkers, "SClusterMarkersPos21072021-res0.9.csv")
# save the top 20 positive markers
write.csv(top20, "SClusterMarkersTOP20_Pos21072021-res0.9.csv")

MBO.combined <- ScaleData(MBO.combined, verbose = FALSE)
png("HeatMaptop5Markersres0.9.png")
DoHeatmap(MBO.combined, features = top5$gene, size=3, angle =90, group.bar.height = 0.02)
dev.off()

# make heat map of selected markers
feature_list = c("DLX2","SOX2","PAX6","SOX9","HES1","ITGA6","MAP2","NCAM1","CD24","GRIA2","GRIN2B","GABBR1","TH","CORIN","CALB1","KCNJ6","CXCR4","ITGA6","NES","SLC1A3","CD44","AQP4","S100B","GFAP","ALDH1L1","EAAT1","EAAT2", "PDGFRA","CLDN11","VIM","VCAM1")
ft2 = c("ACE2","ACE1","SOX9","HES1","ITGA6","MAP2","NCAM1","CD24","GRIA2","GRIN2B","GABBR1","TH","CORIN","CALB1","KCNJ6","CXCR4","ITGA6","NES","SLC1A3","CD44","AQP4","S100B", "PDGFRA","CLDN11","VIM","VCAM1")

png("HeatMapMarkerlistres0.9.png")
DoHeatmap(MBO.combined, group.by = "seurat_clusters", features = feature_list)
dev.off()


# resolution 1.2

DefaultAssay(MBO.combined) <- "integrated"
i = 1.2

MBO.combined <- FindClusters(MBO.combined, resolution = i)
DefaultAssay(MBO.combined) <- "RNA"
# I want to get the cluster markers from the RNA not just the values used for integration
ClusterMarkers <- FindAllMarkers(MBO.combined, logfc.threshold = 0.1, return.thresh = 0.05, only.pos = TRUE)
top20 <-ClusterMarkers %>% group_by(cluster) %>% top_n(n=20, wt = avg_log2FC)
top5 <- ClusterMarkers %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC)
# save all the positive markers
write.csv(ClusterMarkers, "SClusterMarkersPos21072021-res1.2.csv")
# save the top 20 positive markers
write.csv(top20, "SClusterMarkersTOP20_Pos21072021-res1.2.csv")

MBO.combined <- ScaleData(MBO.combined, verbose = FALSE)
png("HeatMaptop5Markersres1.2.png")
DoHeatmap(MBO.combined, features = top5$gene, size=3, angle =90, group.bar.height = 0.02)
dev.off()

# make heat map of selected markers
feature_list = c("DLX2","SOX2","PAX6","SOX9","HES1","ITGA6","MAP2","NCAM1","CD24","GRIA2","GRIN2B","GABBR1","TH","CORIN","CALB1","KCNJ6","CXCR4","ITGA6","NES","SLC1A3","CD44","AQP4","S100B","GFAP","ALDH1L1","EAAT1","EAAT2", "PDGFRA","CLDN11","VIM","VCAM1")
ft2 = c("ACE2","ACE1","SOX9","HES1","ITGA6","MAP2","NCAM1","CD24","GRIA2","GRIN2B","GABBR1","TH","CORIN","CALB1","KCNJ6","CXCR4","ITGA6","NES","SLC1A3","CD44","AQP4","S100B", "PDGFRA","CLDN11","VIM","VCAM1")

png("HeatMapMarkerlistres1.2.png")
DoHeatmap(MBO.combined, group.by = "seurat_clusters", features = feature_list)
dev.off()

# resolution 1.8
DefaultAssay(MBO.combined) <- "integrated"
i = 1.8

MBO.combined <- FindClusters(MBO.combined, resolution = i)
DefaultAssay(MBO.combined) <- "RNA"
# I want to get the cluster markers from the RNA not just the values used for integration
ClusterMarkers <- FindAllMarkers(MBO.combined, logfc.threshold = 0.1, return.thresh = 0.05, only.pos = TRUE)
top20 <-ClusterMarkers %>% group_by(cluster) %>% top_n(n=20, wt = avg_log2FC)
top5 <- ClusterMarkers %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC)
# save all the positive markers
write.csv(ClusterMarkers, "SClusterMarkersPos21072021-res1.8.csv")
# save the top 20 positive markers
write.csv(top20, "SClusterMarkersTOP20_Pos21072021-res1.8.csv")

MBO.combined <- ScaleData(MBO.combined, verbose = FALSE)
png("HeatMaptop5Markersres1.8.png")
DoHeatmap(MBO.combined, features = top5$gene, size=3, angle =90, group.bar.height = 0.02)
dev.off()

# make heat map of selected markers
feature_list = c("DLX2","SOX2","PAX6","SOX9","HES1","ITGA6","MAP2","NCAM1","CD24","GRIA2","GRIN2B","GABBR1","TH","CORIN","CALB1","KCNJ6","CXCR4","ITGA6","NES","SLC1A3","CD44","AQP4","S100B","GFAP","ALDH1L1","EAAT1","EAAT2", "PDGFRA","CLDN11","VIM","VCAM1")
ft2 = c("ACE2","ACE1","SOX9","HES1","ITGA6","MAP2","NCAM1","CD24","GRIA2","GRIN2B","GABBR1","TH","CORIN","CALB1","KCNJ6","CXCR4","ITGA6","NES","SLC1A3","CD44","AQP4","S100B", "PDGFRA","CLDN11","VIM","VCAM1")

png("HeatMapMarkerlistres1.8.png")
DoHeatmap(MBO.combined, group.by = "seurat_clusters", features = feature_list)
dev.off()


























  
print("Script complete - check for your data : ) ")
