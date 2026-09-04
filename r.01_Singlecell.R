
rm(list = ls())
gc()
options(stringsAsFactors = F) 

.libPaths(c("/home/data/t130525/anaconda3/envs/r-4.2.0/lib/R/library/",.libPaths()))
# a <- .libPaths()
# .libPaths(.libPaths()[-3])

library(scales)
library(ggplot2)
library(cowplot)
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(Seurat)
library(harmony)
# library(kableExtra)
library(ggpubr)
# library(clustree)
# library(scRNAstat) 

# search()

# 1.1 import --------------------------------------------------------------

wd <- "../GSE174574_RAW/"
project_name <- list.files(wd)
raw_data_dir <- paste0(wd, project_name)

creat_scRNA_harmonyect <- function(raw_data_dir,project_name){
  raw_data <- Read10X(raw_data_dir)
  dat_object <- CreateSeuratObject(counts = raw_data, project = project_name, min.cells = 3, min.features = 200)
  cat(paste0(raw_data_dir,': Number of genes: ',dim(dat_object)[1]),file = 'number_log.log',append = T,sep = '\n')
  cat(paste0(raw_data_dir,': Number of cells: ',dim(dat_object)[2]),file = 'number_log.log',append = T,sep = '\n')
  return(dat_object)
}

rawdata <- list()
for (i in 1:length(raw_data_dir)) {
  rawdata[[project_name[i]]] <- creat_scRNA_harmonyect(raw_data_dir[i],project_name[i])
}

scRNA <- merge(x = rawdata[[1]], y = rawdata[-1], project = "GSE174574")
# 18676 features across 58523 samples within 1 assay 

# 1.2 QC_VlnPlot---------------

scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")
QC_pic <- VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster=FALSE, group.by = 'orig.ident',pt.size = 0)
dir.create("01_QC")
pdf('01_QC/01_QC_VlnPlot_before.pdf',height = 8, width = 20)
print(QC_pic)
dev.off()

### QC

scRNA <- subset(scRNA, subset = nCount_RNA <= 15000 & nFeature_RNA <= 4000 & percent.mt <= 20)
# 18676 features across 58188 samples within 1 assay 
QC_pic <- VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster=FALSE, group.by = 'orig.ident',pt.size = 0)
pdf('01_QC/02_QC_VlnPlot_after.pdf',height = 8, width = 20)
print(QC_pic)
dev.off()

# 标准化和特征选择
# scRNA0 <- scRNA
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()

### high_varible_genes

dir.create("02_high_varible_genes")
name = as.character(scRNA@project.name)

hvg_info <- HVFInfo(scRNA)
hvg_info <- hvg_info[order(hvg_info$variance.standardized, decreasing = TRUE), ]
top10 <- rownames(hvg_info)[1:10]

# top10 <- head(VariableFeatures(scRNA), 10)
plot6 <- VariableFeaturePlot(scRNA, pt.size = 2, cols = c("black", "blue"))
plot7 <- LabelPoints(plot = plot6, points = top10, repel = T)+ NoLegend()
pdf(paste0('02_high_varible_genes/',name," visualize high varible genes.pdf"),width=7,height =7)
# p <- CombinePlots(plots = list(plot6, plot7), legend = "bottom")
print(plot7)
dev.off()

# 1.3 Harmony------------------------------------------------------------------

### PCA
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)

# Integrated with harmony
scRNA_harmony <- IntegrateLayers(object = scRNA,
                                 method = HarmonyIntegration, # 需要什么方法改变这里即可
                                 orig.reduction = 'pca', # 这里降维必须选择pca
                                 new.reduction = 'harmony') # 储存在新的降维结果integrated.cca中
scRNA_harmony[['RNA']] <- JoinLayers(scRNA_harmony[['RNA']])

dir.create("03_PCA")
pdf(file = "03_PCA/01_PCA_sample.pdf",width = 7, height = 6)
print(DimPlot(scRNA_harmony, dims = 1:2, reduction = "pca",group.by = 'orig.ident') )
dev.off()

pdf(file = "03_PCA/02_Examine and visualize PCA singlecells with ElbowPlot.pdf", width = 9, height = 9)
print(ElbowPlot(scRNA_harmony, ndims = 50, reduction = "pca"))
dev.off()

# 1.4 降维聚类-----

Idents(scRNA_harmony)
pc.num = 1:10
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = pc.num)
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = pc.num)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = pc.num)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.5)
save(scRNA_harmony, file = "01_scRNA_harmony_no_annotion.rda")

my36colors <-c( '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3','#476D87', '#E95C59', '#E59CC4', 
                '#AB3282', '#23452F', '#BD956A','#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', 
                '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35','#C1E6F3', '#6778AE', '#91D0BE', 
                '#B53E2B', '#712820', '#DCC1DD','#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

# ## 添加分组
sample_case <- c("GSM5319990","GSM5319991","GSM5319992")
scRNA_harmony@meta.data$group <- ifelse(scRNA_harmony@meta.data$orig.ident %in% sample_case, "Case", "Control")

plot_umap_sample <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident",pt.size = 1)+ggtitle('')
plot_umap_group <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "group",pt.size = 1)+ggtitle('')
plot_umap_cluster <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 1,cols = my36colors,
                             label.size = 6)
# plot_sample_cluster <- plot_umap_sample+plot_umap_cluster
pdf(file = "03_PCA/03_umap_cluster.pdf", width = 10, height = 9)
print(plot_umap_sample)
print(plot_umap_group)
print(plot_umap_cluster)
dev.off()

plot_tsne_sample <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "orig.ident",pt.size = 1)+ggtitle('')
plot_tsne_group <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "group",pt.size = 1)+ggtitle('')
plot_tsne_cluster <- DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE, repel = TRUE,pt.size = 1,cols = my36colors,
                             label.size = 6)
# plot_sample_cluster <- plot_tsne_sample+plot_tsne_cluster
pdf(file = "03_PCA/04_tsne_cluster.pdf", width = 10, height = 9)
print(plot_tsne_sample)
print(plot_tsne_group)
print(plot_tsne_cluster)
dev.off()


# # 1.5 FindAllMarkers-------------
# 
# dir.create("06_Findmarker")
# all.markers <- FindAllMarkers(scRNA_harmony, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
# all.markers <- all.markers[order(all.markers$cluster,all.markers$avg_log2FC,decreasing = T),]
# write.csv(all.markers, "06_Findmarker/01_marker_genes.csv")
# # all.markers1 <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = 'roc')
# # all.markers1 <- all.markers1[order(all.markers1$cluster,all.markers1$avg_diff,decreasing = T),]
# # write.csv(all.markers1, "各簇marker genes roc.csv")
# 
# # all.markers %>% group_by(cluster) %>% top_n(n = 9, wt = avg_log2FC) %>%
# #   kable("html",caption = '<center>**表10. 各簇marker genes top9**</center>') %>%
# #   kable_styling(bootstrap_options = "striped", full_width = F)%>%
# #   scroll_box(width = "1000px", height = "500px")
# 
# write.csv(all.markers %>% group_by(cluster) %>% top_n(n = 9, wt = avg_log2FC),file='06_Findmarker/02_marker_genes_top9.csv')
# # write.csv(all.markers1 %>% group_by(cluster) %>% top_n(n = 9, wt = avg_diff),file='各簇marker genes top9 roc.csv')

# 1.6 annotion----------------

dir.create("04_annotion")

features <- c('S100a8',
              'Pf4',
              'Sspo',
              'Kcnj8',
              'Lat',
              'Lum',
              'Ccl17',
              'Aldoc',
              'Hexb',
              'Ttr',
              'Ly6d',
              'Ly6c2',
              'Cldn5',
              'Acta2',
              "Vtn",
              "Rgs5",
              "Tmem212",
              "Dcn",
              "Pdgfra")

##1.6.1 提取表达矩阵和元数据 ----------------------
expr_mat <- GetAssayData(scRNA_harmony, assay = "RNA", layer = "count")
write.csv(expr_mat, "data1.csv")
# 提取元数据（包含细胞类型信息）
meta_data <- scRNA_harmony@meta.data
head(expr_mat)

library(dplyr)
library(tidyr) 
library(tibble)

# 转置为长格式，方便ggplot2绘图
expr_long <- expr_mat[features, ] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("cellname") %>%
  left_join(meta_data %>% rownames_to_column("cellname"), by = "cellname") %>%
  pivot_longer(cols = all_of(features),
               names_to = "gene",
               values_to = "expression")

## 1.6.2. 绘制小提琴图 ----------------------
# 设定细胞类型顺序（和你图中一致）

expr_long$celltype <- factor(expr_long$RNA_snn_res.0.5, levels = 0:19)
expr_long$gene <- factor(expr_long$gene, levels = rev(features))
range(expr_long$expression)

pdf(file = "04_annotion/01_VlnPlot_cluster.pdf", width = 10, height = 9)
ggplot(expr_long, aes(x = celltype, y = expression, fill = celltype)) +
  geom_violin(scale = "width", adjust = 0.5) +
  facet_wrap(~gene, ncol = 1, strip.position = "left") +
  scale_y_log10(limits = c(0.01, 7300), breaks = c(0.1, 1, 10, 100),
                labels = c("0", "1", "2", "3"), name = "Log Expression Level") +
  # scale_fill_manual(values = c(
  #   "SMC" = "#FF66CC", "OLG" = "#CC66FF", "PC" = "#9933FF",
  #   "mix" = "#663399", "MG" = "#3366FF", "LYM" = "#33CCCC",
  #   "MdC" = "#33CC99", "NEUT" = "#33CC66", "EPC" = "#00CC99",
  #   "FB" = "#66CC33", "EC" = "#99CC33", "DC" = "#CCCC33",
  #   "CAM" = "#FFCC33", "ASC" = "#FF6666"
  # )) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    strip.text.y.left = element_text(angle = 0, vjust = 0.5),
    strip.background = element_blank(),
    legend.position = "none"
  )

dev.off()


# DotPlot(scRNA_harmony, features = features, group.by = 'seurat_clusters',assay='RNA') +
#   theme_bw()+
#   theme(panel.grid = element_blank(), 
#         axis.text.x=element_text(hjust = 1,vjust=1,size = 15,angle = 45),
#         axis.text.y=element_text(hjust = 1,vjust=0.5,size = 15),
#   )+
#   labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
#   scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

clustername <- c('EC',
                 'MG',
                 'EC',
                 'MG',
                 'EPC1',
                 'ASC',
                 'SMC',
                 'MdC',
                 'CAM',
                 'OLG',
                 'PC',
                 'MG',
                 'EC',
                 'mix',
                 'MEUT',
                 'EPC2',
                 'DC',
                 'FB',
                 'SMC',
                 'EPC1'
)

Idents(scRNA_harmony) <- scRNA_harmony@meta.data$seurat_clusters
new.cluster.ids <- clustername
names(new.cluster.ids) <- levels(scRNA_harmony)
# Idents(scRNA_harmony) <- scRNA_harmony$seurat_clusters
scRNA_harmony <- RenameIdents(scRNA_harmony, new.cluster.ids)
scRNA_harmony$cell <- Idents(scRNA_harmony)

meta_data <- scRNA_harmony@meta.data
expr_long <- expr_mat[features, ] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("cellname") %>%
  left_join(meta_data %>% rownames_to_column("cellname"), by = "cellname") %>%
  pivot_longer(cols = all_of(features),
               names_to = "gene",
               values_to = "expression")

expr_long$gene <- factor(expr_long$gene, levels = rev(features))

pdf(file = "04_annotion/02_VlnPlot_cell.pdf", width = 10, height = 9)
ggplot(expr_long, aes(x = cell, y = expression, fill = cell)) +
  geom_violin(scale = "width", adjust = 0.5) +
  facet_wrap(~gene, ncol = 1, strip.position = "left") +
  scale_y_log10(limits = c(0.01, 7300), breaks = c(0.1, 1, 10, 100),
                labels = c("0", "1", "2", "3"), name = "Log Expression Level") +
  # scale_fill_manual(values = c(
  #   "SMC" = "#FF66CC", "OLG" = "#CC66FF", "PC" = "#9933FF",
  #   "mix" = "#663399", "MG" = "#3366FF", "LYM" = "#33CCCC",
  #   "MdC" = "#33CC99", "NEUT" = "#33CC66", "EPC" = "#00CC99",
  #   "FB" = "#66CC33", "EC" = "#99CC33", "DC" = "#CCCC33",
  #   "CAM" = "#FFCC33", "ASC" = "#FF6666"
  # )) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    strip.text.y.left = element_text(angle = 0, vjust = 0.5),
    strip.background = element_blank(),
    legend.position = "none"
  )
dev.off()

pdf(file = "04_annotion/03_umap_cluster_cell.pdf", width = 11, height = 9)
DimPlot(scRNA_harmony, reduction = "umap", label = TRUE, repel = F,pt.size = 1,cols = my36colors, label.size = 6)
dev.off()

pdf(file = "04_annotion/04_tsne_cluster_cell.pdf", width = 11, height = 9)
DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE, repel = F,pt.size = 1,cols = my36colors, label.size = 6)
DimPlot(scRNA_PC, reduction = "tsne", label = TRUE, repel = F,pt.size = 1,cols = my36colors, label.size = 6)
dev.off()

save(scRNA_harmony, file='02_scRNA_harmony_annotion.rda')
load("02_scRNA_harmony_annotion.rda")


# 1.7 score ---------------------------------------------------------------

dir.create("05_Score")

ROCK_genes <- readLines("ROCK_Genes.txt") %>% list()
scRNA_harmony <- AddModuleScore(scRNA_harmony,
                                features = ROCK_genes,
                                ctrl = 100,
                                name = "ROCK_genes")
head(scRNA_harmony@meta.data)
#这里就得到了基因集评分结果，但是注意列名为 WNT_features1
colnames(scRNA_harmony@meta.data)[ncol(scRNA_harmony@meta.data)] <- 'ROCK_Score'

pdf(file = "05_Score/01_Score_cell.pdf", width = 10, height = 6)
VlnPlot(scRNA_harmony,features = 'ROCK_Score', 
        pt.size = 0, group.by = "cell")
dev.off()

score_dat <- scRNA_harmony@meta.data[,8:10]
score_dat$group <- factor(score_dat$group, levels = c("Control", "Case"))

pdf(file = "05_Score/02_Score_cell_group.pdf", height = 6,width = 15)
ggplot(score_dat,aes(x = cell, y = ROCK_Score, fill = group)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(values=c("#53A85F","#E95C59"))+
  # ylim(min(score_dat$exp), (max(score_dat$exp)*1.1))+
  # scale_x_discrete()+
  labs(x="",y="", title = "ROCK_Score")+
  theme_bw() + 
  theme(#plot.margin=unit(rep(3,4),'lines'), 
    legend.position = 'top', # 'top','bottom','right','left'
    legend.title=element_blank(),
    text = element_text(size = 20), 
    axis.line = element_line(color = "black"),
    axis.title = element_text(face="bold",size = 20), 
    axis.text.x=element_text(face="bold",size = 18),
    axis.text.y=element_text(size = 12),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank()) + 
  stat_compare_means(label = "p.signif", method = "wilcox.test", 
                     size = 6, hide.ns = T, label.y = (max(score_dat$ROCK_Score)*0.95)) 
dev.off()

# 1.6 PC细胞 ---------------------------------------------------------------

dir.create("06_Findmarker")

# 提取 T 细胞子集
scRNA_PC <- subset(scRNA_harmony, subset = cell == "PC")

# expr_mat1 <- GetAssayData(scRNA_PC, assay = "RNA", layer = "scale.data")
# write.csv(expr_mat1, "data2.csv")
# meta_data1 <- scRNA_PC@meta.data
# write.csv(meta_data1, "data3.csv")


# 设置年龄分组为分组标识
Idents(scRNA_PC) <- scRNA_PC$group

# # 寻找 old 和 young 组间的差异基因
# de_genes <- FindMarkers(scRNA_PC, ident.1 = "Case", ident.2 = "Control", 
#                         min.pct = 0.25, logfc.threshold = 0.5)
# write.csv(de_genes, "06_Findmarker/01_DE_genes.25.csv")
# 
# de_genes <- read.csv("06_Findmarker/01_DE_genes.csv", row.names = 1)
# 
# 
# p_choose <- "p_val"
# adj.P.Val <- 0.05
# logFoldChange <- 0.25
# DEGs_file <- de_genes %>% as.data.frame() 
# DEGs_file$color <- ifelse(DEGs_file[,p_choose] < adj.P.Val & abs(DEGs_file$avg_log2FC) >= logFoldChange,
#                           ifelse(DEGs_file$avg_log2FC > logFoldChange, "Up expression", "Down expression"), 
#                           "Non significant")
# DEGs_file$color <- factor(DEGs_file$color, levels = c("Up expression", "Down expression", "Non significant"))
# colnames(DEGs_file)[colnames(DEGs_file) == p_choose] <- "p_choose"
# table(DEGs_file$color)
# label <- c(rownames(DEGs_file[DEGs_file$color != "Non significant",] %>% slice_max(avg_log2FC, n = 5)), rownames(DEGs_file[DEGs_file$color != "Non significant",] %>% slice_min(avg_log2FC, n = 5)))
# color <- c("Up expression" = "#E63863", "Non significant" = "gray", "Down expression" = "#68A180")
# shape <- c("Up expression" = 24, "Non significant" = 16, "Down expression" = 25)
# p <- ggplot(DEGs_file, aes(avg_log2FC, -log10(p_choose), col = color)) +
#   geom_point(aes(shape = color, fill = color), size=2) +
#   theme_bw() +
#   scale_color_manual(values = color) +
#   scale_shape_manual(values = shape) +
#   labs(x = paste0("log2(Fold-change Case VS Control)"), y = paste0("-log10(", p_choose, ")")) +
#   # geom_label_repel(
#   #   data = DEGs_file[label,],
#   #   aes(label = label),
#   #   size = 5,
#   #   # fill = "darkred", color = "white",
#   #   box.padding = unit(0.35, "lines"),
#   #   point.padding = unit(0.3, "lines"),
#   #   show.legend = F
#   # ) +
#   geom_hline(yintercept = -log10(adj.P.Val), lty = 4, col = "#E39A35", lwd = 1) + 
#   geom_vline(xintercept = c(-logFoldChange, logFoldChange), lty = 4, col = "#E39A35", lwd = 1) +
#   # annotate("text", x = 2, y = 0, label = "|logFC|>1", size = 5, col = "black") +
#   # annotate("text", x = 4, y = 2, label = "adj.p<0.05", size = 5, col = "black") +
#   theme(
#     # plot.margin = unit(rep(3, 4), "lines"), 
#     plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
#     text = element_text(size = 12),
#     axis.title.x = element_text(face = "bold", size = 16),
#     axis.title.y = element_text(face = "bold", size = 16),
#     axis.text.x = element_text(size = 14, angle = 0),
#     axis.text.y = element_text(size = 14, angle = 0),
#     legend.text = element_text(face = "bold", size = 13),
#     legend.title = element_blank(),
#     panel.grid = element_blank(),
#     axis.line = element_line(linewidth = 1, colour = "black"),
#     panel.border = element_blank(),
#     legend.position = "bottom"
#   )
# 
# pdf(file = paste0("06_Findmarker/02_DEGs_volcano.pdf"), height = 8, width = 8, onefile = FALSE)
# print(p)
# dev.off()


scRNA_PC <- FindClusters(scRNA_PC, resolution = 0.3)

clustername <- c("PC1", "PC2", "PC3", "PC4")

Idents(scRNA_PC) <- scRNA_PC@meta.data$seurat_clusters
new.cluster.ids <- clustername
names(new.cluster.ids) <- levels(scRNA_PC)
# Idents(scRNA_harmony) <- scRNA_harmony$seurat_clusters
scRNA_PC <- RenameIdents(scRNA_PC, new.cluster.ids)
scRNA_PC$cell <- Idents(scRNA_PC)

plot_tsne_cluster <- DimPlot(scRNA_PC, reduction = "tsne", label = TRUE, repel = TRUE,pt.size = 1,cols = my36colors,
                             label.size = 6)
# plot_sample_cluster <- plot_tsne_sample+plot_tsne_cluster
pdf(file = "06_Findmarker/03_tsne_cluster.pdf", width = 10, height = 9)
print(plot_tsne_cluster)
dev.off()

scRNA_PC <- AddModuleScore(scRNA_PC,
                                features = ROCK_genes,
                                ctrl = 100,
                                name = "ROCK_genes")
head(scRNA_PC@meta.data)
#这里就得到了基因集评分结果，但是注意列名为 WNT_features1
Idents(scRNA_PC) <- scRNA_PC@meta.data$cell

pdf(file = "06_Findmarker/04_Score_cell.pdf", width = 8, height = 6)
VlnPlot(scRNA_PC,features = 'ROCK_Score', 
        pt.size = 0, group.by = "cell")
dev.off()

score_dat <- scRNA_PC@meta.data[,c(8,9,12)]
colnames(score_dat) <- c("group", "cell", "ROCK_Score")
score_dat$group <- factor(score_dat$group, levels = c("Control", "Case"))

pdf(file = "06_Findmarker/05_Score_cell_group.pdf", height = 6,width = 10)
ggplot(score_dat,aes(x = cell, y = ROCK_Score, fill = group)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(values=c("#53A85F","#E95C59"))+
  # ylim(min(score_dat$exp), (max(score_dat$exp)*1.1))+
  # scale_x_discrete()+
  labs(x="",y="", title = "ROCK_Score")+
  theme_bw() + 
  theme(#plot.margin=unit(rep(3,4),'lines'), 
    legend.position = 'top', # 'top','bottom','right','left'
    legend.title=element_blank(),
    text = element_text(size = 20), 
    axis.line = element_line(color = "black"),
    axis.title = element_text(face="bold",size = 20), 
    axis.text.x=element_text(face="bold",size = 18),
    axis.text.y=element_text(size = 12),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank()) + 
  stat_compare_means(label = "p.signif", method = "wilcox.test", 
                     size = 6, hide.ns = T, label.y = (max(score_dat$ROCK_Score)*0.95)) 
dev.off()

head(scRNA_PC@meta.data)
scRNA_PC <- subset(scRNA_PC, subset = `RNA_snn_res.0.3` == 1)

# 设置年龄分组为分组标识
Idents(scRNA_PC) <- scRNA_PC$group

# 寻找 old 和 young 组间的差异基因
de_genes <- FindMarkers(scRNA_PC, ident.1 = "Case", ident.2 = "Control", 
                        min.pct = 0.25, logfc.threshold = 0.25)
de_genes["Msn",]
write.csv(de_genes, "06_Findmarker/06_DE_genes.25.csv")

p_choose <- "p_val"
adj.P.Val <- 0.1
logFoldChange <- 0.25
DEGs_file <- de_genes %>% as.data.frame() 
DEGs_file$color <- ifelse(DEGs_file[,p_choose] < adj.P.Val & abs(DEGs_file$avg_log2FC) >= logFoldChange,
                          ifelse(DEGs_file$avg_log2FC > logFoldChange, "Up expression", "Down expression"), 
                          "Non significant")
DEGs_file$color <- factor(DEGs_file$color, levels = c("Up expression", "Down expression", "Non significant"))
colnames(DEGs_file)[colnames(DEGs_file) == p_choose] <- "p_choose"
table(DEGs_file$color)
label <- c(rownames(DEGs_file[DEGs_file$color != "Non significant",] %>% slice_max(avg_log2FC, n = 5)), rownames(DEGs_file[DEGs_file$color != "Non significant",] %>% slice_min(avg_log2FC, n = 5)))
color <- c("Up expression" = "#E63863", "Non significant" = "gray", "Down expression" = "#68A180")
shape <- c("Up expression" = 24, "Non significant" = 16, "Down expression" = 25)
p <- ggplot(DEGs_file, aes(avg_log2FC, -log10(p_choose), col = color)) +
  geom_point(aes(shape = color, fill = color), size=2) +
  theme_bw() +
  scale_color_manual(values = color) +
  scale_shape_manual(values = shape) +
  labs(x = paste0("log2(Fold-change Case VS Control)"), y = paste0("-log10(", p_choose, ")")) +
  # geom_label_repel(
  #   data = DEGs_file[label,],
  #   aes(label = label),
  #   size = 5,
  #   # fill = "darkred", color = "white",
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines"),
  #   show.legend = F
  # ) +
  geom_hline(yintercept = -log10(adj.P.Val), lty = 4, col = "#E39A35", lwd = 1) + 
  geom_vline(xintercept = c(-logFoldChange, logFoldChange), lty = 4, col = "#E39A35", lwd = 1) +
  # annotate("text", x = 2, y = 0, label = "|logFC|>1", size = 5, col = "black") +
  # annotate("text", x = 4, y = 2, label = "adj.p<0.05", size = 5, col = "black") +
  theme(
    # plot.margin = unit(rep(3, 4), "lines"), 
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    text = element_text(size = 12),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14, angle = 0),
    axis.text.y = element_text(size = 14, angle = 0),
    legend.text = element_text(face = "bold", size = 13),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.border = element_blank(),
    legend.position = "bottom"
  )

pdf(file = paste0("06_Findmarker/07_DEGs_volcano.pdf"), height = 8, width = 8, onefile = FALSE)
print(p)
dev.off()

all.markers <- FindAllMarkers(scRNA_PC, only.pos = F, min.pct = 0.25, logfc.threshold = 0.5)
all.markers <- all.markers[order(all.markers$cluster,all.markers$avg_log2FC,decreasing = T),]
write.csv(all.markers, "06_Findmarker/08_marker_genes.csv")
# all.markers1 <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = 'roc')
# all.markers1 <- all.markers1[order(all.markers1$cluster,all.markers1$avg_diff,decreasing = T),]
# write.csv(all.markers1, "各簇marker genes roc.csv")

# all.markers %>% group_by(cluster) %>% top_n(n = 9, wt = avg_log2FC) %>%
#   kable("html",caption = '<center>**表10. 各簇marker genes top9**</center>') %>%
#   kable_styling(bootstrap_options = "striped", full_width = F)%>%
#   scroll_box(width = "1000px", height = "500px")

write.csv(all.markers %>% group_by(cluster) %>% top_n(n = 9, wt = avg_log2FC),file='06_Findmarker/09_marker_genes_top9.csv')

features <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
features <- features$gene

expr_mat <- GetAssayData(scRNA_PC, assay = "RNA", layer = "data")
meta_data <- scRNA_PC@meta.data
head(expr_mat)

# 转置为长格式，方便ggplot2绘图
expr_long <- expr_mat[features, ] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("cellname") %>%
  left_join(meta_data %>% rownames_to_column("cellname"), by = "cellname") %>%
  pivot_longer(cols = all_of(features),
               names_to = "gene",
               values_to = "expression")

expr_long$celltype <- factor(expr_long$cell)
expr_long$gene <- factor(expr_long$gene, levels = rev(features))
range(expr_long$expression)

pdf(file = "06_Findmarker/10_VlnPlot_pc_cluster_top5.pdf", width = 10, height = 9)
ggplot(expr_long, aes(x = celltype, y = expression, fill = celltype)) +
  geom_violin(scale = "width", adjust = 0.5) +
  facet_wrap(~gene, ncol = 1, strip.position = "left") +
  # scale_y_log10(limits = c(0.01, 6), breaks = c(0.1, 1, 10, 100),
  #               labels = c("0", "1", "2", "3"), name = "Log Expression Level") +
  # scale_fill_manual(values = c(
  #   "SMC" = "#FF66CC", "OLG" = "#CC66FF", "PC" = "#9933FF",
  #   "mix" = "#663399", "MG" = "#3366FF", "LYM" = "#33CCCC",
  #   "MdC" = "#33CC99", "NEUT" = "#33CC66", "EPC" = "#00CC99",
  #   "FB" = "#66CC33", "EC" = "#99CC33", "DC" = "#CCCC33",
  #   "CAM" = "#FFCC33", "ASC" = "#FF6666"
  # )) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    strip.text.y.left = element_text(angle = 0, vjust = 0.5),
    strip.background = element_blank(),
    legend.position = "none"
  )

dev.off()




# 1.6 细胞异质性 ---------------------------------------------------------------

dir.create("05_Heterogeneity")

a <- table(scRNA_harmony$orig.ident,scRNA_harmony$cell)
write.csv(a, "05_Heterogeneity/01_各类细胞数量.csv")

a <- round(a/rowSums(a),4)
Cellratio <- as.data.frame(a)
write.csv(a, "05_Heterogeneity/02_各类细胞比例.csv")
# Cellratio$group <- ifelse(Cellratio$Var1 %in% c('LZ002','LZ004 ','LZ010','LZ012','LZ017','LZ018','LZ019'),'NOA','NORMAL')

colnames(Cellratio)[1:2] <- c("Sample", "cell")
# Cellratio$Sample <- factor(Cellratio$Sample, levels = c("IP-1","IP-1N", "IP3", "IP-4", "IP-6", "IP-7"))
p <- 
  ggboxplot(Cellratio, x = "cell", y = "Freq",
            color = "Sample", ylab = "cell ratio",
            xlab = "", add = "jitter", size = 1, axis.line = 2, palette=my36colors) +
  # scale_fill_brewer(palette="Set1") +
  # scale_fill_manual(values = my36colors) +
  geom_vline(xintercept = seq(from = 1.5, to = 14.5, by = 1), lty = 4, col = "#E39A35", lwd = 1) +
  # stat_compare_means(label = "p.signif", label.y = max(Cellratio[, 3]), aes(group = Var1),method = 'kruskal.test') + 
  rotate_x_text(45)

pdf(file = "05_Heterogeneity/03_各类群细胞含量比率统计.pdf", width = 12, height = 6)
print(p)
dev.off()
# Cellratio$group <- factor(Cellratio$group,levels = unique(Cellratio$group),ordered = T)
# Cellratio <- Cellratio[order(Cellratio$group),]
# Cellratio$Var1 <- factor(Cellratio$Var1,levels = unique(Cellratio$Var1),ordered = T)
p1 <- ggplot(data = Cellratio,aes(x=Sample,y=Freq,fill=cell))+
  geom_bar(stat="identity",position = "fill")+
  ##添加辅助线
  geom_hline(yintercept=0.25,linetype=2,size=1)+
  geom_hline(yintercept=0.50,linetype=2,size=1)+
  geom_hline(yintercept=0.75,linetype=2,size=1)+
  ##旋转坐标轴
  # coord_flip()+
  xlab("")+
  # scale_fill_brewer(value = my36colors)+
  scale_fill_manual(values = my36colors)+
  theme_classic()+
  theme(
    ##设置图例上方，字体，大小，颜色
    legend.position = "top",legend.title = element_blank(),
    legend.text=element_text(size = 14,colour = "black"),
    ####设置边框线粗细，颜色，类型
    line = element_line(colour = "black", size = 1, linetype = 1), 
    ###设置坐标轴标题字体，颜色，大小
    axis.title=element_text(size = 16,face="bold",colour = "black"),
    ###设置坐标轴标签字体，颜色，大小
    axis.text.y = element_text(size = 14,colour = "black"),
    axis.text.x = element_text(size = 12,colour = "black"))

pdf(file = "05_Heterogeneity/04_各类群细胞含量堆积图.pdf", width = 12, height = 6)
print(p1)
dev.off()


box <- Cellratio
box$Sample <- ifelse(substr(box$Sample, 1, 1) == "Y", "Young", "Old")
box$Sample <- factor(box$Sample, levels = c("Young", "Old"))
box$cell <- factor(box$cell, levels = unique(box$cell))

pdf(file = "05_Heterogeneity/05_box.pdf", height = 6,width = 10)
ggplot(box,aes(x = cell,y = Freq, fill = Sample)) + 
  geom_boxplot(alpha=0.7) + 
  scale_fill_manual(values=c("#53A85F","#E95C59"))+ 
  ylim(min(box$Freq), (max(box$Freq)*1.1))+
  #scale_x_discrete()+
  labs(x="",y="Percent")+
  theme_bw() + 
  theme(#plot.margin=unit(rep(3,4),'lines'), 
    legend.position = 'top', # 'top','bottom','right','left'
    legend.title=element_blank(),
    text = element_text(size = 20), 
    axis.line = element_line(color = "black"),
    axis.title = element_text(face="bold",size = 20), 
    axis.text.x=element_text(face="bold",size = 18,angle=45,vjust = 1,hjust = 1),
    axis.text.y=element_text(size = 12),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank()) + 
  stat_compare_means(method = "wilcox", 
                     label = "p.signif",
                     # label = "p.format",
                     size = 8, 
                     hide.ns = T, 
                     label.y = (max(box$Freq)*1.05)) 
dev.off()




# 1.7 hub gene----------------

scRNA_harmony@meta.data$cell_group <- paste0(scRNA_harmony@meta.data$cell, " - ", scRNA_harmony@meta.data$group)
features <- c("BNC2", "MRAP2", "PLEKHA2")

dir.create("06_hub_gene")

pdf('06_hub_gene/01_hub-DotPlot-cell.pdf',height = 5, width = 8)
DotPlot(scRNA_harmony, features = features, group.by = "cell", assay='RNA') + RotatedAxis()
dev.off()

pdf('06_hub_gene/02_hub-DotPlot-group.pdf',height = 3, width = 6)
DotPlot(scRNA_harmony, features = features, group.by = "group", assay='RNA') + RotatedAxis()
dev.off()

pdf('06_hub_gene/03_hub-DotPlot-cell_dis.pdf',height = 8, width = 8)
DotPlot(scRNA_harmony, features = features, group.by = "cell_group", assay='RNA') + RotatedAxis()
dev.off()

# pdf('06_hub_gene/04_hub-VlnPlot-cell.pdf',height = 12, width = 18)
# VlnPlot(scRNA, features = features, group.by = "cell", assay='RNA', ncol = 3, raster=FALSE)
# dev.off()
# 
# pdf('06_hub_gene/hub-VlnPlot-cell_dis.pdf',height = 12, width = 21)
# VlnPlot(scRNA, features = features, group.by = "cell_group", assay='RNA', ncol = 3, raster=FALSE)
# dev.off()

DefaultAssay(scRNA) <- "RNA"

pdf('06_hub_gene/04_hub-FeaturePlot-cell.pdf',height = 6, width = 19)
FeaturePlot(scRNA_harmony, features = features, ncol = 3, pt.size = 0.1, raster=FALSE) 
# scale_colour_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
dev.off()


# 1.8 cell-phone--------

library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
options(stringsAsFactors = FALSE)


dir.create("07_cellphone")

# 创建cellchat对象

scRNA_harmony[["RNA4"]] <- as(object = scRNA_harmony[["RNA"]], Class = "Assay")
DefaultAssay(scRNA_harmony) <- "RNA"

table(scRNA_harmony@meta.data$group)

# 按分组拆分对象
cellchat_group <- SplitObject(scRNA_harmony, split.by = 'group')

# 循环处理每个分组
i <- "Control"
for (i in names(cellchat_group)) {
  # 获取当前分组的对象
  object <- cellchat_group[[i]]
  
  # 创建 CellChat 对象
  cellchat <- createCellChat(object@assays$RNA4@data)
  
  # 添加元数据
  meta <- data.frame(cellType = object$cell, row.names = Cells(object))
  cellchat <- addMeta(cellchat, meta = meta, meta.name = "cell")
  cellchat <- setIdent(cellchat, ident.use = "cell") # set "labels" as default cell identity
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  
  # 加载与设定需要的 CellChatDB 数据库
  CellChatDB <- CellChatDB.mouse ## CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted signaling for cell-cell communication analysis
  cellchat@DB <- CellChatDB.use # set the used database in the object
  
  # 预处理表达数据以进行细胞间相互作用分析
  cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  cellchat <- identifyOverExpressedGenes(cellchat) 
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- smoothData(cellchat, adj = PPI.mouse)
  
  # 根据表达值推测细胞互作的概率（cellphonedb 是用平均表达值代表互作强度）。
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) #如果不想用上一步 PPI 矫正的结果，raw.use = TRUE 即可。
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  write.csv(df.net, paste0("07_cellphone/", i, "_net_lr.csv"))
  
  cellchat <- computeCommunProbPathway(cellchat)
  df.netp <- subsetCommunication(cellchat, slot.name = "netP")
  write.csv(df.netp, paste0("07_cellphone/", i, "_net_pathway.csv"))
  
  # 统计细胞和细胞之间通信的数量（有多少个配体 - 受体对）和强度（概率）
  cellchat <- aggregateNet(cellchat)
  # 计算每种细胞各有多少个
  groupSize <- as.numeric(table(cellchat@idents))
  
  # 定义一个非常小的非零值
  epsilon <- 1e-6
  
  # 将值为 0 的数据替换为 epsilon
  cellchat@net$count[cellchat@net$count == 0] <- epsilon
  cellchat@net$weight[cellchat@net$weight == 0] <- epsilon
  
  # 绘图
  # par(mfrow = c(1, 2), xpd = TRUE)
  p1 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                         label.edge = T, title.name = "Number of interactions")
  p2 <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                         label.edge = F, title.name = "Interaction weights/strength")
  p3 <- plot_grid(p1, p2, ncol = 2)
  pdf(paste0("07_cellphone/", i, "_cellphone.pdf"), width = 14, height = 7)
  print(p3)
  # ggsave(paste0("07_cellphone/", i, "_cellphone.pdf"), width = 14, height = 7)
  dev.off()
  
  print(i)
}





















# 1.7 hub box---------------------------------------------------------------------

mymatrix <- as.data.frame(scRNA_harmony@assays$RNA4@data)
mymatrix2<-t(mymatrix)%>%as.data.frame()
mymatrix2[,1]<-scRNA_harmony$cell
mymatrix3 <- data.frame(cell = scRNA_harmony$cell,
                        group = scRNA_harmony$group,
                        mymatrix2[,features])

pic <- list()

for (i in 3:4) {
  box <- mymatrix3[,c(1:2,i)]
  colnames(box)[3] <- "exp"
  pic[[colnames(mymatrix3)[i]]] <- 
    ggplot(box,aes( x = cell, y = exp, fill = group)) + 
    geom_boxplot(alpha=0.7) + 
    scale_fill_manual(values=c("#53A85F","#E95C59"))+ 
    #scale_x_discrete()+
    labs(x="",y="Expression level", title = colnames(mymatrix3)[i])+
    theme_bw() + 
    theme(#plot.margin=unit(rep(3,4),'lines'), 
      legend.position = 'top',
      legend.title=element_blank(),
      text = element_text(size = 20), 
      axis.line = element_line(color = "black"),
      axis.title = element_text(face="bold",size = 22), 
      axis.text.x=element_text(face="bold",size = 18,vjust = 1, hjust = 1, angle = 45),
      axis.text.y=element_text(size = 12),
      panel.border = element_blank(),
      panel.background = element_blank(),
      panel.grid.major =element_blank(),
      panel.grid.minor = element_blank()) + 
    stat_compare_means(label = "p.signif", method = "t.test", size = 10,
                       hide.ns = T, label.y = (max(box[,3])*0.95)) 
}

pdf(paste0('06_hub_gene/05.hub.box.pdf'),height = 7,width = 10)
print(pic)
dev.off()

# 1.8 time ---------------------------------------------------------------------

library(monocle)

rm(list = ls())
gc()

load("02_scRNA_harmony_annotion.rda")
dir.create(("08_time"))

## 样本拆分

# object_time <- scRNA_harmony[,which(Idents(scRNA_harmony) %in% c( "NK cells"))]
object_time <- scRNA_PC
object_time[["RNA4"]] <- as(object = object_time[["RNA"]], Class = "Assay")
DefaultAssay(object_time) <- "RNA4"

# object_time@meta.data$cellname <- object_time@assays$RNA4@data@Dimnames[[2]]
# sample <- object_time@meta.data$cellname 
# sample <- unique(sample)

# set.seed(438)
# index <- caTools::sample.split(sample, SplitRatio = 0.3)
# object_time <- object_time[,index]

table(object_time$cell)
# rm(sample_all)
# gc()

##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
exprMat <- as(as.matrix(object_time@assays$RNA4@counts), 'sparseMatrix')
##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
p_data <- object_time@meta.data
#p_data$celltype <- object_time@active.ident ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
##提取基因信息 如生物类型、gc含量等
f_data <- data.frame(gene_short_name = row.names(object_time@assays$RNA4),row.names = row.names(object_time@assays$RNA4))
##exprMat的行数与f_data的行数相同(gene number), exprMat的列数与p_data的行数相同(cell number)

pd <- new('AnnotatedDataFrame', data = p_data)

fd <- new('AnnotatedDataFrame', data = f_data)

cds <- newCellDataSet(exprMat, phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1) #过滤基因，将会在fData(cds)中添加一列num_cells_expressed
# expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10)) #过滤掉在小于10个细胞中表达的基因

#使用2000高边基因
# DefaultAssay(object_time) <- "integrated"
express_genes <- VariableFeatures(object_time)
cds<- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)

#降维
cds <- reduceDimension(cds,max_components = 2,method = 'DDRTree')
saveRDS(cds,"object_time_monocle2.rds")
cds <- readRDS("object_time_monocle2.rds")

### 重启，仅加载igraph2.0.3和monocle

#拟时序轴轨迹构建并排列

# install.packages("https://mirrors.westlake.edu.cn/CRAN/src/contrib/Archive/igraph/igraph_2.0.3.tar.gz",repos = NULL,type = "source")
cds <- orderCells(cds) ### long times
saveRDS(cds,"object_time_monocle2.rds")

# cds <- readRDS("object_time_monocle2.rds")
#按照细胞状态
pdf(file = "08_time/01.cell_object_time_time.pdf", width = 5, height = 4.5)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE)+ scale_colour_gradient(low = "#4DBBD5FF", high = "#DC0000FF")
dev.off()

pdf(file = "08_time/02.cell_object_time_State.pdf", width = 5, height = 4.5)
plot_cell_trajectory(cds, color_by = "State", cell_size =0.7,show_backbone=TRUE)
dev.off()

pdf(file = "08_time/03.cell_object_time_cell.pdf", width = 5, height = 4.5)
plot_cell_trajectory(cds, color_by = "cell", cell_size =0.7,show_backbone=TRUE)
dev.off()

pdf(file = "08_time/04.cell_object_time_Group.pdf", width = 5, height = 4.5)
plot_cell_trajectory(cds, color_by = "group", cell_size =0.7)+
  scale_colour_manual(values = c("#91D1C2FF","#FF7F00")
                      #"Control" = "#91D1C2FF","IgAN" = "#FF7F00"
  )   #,cols =c("Control"="#91D1C2FF","IgAN"="#FF7F00")
dev.off()

hub <- c("Rhoa","Rock1","Rock2")
hub <- c("Myl9","Acta2","Msn")
cds_subset <- cds[hub,]
##可视化：以state/celltype/pseudotime进行
p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "cell")
p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
plotc <- p1|p2|p3
ggsave("08_time/05.Genes_pseudotimeplot1.pdf", plot = p3, width = 8, height = 8)
ggsave("08_time/05.Genes_Stateplot1.pdf", plot = p1, width = 8, height = 8)
ggsave("08_time/05.Genes_cellplot1.pdf", plot = p2, width = 8, height = 8)




# 1.11 SCENIC -------------------------------------------------------

dir.create("10_SCENIC")
setwd("10_SCENIC")

setwd("/home/tom/AnTao/Project/SingleCell_scMetabolism/SCENIC/cisTarget_databases")
dir.create("SCENIC")
setwd("SCENIC")

library(SCENIC)
library(SCopeLoomR)
library(AUCell)
library(arrow)

defaultDbNames

# https://resources.aertslab.org/cistarget/
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")
# mc9nr: Motif collection version 9: 24k motifs
dbFiles
for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
  #  (1041.7 MB)
  # 
}
setwd("..")

load("../02_scRNA_harmony_annotion.rda")

scRNA_PC[["RNA4"]] <- as(object = scRNA_PC[["RNA"]], Class = "Assay")
DefaultAssay(scRNA_PC) <- "RNA4"

scRNA_PC@meta.data$cellname <- scRNA_PC@assays$RNA4@data@Dimnames[[2]]
sample <- scRNA_PC@meta.data$cellname
sample <- unique(sample)

set.seed(438)
index <- caTools::sample.split(sample, SplitRatio = 0.1)
scRNA_PC <- scRNA_PC[,index]


exprMat <- as.matrix(scRNA_PC@assays$RNA4@data)
dim(exprMat)
exprMat[1:4,1:4] 
cellInfo <-  scRNA_PC@meta.data[,c(9,3,2)]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
head(cellInfo)
table(cellInfo$CellType)

library(parallel)
detectCores()
# data(package = "RcisTarget")

options(RcisTarget_rankingsFeatureCol = "motifs")
data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9

# Sys.setenv(LIBARROW_MINIMAL = "false")
# Sys.setenv(ARROW_WITH_ZSTD = "ON") 
# install.packages("arrow", repos = c(arrow = "https://apache.r-universe.dev"))
# devtools::install_version("arrow", version = "14.0.0.2")
install.packages("../SCENIC_1.1.2.tar.gz", repos=NULL)

## 加载4.2,library(arrow)

scenicOptions <- initializeScenic(org="mgi", 
                                  dbDir="cisTarget_databases", nCores=16) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
scenicOptions <- readRDS("int/scenicOptions.Rds")

### Co-expression network
# genesKept <- geneFiltering(exprMat, scenicOptions)
genesKept <- geneFiltering(exprMat, scenicOptions = scenicOptions,
                           minCountsPerGene = 3*.01*ncol(exprMat),
                           minSamples = ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings


regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"])
rssPlot <- plotRSS(rss)

pdf(file = "plotRSS.pdf", width = 6, height = 9)
print(rssPlot$plot)
dev.off()



##导入原始regulonAUC矩阵
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC

##导入二进制regulonAUC矩阵
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN




library(pheatmap)
cellInfo <- readRDS("int/cellInfo.Rds")
subtype = cellInfo[order(cellInfo$CellType),]
# subtype = subset(cellInfo,select = 'subtype')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
#挑选部分感兴趣的regulons
# my.regulons <- c('CD59_extended_34g','CEBPD_extended_14g','CHD2_extended_668g','TEAD1_extended_68g','TEAD1_61g',
#                  'JUNB_extended_24g','CREB3L4_37g','SPDEF_56g','CREB3L4_extended_124g','IRF1_31g','IRF1_extended_40g')
# myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
# myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]
myAUCmatrix <- AUCmatrix[,rownames(subtype)]
myBINmatrix <- BINmatrix[,rownames(subtype)]

#使用regulon原始AUC值绘制热图

pdf(file = "AUCmatrix.pdf", width = 8, height = 9)
pheatmap(myAUCmatrix, show_colnames=F, annotation_col=subtype, cluster_cols = F)
dev.off()
#使用regulon二进制AUC值绘制热图
pdf(file = "BINmatrix.pdf", width = 8, height = 9)
pheatmap(myBINmatrix, show_colnames=F, annotation_col=subtype, cluster_cols = F,
         color = colorRampPalette(colors = c("white","black"))(100))
dev.off()

