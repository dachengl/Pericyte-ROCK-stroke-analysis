## Plot
pacman::p_load(ggplot2, ggpubr, aplot, ggvenn, psych, ggcorrplot, corrplot, pROC, UpSetR, patchwork, VennDiagram, eulerr, ggcor, ggunchained, ggrepel) 
## Machine learning
pacman::p_load(glmnet, e1071, caret, doParallel, parallel, kableExtra, Boruta, neuralnet, NeuralNetTools, caret) 

file_name <- "02_Expression"; dir.create(file_name); setwd(file_name)

## 5.4 key gene--------------------------------

DEG_PC <- readLines("01_DEG_PC2.txt")
ROCK <- readLines("../../01_rawdata/02_gene/ROCK_Genes.txt")

hub <- Reduce(intersect, list(DEG_PC, ROCK))
# hub <- Reduce(intersect, list(sel, fea_gene))
write.table(hub, file = "02_hub_gene.txt", row.names = F, col.names = F, quote=F)

upset_list <- list(
  DEGs_PC = DEG_PC,
  ROCK_Genes = ROCK
)
pdf(file = "02_hub_venn.pdf",height=4,width=5,onefile=FALSE)
ggvenn(upset_list, 
       stroke_linetype = 2, 
       stroke_size = 1,
       text_size = 4,
       text_color = "black",
       # set_name_color = "red", 
       set_name_size = 4,
       fill_color = c("#91D0BE", "gold","pink"))
dev.off()

## 5.5 ROC ---------------------------------------------------------------------------------------------------------

GSE_accession <- "GSE23160"
exp_xl <- read.csv("../00_cleandata/expression_GSE23160.csv", row.names = 1, check.names = F)
group_xl <- read.csv("../00_cleandata/group_GSE23160.csv")
rownames(group_xl) <- group_xl$sample
group_xl$group <- factor(group_xl$group, levels = c("Control", "Case"))
exp_roc <- exp_xl
group_roc <- group_xl

GSE_accession <- "GSE104036"
exp_yz <- read.csv("../00_cleandata/expression_GSE104036.csv", row.names = 1, check.names = F)
group_yz <- read.csv("../00_cleandata/group_GSE104036.csv")
rownames(group_yz) <- group_yz$sample
group_yz$group <- factor(group_yz$group, levels = c("Control", "Case"))
# exp_roc <- exp_yz
exp_roc <- log2(exp_yz + 1)
group_roc <- group_yz

hub <- c("Myl9",
        "Acta2",
        "Msn"
)

# hub <- hub[c(2,3,5,6)]

roc_data <- data.frame(group = group_roc$group, t(exp_roc[hub,] %>% na.omit()), check.names = F)
# roc_data <- data.frame(group = group_roc$group, t(log2(exp_roc[hub,] + 1) %>% na.omit()), check.names = F)
write.csv(roc_data, file = paste0("03_roc_data_", GSE_accession, ".csv"), row.names = F)

## expression

box <- roc_data %>% gather(gene, exp, -group)
box$gene <- factor(box$gene, levels = hub)

Controlnumber <- table(roc_data$group)[1]
d <- c()
e <- c()
for(i in colnames(roc_data)[-1]){
  a = median(roc_data[c(1:Controlnumber),i])
  b = median(roc_data[-c(1:Controlnumber),i])
  pa <-  wilcox.test(roc_data[c(1:Controlnumber),i],roc_data[-c(1:Controlnumber),i])
  p = pa$p.value
  c = if_else(p>0.05,"black",if_else(a>b,"#68A180","#8C549C"))
  d <- c(d,c)
  e <- c(e,c,c)
}

pdf(file = paste0("03_hub_box_", GSE_accession, ".pdf"), height = 5,width = 8)
ggplot(box,aes(x = gene,y = exp, fill = group)) + 
  geom_boxplot(alpha=0.7) + 
  scale_fill_manual(values=c("#53A85F","#E95C59"))+ 
  ylim(min(box$exp), (max(box$exp)*1.1))+
  #scale_x_discrete()+
  labs(x="",y="Expression level", title = GSE_accession)+
  theme_bw() + 
  theme(#plot.margin=unit(rep(3,4),'lines'), 
    legend.position = 'top', # 'top','bottom','right','left'
    legend.title=element_blank(),
    text = element_text(size = 20), 
    axis.line = element_line(color = "black"),
    axis.title = element_text(face="bold",size = 20), 
    axis.text.x=element_text(face="bold",size = 18,angle=0,vjust = 0.5),
    axis.text.y=element_text(size = 12),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank()) + 
  stat_compare_means(label = "p.signif", method = "wilcox", 
                     size = 8, hide.ns = T, label.y = (max(box$exp)*1)) 
dev.off()



setwd("..")
