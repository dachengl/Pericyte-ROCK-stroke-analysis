
## GSEA
pacman::p_load(msigdbr, enrichplot, clusterProfiler, GSVA) 
pacman::p_load(clusterProfiler, org.Hs.eg.db, org.Mm.eg.db, treemap, GOplot, tidyverse, ggraph, ccgraph, tidygraph, ggsci, scales) 

# 7. GSEA --------------------------------------------------------------------

file_name <- "03_GSEA"
dir.create(file_name)
setwd(file_name)

df.expr <- exp_roc
all.gene <- AnnotationDbi::select(org.Mm.eg.db, rownames(df.expr), "ENTREZID", "SYMBOL")
exprSet <- merge(df.expr, all.gene, by.x = "row.names", by.y  = "SYMBOL") %>% 
  dplyr::distinct(ENTREZID, .keep_all = T) %>% na.omit()
rownames(exprSet) <- NULL
exprSet <- exprSet %>% tibble::column_to_rownames(var = "ENTREZID") %>% select(-Row.names)

hub_gene <- hub
hub_gene = AnnotationDbi::select(org.Mm.eg.db, hub_gene, "ENTREZID", "SYMBOL")
entrez.gene <- hub_gene$ENTREZID
df.gene <- exprSet[entrez.gene, ] %>% as.data.frame()
df.gene2 <- cbind(GeneID = hub_gene$SYMBOL, df.gene)
write.csv(df.gene2, file = "01_hub_gene_expr.csv", row.names = F)

gsea.plot = function(res.kegg, top.hall, gene){
  gsdata <- do.call(rbind, lapply(top.hall, enrichplot:::gsInfo, object = res.kegg))
  gsdata$Description = factor(gsdata$Description, levels = top.hall)
  p1 = ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(14) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_brewer(palette = "Set1") +
    ggtitle(gene) +
    geom_hline(yintercept = 0, color = "black", size = 0.8) +
    geom_line(aes_(y = ~runningScore, color = ~Description), size = 1) +
    theme(legend.position = "right", legend.title = element_blank(), legend.background = element_rect(fill = "transparent")) +
    ylab("Running Enrichment Score") + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank(),
          text = element_text(face = "bold", family = "Times"))
  i = 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1 }
  p2 = ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
    theme_classic(14) + theme(legend.position = "none", 
                              axis.ticks = element_blank(), 
                              axis.text = element_blank(), 
                              axis.line.x = element_blank(),
                              text = element_text(face = "bold", family = "Times")) + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    scale_color_brewer(palette = "Set1")
  p = aplot::insert_bottom(p1, p2, height = 0.15)
  return(p)
}

df.m = msigdbr(species = "Mus musculus")
df.kegg = subset(df.m, gs_subcat == "CP:KEGG")[c(3,5)]
df.exp <- exprSet %>% t %>% as.data.frame()

if (!dir.exists("02_GSEA")) {dir.create("02_GSEA")}
if (!dir.exists("03_GSEA_Res")) {dir.create("03_GSEA_Res")}
lapply(1:nrow(hub_gene), function(i){
  hub = rownames(df.gene2)[i] %>% as.character()
  hub.exp = df.exp[[hub]]
  hub.cor = cor(df.exp, hub.exp, method = "spearman") %>% as.data.frame %>% na.omit
  hub.coreff = hub.cor[[1]]
  names(hub.coreff) = rownames(hub.cor)
  hub.coreff = hub.coreff[order(hub.coreff, decreasing = T)]
  res.kegg = GSEA(hub.coreff, TERM2GENE = df.kegg, pvalueCutoff = 0.2, seed = 1, pAdjustMethod = "none", eps = 0)
  write.table(res.kegg, file = paste0("03_GSEA_Res/0", i, ".", hub_gene$SYMBOL[which(hub_gene$ENTREZID==hub)], ".res.xls"),
              sep = "\t", row.names = T, col.names = T, quote = F)
  top.kegg = res.kegg@result
  top.kegg = top.kegg[order(top.kegg$p.adjust, decreasing = F),]$Description[1:5]
  p = gsea.plot(res.kegg, top.kegg, hub_gene$SYMBOL[which(hub_gene$ENTREZID==hub)])
  fn1 = paste0("02_GSEA/", sprintf("%02d",i), ".", hub_gene$SYMBOL[which(hub_gene$ENTREZID==hub)], ".png")
  fn2 = paste0("02_GSEA/", sprintf("%02d",i), ".", hub_gene$SYMBOL[which(hub_gene$ENTREZID==hub)], ".pdf")
  ggsave(fn1, p, width = 12, height = 6, units = "in", limitsize = 300)
  ggsave(fn2, p, width = 12, height = 6, units = "in", limitsize = 300)
  return(0)
})

setwd("..")
