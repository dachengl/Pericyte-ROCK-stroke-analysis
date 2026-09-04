## Immune--ssGSEA
pacman::p_load(GSEABase, GSVA, RColorBrewer, tibble, cowplot, tidyr, msigdbr, xCell, immunedeconv, tidyverse) 

# 8. Immune infiltration-ssGESA ----------------------------------------------------------------------------------------------------------

file_name <- "04_Immune"
dir.create(file_name)
setwd(file_name)

# data

exp_xl <- exp_roc
exp_xl <- cbind(symbol = rownames(exp_xl), exp_xl)
exp_xl$symbol <- toupper(exp_xl$symbol)
exp_xl <- aggregate(. ~ symbol, data = exp_xl, max)
rownames(exp_xl) <- exp_xl$symbol
exp_xl <- exp_xl[,-1]
group_xl <- group_roc

range(exp_xl)

data <- exp_xl
Group <- group_xl
gmtFile="D:/00_Laborer/00_1741/01_Database/00_code/Immunity/mmc3.gmt"
geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

## 7.1 Start ssGSEA analysis --------

ssgseaScore=gsva(as.matrix(data), geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=T) #"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts

range(ssgseaScore)
ssgseaOut <- data.frame(t(ssgseaScore), check.names = F)
ssgseaOut <- data.frame(group = Group$group, ssgseaOut, check.names = F)

## 7.2 Organize data format-----

result1 <- ssgseaOut
write.csv(result1,"00_ssGSEA_result.csv")
# result1 <- read.csv("00_ssGSEA_result.csv", row.names = 1)

re1 <- result1
colnames(re1)[1] <- "Type"
table(re1$Type)  

re2 = re1[,-1]    
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
dat_cell <- re2 %>% as.data.frame() %>%rownames_to_column("Sample") %>%gather(key = Cell_type,value = Proportion,-Sample)
dat_group = gather(re1,Cell_type,Proportion,-Type )
dat = cbind(dat_cell,dat_group$Type)
colnames(dat)[4] <- "Type"

## 7.3 box###################

box=dat
box$Cell_type <- factor(box$Cell_type, levels = colnames(re1)[-1])
Controlnumber <- table(re1$Type)[1]
d <- c()
e <- c()
for(i in colnames(re1)[-1]){
  a = mean(re1[c(1:Controlnumber),i])
  b = mean(re1[-c(1:Controlnumber),i])
  pa <-  t.test(re1[c(1:Controlnumber),i],re1[-c(1:Controlnumber),i])
  p = pa$p.value
  c = if_else(p>0.05,"black",if_else(a>b,"#53A85F","#E95C59"))
  d <- c(d,c)
  e <- c(e,c,c)
}

pdf(file="01_ssGSEA_box.pdf",width=15,height=8,onefile=FALSE)
ggplot(box,aes( x = Cell_type,y = Proportion, fill = Type)) + 
  geom_boxplot(alpha=0.7) + 
  scale_fill_manual(values=c("#53A85F","#E95C59"))+ 
  #scale_x_discrete()+
  labs(x = "Cell Type", y = "Immune cell Score")+
  theme_bw() + 
  theme(
    #plot.margin=unit(rep(3,4),'lines'), 
    legend.position = 'top',
    text = element_text(size = 20), 
    axis.line = element_line(color = "black"),
    axis.title = element_text(face="bold",size = 22), 
    axis.text.x=element_text(size = 16, vjust = 1, hjust = 1, angle = 45,colour = d),
    axis.text.y=element_text(size = 12),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank()) + 
  stat_compare_means(label = "p.signif", method = "t.test", hide.ns = T, size = 6)
dev.off()

## 7.4 免疫细胞相关性----------------

cor <- corr.test(re2, method = "spearman")

m = par(no.readonly = T)
pdf(file="02_cor_imm.pdf",width = 10,height = 10)
# corrplot(corr = cor$r,
#          p.mat = cor$p,
#          diag = F,
#          method = "number",
#          type = "lower",
#          cl.pos = "r",
#          tl.pos = "l",
#          insig="blank",
#          col = colorRampPalette(c("#476D87", "pink", "#E95C59"))(100)
#          # col = "black"
# )
corrplot(corr = cor$r,
         p.mat = cor$p,
         add = F,
         diag = T,
         method = "square",
         type = "upper",
         sig.level = c(0.001, 0.01, 0.05),
         insig = 'label_sig',
         # cl.pos = "n",
         tl.pos = "td",
         pch.cex = 1,
         col = colorRampPalette(c("#476D87", "pink", "#E95C59"))(100))
dev.off()

## 7.4 cells and genes-------------------------------------------------------------------------

rownames(re2) == group_xl$sample

cor_hub_imm <- corr.test(re2[,d != "black"], t(exp_roc[hub,]), method = "spearman", adjust="BH")
write.csv(cor_hub_imm$r, file = "02_hub_imm_cor.csv")
write.csv(cor_hub_imm$p, file = "02_hub_imm_pvalue.csv")


max(cor_hub_imm$r) # 0.8283
min(cor_hub_imm$r) # -0.483

exp_imm <- re2[,d != "black"]
exp_imm <- re2
exp_hub <- t(exp_roc[hub,])

pdf('03_cor_hub_imm.pdf',height = 6,width = 20)
quickcor(exp_hub,exp_imm,method="spearman",cor.test = T)+ 
  #geom_circle2()+ 
  geom_star(n = 5)+
  theme(text = element_text(size = 13),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        #legend.position=c(0.9,0.9),
        legend.position="right",
        legend.text=element_text(size = 13))+
  scale_fill_gradient2(low = "#6D9EC1", high = "#FF7F00")+ 
  # geom_mark(r = NA,size = 5) + #框内显著性
  geom_mark(size = 5) + #框内显著性
  #geom_cross()+ 
  scale_size_manual(values = c(2, 0.5, 2)) +
  scale_colour_manual(values = c("#FF9933", "#6D9EC1", "white"))
dev.off()

setwd("..")


