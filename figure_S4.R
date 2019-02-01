
#read in the brca1
if(!exists("allBRCA1")){
  allBRCA1 <- fread("~/mahmood/binom/analyzed data/analyzed_allBRCA1diffmodels.csv")
}

#read in the TP53
if(!exists("allTP53")){
  allTP53 <- fread("~/mahmood/binom/analyzed data/analyzed_allTP53diffmodels.csv")
}
allTP53[CancerModel == "Lung Cancer", CancerModel := "Lung Adenocarcinoma"]
#read in the TP53
if(!exists("allMart")){
  allMart <- fread("~/mahmood/binom/analyzed data/analyzed_allMartelottodiffmodels.csv")
}

allBRCA1$label <- factor(allBRCA1$label, levels = c("Neutral","Deleterious"))


all_B1 <- ggplot(allBRCA1, aes(label, log10(mutability), fill = label)) + 
  geom_boxplot() + facet_grid(. ~ CancerModel) + scale_fill_manual(values = c("#ffeda0","#f03b20")) + theme(axis.title = element_blank()) +
  theme(strip.text= element_blank(),axis.title = element_blank(), axis.text.x = element_blank()) +  scale_y_continuous(breaks=c(-7, -6,-5,-4),labels = c("7","6","5","4")) +
  theme(legend.position = "none" ) +
  theme(plot.margin=unit(c(0,0,0,1.2),"cm"))


allMart[,label := as.character(label)]
allMart$label <- factor(allMart$label, levels = c("0","1"))


allMart_plot <- ggplot(allMart, aes(label, log10(mutability), fill = label)) + 
  geom_boxplot() + facet_grid(. ~ CancerModel) + scale_fill_manual(values = c("#ffeda0","#f03b20")) + theme(axis.title = element_blank()) +
  theme(strip.text= element_blank(),axis.title = element_blank(), axis.text.x = element_blank()) +  scale_y_continuous(breaks=c(-7, -6,-5,-4),labels = c("7","6","5","4")) +
  theme(legend.position = "none" ) +
  theme(plot.margin=unit(c(0,0,0,1.2),"cm"))






allTP53$label <- factor(allTP53$label, levels = c("functional","partially-functional", "non-functional"))


all_T <- ggplot(allTP53, aes(label, log10(mutability), fill = label)) + 
  geom_boxplot() + facet_grid(. ~CancerModel) + scale_fill_brewer(palette = "YlOrRd",labels = c("Functional","Partially Functional", "Damaging")) + theme(axis.title = element_blank()) +
  theme(axis.title = element_blank(), axis.text.x = element_blank()) + 
  scale_y_continuous(breaks=c(-7, -6,-5, -4),labels = c("7","6","5","4")) +
  theme(legend.position = "none") + 
  scale_x_discrete(labels = c("Functional","Partially-Functional", "Damaging")) +
  theme(strip.text= element_text(size = 9))
  theme(plot.margin=unit(c(0,0,0,1.2),"cm"))
  

fig6 <- ggarrange(all_T, all_B1,allMart_plot, nrow = 3, labels = c("A","B","C"),align = "hv")

fig6 <- annotate_figure(fig6,left = text_grob(expression(bold(paste("-Log"[10], "Mutability"))),rot = 90))
fig6


breast<- readPNG("~/Desktop/9032.png")
lung <- readPNG("~/Desktop/9026.png")
pancan <- readPNG("~/Desktop/33999.png")
skin <- readPNG("~/Desktop/9025.png")
snp <- readPNG("~/Desktop/34275.png")

rH = 0.095
rW = 0.15 
disp = 0.18
first = -.31
picy = .435

ggdraw() + draw_plot(fig6,0, 0, 1, .90) + 
  draw_plot(rasterGrob(breast,height = rH,width = rW),x = first, y = picy) +
  draw_plot(rasterGrob(lung,height = rH,width = rW),x = first + disp, y = picy) +
  draw_plot(rasterGrob(skin,height = rH,width = rW),x = first + 2*disp, y = picy) +
  draw_plot(rasterGrob(pancan,height = rH,width = rW),x = first + 3*disp, y = picy) +
  draw_plot(rasterGrob(snp,height = rH,width = rW),x = first + 4*disp, y = picy) 


brca_p <- data.table()  
for(m in unique(allBRCA1[,CancerModel])){
  print(m)
  test <- (wilcox.test(allBRCA1[label == "Neutral" & CancerModel == m, mutability],allBRCA1[label != "Neutral" & CancerModel == m, mutability]))
  brca_p <- rbind(brca_p, cbind(m,test$p.value,test$statistic))
  
}

names(brca_p) <- c("Model", "Pvalue", "W Statistic")

p53_p <- data.table()  

for(m in unique(allTP53[,CancerModel])){
  print(m)
  test <- dunn.test(allTP53[CancerModel == m, mutability],allTP53[CancerModel == m, label],method = "by")
  print(test)
  p53_p <- rbind(p53_p, cbind(as.data.table(test),rep(m,3)),fill = T)
}


martelotto_p <- data.table()  
for(m in unique(allMart[,CancerModel])){
  print(m)
  test <- (wilcox.test(allMart[label == "0" & CancerModel == m, mutability],allMart[label == "1" & CancerModel == m, mutability]))
  martelotto_p <- rbind(martelotto_p, cbind(m,test$p.value,test$statistic))
  
}

names(martelotto_p) <- c("Model", "Pvalue", "W Statistic")


fwrite(p53_p, "~/mahmood/binom/analyzed data/TP53 statistical comparisons Diff Background.csv")
fwrite(brca_p, "~/mahmood/binom/analyzed data/BRCA1 statistical comparisons Diff Background.csv")
fwrite(martelotto_p, "~/mahmood/binom/analyzed data/Martelotto statistical comparisons Diff Background.csv")
