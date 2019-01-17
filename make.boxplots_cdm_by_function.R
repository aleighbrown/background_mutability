####make sure that the algo for p53 has run, and that the algo for the brcas has run
#source('~/mahmood/functional data/process_brca2.R')

#read in BRCA1
if(!exists("keepBRCA")){
  keepBRCA <- fread("~/mahmood/binom/analyzed data/keepBRCA2.csv")
}

if(!exists("keepMo")){
  keepMo <- fread("~/mahmood/binom/analyzed data/keepMartelotto2.csv")
}
#read in p53
if(!exists("keepP53")){
  keepP53 <- fread("~/mahmood/binom/analyzed data/keepP532.csv")
}
#read in the combined data set
if(!exists("full_combined_dataset")){
  full_full_combined_dataset <- fread("~/mahmood/binom/analyzed data/final_benchmark.csv")
}


alp = 0.5
liney = -5
# #set the order of the Transactivation class so that it's going to plot in the order I like
# keepP53$TransactivationClass <- factor(keepP53$TransactivationClass, levels = c("functional","partially-functional", "non-functional"))
# #what am I going to compare
# 
# p <- ggplot(keepP53, aes(TransactivationClass, log10(mutability), fill = TransactivationClass)) + 
#   geom_boxplot() + 
#   geom_jitter(aes(TransactivationClass, log10(mutability)),position=position_jitter(width=0.1,height=0),size = 2, alpha = alp) + 
#   scale_fill_brewer(palette = "YlOrRd") +
#   ylab("Mutability") +
#   theme(axis.title = element_blank()) +
#   theme(legend.position = "none") +
#   scale_x_discrete(labels = c("Functional","Partially-Functional", "Non-functional")) +
#   scale_y_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) 
# 
# 
# 
# ####adding lines for significance
# p <- p +  geom_segment(aes(x = "functional", y = -4.95, xend = "non-functional", yend = -4.95),size = 1) +
#   annotate("text", x = "partially-functional", y = -4.95, label = "**", size = 12) +
#   geom_segment(aes(x = 2.1, y = -5.1, xend = 2.9, yend = -5.1),size = 1) +
#   annotate("text", x = 2.5, y = -5.1, label = "**", size = 12)
# 
# ####brca1 is from the deep mutational scanning experiment
# 
# keepBRCA$Annotation <- factor(keepBRCA$Annotation, levels = c("Benign","Deleterious"))
# 
# 
# b1 <- ggplot(keepBRCA, aes(Annotation, log10(mutability), fill = Annotation)) + 
#   geom_boxplot() + 
#   geom_jitter(aes(Annotation, log10(mutability)),position=position_jitter(width=0.1,height=0),size = 2, alpha = alp) + 
#   scale_fill_manual(values = c("#ffeda0","#f03b20")) +
#   theme(axis.title=element_blank()) +
#   theme(axis.title = element_blank()) +
#   theme(legend.position = "none") +
#   scale_y_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) 
# 

# b1 <- b1 + geom_segment(aes(x = 1.1, y = liney, xend = 1.9, yend = liney),size = 1) +
#   annotate("text", x = 1.5, y = -4.95, label = "*", size = 12)
# 
# ####martelotto
# 
# #set order for plotting
# keepMo$Functionalclassification <- factor(keepMo$Functionalclassification , levels = c("neutral","non-neutral"))
# 
# mar <- ggplot(keepMo, aes(Functionalclassification, log10(mutability), fill = Functionalclassification)) + 
#   geom_boxplot() + 
#   geom_jitter(aes(Functionalclassification, log10(mutability)),position=position_jitter(width=0.1,height=0),size = 2, alpha = alp) + 
#   scale_fill_manual(values = c("#ffeda0","#f03b20")) +
#   theme(axis.title=element_blank()) +
#   theme(axis.title = element_blank()) +
#   theme(legend.position = "none") +
#   scale_y_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) +
#   scale_x_discrete(labels = c("Neutral","Non-neutral"))
# 
# mar <- mar + geom_segment(aes(x = 1.1, y = liney, xend = 1.9, yend = liney),size = 1) +
#   annotate("text", x = 1.5, y = -4.95, label = "**", size = 12)

#set order for plotting
full_combined_dataset[label == 0,Annotation2 := "Neutral"]
full_combined_dataset[label == 1,Annotation2 := "Non-neutral"]
full_combined_dataset$Annotation2 <- factor(full_combined_dataset$Annotation2 , levels = c("Neutral","Non-neutral"))

combined <- ggplot(full_combined_dataset, aes(Annotation2, log10(Mutability), fill = Annotation2)) + 
  geom_boxplot() + 
  geom_jitter(aes(Annotation2, log10(Mutability)),position=position_jitter(width=0.1,height=0),size = 2, alpha = alp) + 
  scale_fill_manual(values = c("#ffeda0","#f03b20")) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) +
  ylab(expression(bold(paste("-Log"[10], "Mutability")))) +
  scale_x_discrete(labels = c("Neutral","Non-neutral"))

combined <- combined + geom_segment(aes(x = 1.1, y = liney, xend = 1.9, yend = liney),size = 1) +
  annotate("text", x = 1.5, y = -4.955, label = "**", size = 12)
#add some extra space so things align properly 

group_annotation = ggplot(full_combined_dataset, aes (x = CosmicGroup, y = log10(Mutability), fill = CosmicGroup)) + 
  facet_grid(rows = "Annotation2") + geom_boxplot() + 
  scale_fill_manual(values = rep("#3C77AF",4)) + 
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) +
  ylab(expression(bold(paste("-Log"[10], "Mutability"))))







fig5 <- ggarrange(combined,group_annotation, labels = c("A", "B"), nrow = 2)




group_annotationMSK = ggplot(full_combined_dataset, aes (x = MSKGroup, y = log10(Mutability), fill = MSKGroup)) + 
  facet_grid(rows = "Annotation2") + geom_boxplot() + 
  scale_fill_manual(values = rep("#B3DBB8",4)) + 
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) +
  ylab(expression(bold(paste("-Log"[10], "Mutability"))))





####dunn test to compare between the three groups for TP53
p53_test = dunn.test(keepP53$mutability,keepP53$TransactivationClass)
### mann u test to compare betwen BRCA1, and Martelotto
#wilcox.test
wilcox.test(mutability ~ Functionalclassification, keepMo)
wilcox.test(mutability ~ Annotation, keepBRCA)
#wilcox test to compare between the combinded dataset
#wilcox.test(mutability ~ label, full_combined_dataset)
