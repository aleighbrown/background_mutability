#read in a csv containing all the mutations for all the census gens
#read in the census table
if(!exists("amino_rebuilt")){
  amino_rebuilt <- fread("~/mahmood/binom/analyzed data/rebuilt_amino_table.csv")
}
#all

####now by oncogene tsg

oncoTSG <- ggplot(subset(amino_rebuilt, BinnedRole %in% c("Oncogene","TSG")),
                  aes(x=BinnedRole, y=log10(AminoMutabilityPan),fill = SubType, alpha = AminoCosmicGroupPan)) +
  geom_boxplot() +
  facet_grid(. ~ SubType) +
  scale_fill_manual(values = c("#00BFC4","#F8766D","#7CAE00")) +
  theme(axis.title = element_blank(), legend.position = "none") +
  scale_y_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) 

xvals <- c(1.28125, 1.09375, 0.90625, 0.71875, 2.28125, 2.09375, 1.90625, 1.71875)
xvals <- sort(xvals)

group_labels <- c("0","1","2","3+","0","1","2","3+")

oncoTSG <- annotate_groups(oncoTSG,xvals, yval = -7, fsize = 5, labels = group_labels)

#version 1
#fig3 <- ggarrange(all, domRec, oncoTSG, ncol = 1, nrow = 3, labels = c("A", "B", "C"))


fig3 <- annotate_figure(oncoTSG,left = text_grob(expression(bold(paste("-Log"[10], " Mutability"))), rot = 90))

fig3
##i guess we can do tests as well 
# #correlation between mutability and frequency
# allamino_census$group <- as.factor(allamino_census$group)
# allamino_census$RoleinCancer <- as.factor(allamino_census$RoleinCancer)
# 
# kruskal.test(allamino_census$mutability, allamino_census$SubType)
# dunn.test::dunn.test(allamino_census$mutability, allamino_census$SubType)

