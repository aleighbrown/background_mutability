
full_combined_dataset <- fread("~/mahmood/binom/analyzed data/final_benchmark.csv")

#setnames(full_combined_dataset, c("fathmScore","CountNew","binomial","CSB","mutability", "vest","candra_score"), c("FatHMM","Frequency","B-Score","CHASM","Mutability", "VEST", "CanDrA"))

scores =  c("FatHMM","Frequency","B-Score","CHASM","Mutability", "VEST", "CanDrA","REVEL","LR")

negsvec =  get_negs_vector(full_combined_dataset, scores)

combined_roc_table <- make_roc_table(full_combined_dataset,scores = scores, negs = negsvec)
 
combined_roc_table <- combined_roc_table[!(Measure == "Frequency" & alpha == 0)]
combined_roc_table <- combined_roc_table[!(Measure == "LR" & alpha == 0)]

combined_pr_table <- make_pr_table(full_combined_dataset,scores = scores ,negs = negsvec)

combined_pr_table <- combined_pr_table[!(Measure == "Frequency" & alpha == 0)]
combined_pr_table <- combined_pr_table[!(Measure == "LR" & alpha == 0)]

one <- plot_roc_table(combined_roc_table)
one_withinset <- make_inset_ROC(one,xmin = 0.3,ymax = 0.65) + theme(legend.position = c(0.7,0.7))
 
one_withinset <- annotate_figure(one_withinset,bottom = text_grob("False positive rate", face = "bold"), left = text_grob("True positive rate", rot = 90,face = "bold"))

two <- plot_pr_table(combined_pr_table) + theme(legend.position = "none")
#two_withinset <- make_inset_ROC(two,zoomXmin = 0,zoomXmax = 0.5, zoomYmin = 0.5, zoomYmax = 1,ab = F,xbreaks = c(0,0.25,0.5), xmin = 0.2,xmax = 0.76,ymin = -0.05,ybreaks = c(0.5, 0.66,.85,1))
two_withinset <- annotate_figure(two,bottom = text_grob("Recall", face = "bold"), left = text_grob("Precision", rot = 90,face = "bold"))
both <- ggarrange(one_withinset, two_withinset, labels = c("A","B"))
both
png(file="mygraphic.png",width=800,height=400,res=72)
both
dev.off()




scores2 =  c("FatHMM","B-Score","CHASM","Mutability", "VEST", "CanDrA","REVEL")

negsvec2 =  get_negs_vector(full_combined_dataset, scores2)
