
full_combined_dataset <- fread("~/mahmood/binom/analyzed data/final_benchmark.csv")
full_combined_dataset[,MSKBScore := binom.test(x = (MSKFreq + 1), p = Mutability, n = 9228, alt = "g")$p.value, by = 1:nrow(full_combined_dataset)]

#setnames(full_combined_dataset, c("fathmScore","CountNew","binomial","CSB","mutability", "vest","candra_score"), c("FatHMM","Frequency","B-Score","CHASM","Mutability", "VEST", "CanDrA"))


scores =  c("FatHMM","Frequency","B-Score","CHASMplus", "VEST", "CanDrA","REVEL","CHASM","Mutability")

negsvec =  get_negs_vector(full_combined_dataset, scores)

combined_roc_table <- make_roc_table(full_combined_dataset,scores = scores, negs = negsvec)
 
combined_roc_table <- combined_roc_table[!(Measure == "Frequency" & alpha == 0)]

combined_pr_table <- make_pr_table(full_combined_dataset,scores = scores ,negs = negsvec)

combined_pr_table <- combined_pr_table[!(Measure == "Frequency" & alpha == 0)]

one <- plot_roc_table(combined_roc_table)
one_withinset <- make_inset_ROC(one,xmin = 0.3,ymax = 0.55,ymin = -0.07) + theme(legend.position = c(0.7,0.7))
 
one_withinset <- annotate_figure(one_withinset,bottom = text_grob("False positive rate", face = "bold", size = 18), left = text_grob("True positive rate", rot = 90,face = "bold", size = 18))

two <- plot_pr_table(combined_pr_table) + theme(legend.position = "none")
#two_withinset <- make_inset_ROC(two,zoomXmin = 0,zoomXmax = 0.5, zoomYmin = 0.5, zoomYmax = 1,ab = F,xbreaks = c(0,0.25,0.5), xmin = 0.2,xmax = 0.76,ymin = -0.05,ybreaks = c(0.5, 0.66,.85,1))
two_withinset <- annotate_figure(two,bottom = text_grob("Recall", face = "bold", size = 18), left = text_grob("Precision", rot = 90,face = "bold", size = 18))
both <- ggarrange(one_withinset, two_withinset, labels = c("A","B"))
both

make_score_table(full_combined_dataset,scores = scores, negs = negsvec, roundto = 2)

scores_msk =  c("FatHMM","MSKFreq","MSKBScore","CHASMplus","Mutability", "VEST", "CanDrA","REVEL","CHASM")
negsvecmsk =  get_negs_vector(full_combined_dataset, scores_msk)

make_score_table(full_combined_dataset[Frequency != 0],scores = scores, negs = negsvec, roundto = 2)

make_score_table(full_combined_dataset[MSKFreq != 0],scores = scores_msk, negs = negsvecmsk, roundto = 4)



###Rare
scores_rare =  c("FatHMM","B-Score","CHASMplus", "VEST", "CanDrA","REVEL","CHASM")
scores_raremsk =  c("FatHMM","MSKBScore","CHASMplus", "VEST", "CanDrA","REVEL","CHASM")

negsvec_rare =  get_negs_vector(full_combined_dataset, scores_rare)
negsvec_raremsk =  get_negs_vector(full_combined_dataset, scores_raremsk)

make_score_table(full_combined_dataset[Frequency == 0],scores = scores_rare, negs = negsvec_rare)
make_score_table(full_combined_dataset[Frequency == 1],scores = scores_rare, negs = negsvec_rare)

make_score_table(full_combined_dataset[MSKFreq == 0],scores = scores_raremsk, negs = negsvec_raremsk)
make_score_table(full_combined_dataset[MSKFreq == 1],scores = scores_raremsk, negs = negsvec_raremsk)
