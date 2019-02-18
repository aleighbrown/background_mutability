##read in the census table
if(!exists("allamino_census_observed")){
  allamino_census <- fread("/Users/browna6/mahmood/binom/analyzed\ data/codon_substitution_table.csv")
  allamino_census_observed = allamino_census[AminoCountCosmicPan !=0]
}

#all genes and all subtypes
all <-  make_small_scatter(dtab = allamino_census_observed,measures = c("s","p")) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
all <- annotate_figure(all,left = text_grob("All Substitutions", rot = 90,size = 11))
#all subtypes only oncogenes
all_onco <- make_small_scatter(dtab = allamino_census_observed[BinnedRole == "Oncogene"],measures = c("s","p"))+ 
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
#all subtypes only TSG
all_TSG <- make_small_scatter(dtab = allamino_census_observed[BinnedRole == "TSG"],measures = c("s","p"))+ 
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 


mOnc <- make_small_scatter(gene = "none",subtype = "Missense",Rvj = 0.1, dtab = allamino_census_observed[BinnedRole == "Oncogene"],frequencyType = "observed_wgs",measures = c("s","p")) + 
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
nOnc <- make_small_scatter(gene = "none",subtype = "Nonsense",Rvj = 0.1, dtab = allamino_census_observed[BinnedRole == "Oncogene"],frequencyType = "observed_wgs",measures = c("s","p"))+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
sOnc <- make_small_scatter(gene = "none",subtype = "Silent",Rvj = 0.1, dtab = allamino_census_observed[BinnedRole == "Oncogene"],frequencyType = "observed_wgs",measures = c("s","p"))+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 

mTSG <- make_small_scatter(gene = "none",subtype = "Missense",Rvj = 0.1, dtab = allamino_census_observed[BinnedRole == "TSG"],frequencyType = "observed_wgs",measures = c("s","p")) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
nTSG <- make_small_scatter(gene = "none",subtype = "Nonsense",Rvj = 0.1, dtab = allamino_census_observed[BinnedRole == "TSG"],frequencyType = "observed_wgs",measures = c("s","p"))+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
sTSG <- make_small_scatter(gene = "none",subtype = "Silent",Rvj = 0.1, dtab = allamino_census_observed[BinnedRole == "TSG"],frequencyType = "observed_wgs",measures = c("s","p")) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 

onco <- plot_grid(all_onco, mOnc,nOnc,sOnc, labels = c("B") ,nrow = 4,align = "hv")
onco <- annotate_figure(onco, top = text_grob("Oncogenes",size = 18))

tsg <- plot_grid(all_TSG, mTSG,nTSG,sTSG, labels = c("C"), nrow = 4,align = "hv")
tsg <- annotate_figure(tsg, top = text_grob("TSG",size = 18))

nonsense <- make_small_scatter(subtype = "Nonsense", dtab = allamino_census_observed,measures = c("s","p")) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
nonsense <- annotate_figure(nonsense,left = text_grob("Nonsense", rot = 90,size = 11))

silent <- make_small_scatter(subtype = "Silent", dtab = allamino_census_observed,measures = c("s","p")) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
silent <- annotate_figure(silent,left = text_grob("Silent", rot = 90,size = 11))

missense <- make_small_scatter(subtype = "Missense", dtab = allamino_census_observed,measures = c("s","p")) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
missense <- annotate_figure(missense,left = text_grob("Missense", rot = 90,size = 11))


bySub <- plot_grid(all, missense,nonsense,silent,nrow = 4, labels = c("A"), align = "hv")
bySub <- annotate_figure(bySub, top =text_grob("All Genes",size = 18))
bysubRole <- plot_grid(bySub,onco,tsg,ncol = 3) + theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
annotate_figure(bysubRole,left = text_grob("Observed Mutation", rot = 90,size = 12,face = "bold"), bottom = text_grob(expression(bold("Mutability"))))
