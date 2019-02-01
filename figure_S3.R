if(!exists("census_84_allnucleot")){
  census_84_allnucleot<- fread("~/mahmood/binom/analyzed data/census_84_allnucleot_exons_take3.csv")
  
}
#read in a csv containing all the mutations for all the census gens
#read in the census table
if(!exists("allamino_census")){
  allamino_census <- fread( "~/mahmood/binom/analyzed data/cancer_census_all_exons_take3.csv")
  allamino_census_observed = allamino_census[CountNew != 0]
}


subDTAB <- copy(allamino_census_observed[MolecularGenetics == "Dominant"])

mDom <-  make_small_scatter(subtype = "Missense", dtab = subDTAB,measures = c("s","p")) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
nDom <- make_small_scatter(subtype = "Nonsense", dtab = subDTAB,measures = c("s","p")) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
sDom <- make_small_scatter(subtype = "Silent", dtab = subDTAB,measures = c("s","p")) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 

all_dom <- make_small_scatter(dtab = subDTAB,measures = c("s","p"))  +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 

subDTAB <- copy(allamino_census_observed[MolecularGenetics == "Recessive"])

mRec <- make_small_scatter(subtype = "Missense", dtab = subDTAB,measures = c("s","p"))  +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
nRec <- make_small_scatter(subtype = "Nonsense", dtab = subDTAB,measures = c("s","p")) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
sRec <- make_small_scatter(subtype = "Silent", dtab = subDTAB,measures = c("s","p")) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 


all_Rec <- make_small_scatter(dtab = subDTAB,measures = c("s","p"))+ 
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 

#split by Domiant Recess  
domRec <- ggplot(subset(allamino_census, MolecularGenetics %in% c("Dominant","Recessive")),
                 aes(x=MolecularGenetics, y=log10(mutability),fill = SubType, alpha = group)) +
  geom_boxplot() +
  facet_grid(. ~ SubType) +
  scale_fill_manual(values = c("#00BFC4","#F8766D","#7CAE00")) +
  theme(axis.title = element_blank(), legend.position = "none") +
  scale_y_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) 

xvals <- c(1.28125, 1.09375, 0.90625, 0.71875, 2.28125, 2.09375, 1.90625, 1.71875)
xvals <- sort(xvals)
group_labels <- c("0","1","2","3+","0","1","2","3+")


domRec <- annotate_groups(domRec,xvals, yval = -7, fsize = 5, labels = group_labels)
domRec <- annotate_figure(domRec,left = text_grob(expression(bold(paste("-Log"[10], " Mutability"))), rot = 90))



dom <- plot_grid(all_dom, mDom,nDom,sDom, labels = c("A"),nrow = 4,align = "hv")
dom <- annotate_figure(dom,top = text_grob("Dominant"))
rec <- plot_grid(all_Rec,mRec,nRec,sRec, labels = c("B"),nrow = 4,align = "hv")
rec <- annotate_figure(rec,top = text_grob("Recessive"))

bysubRole <- plot_grid(dom,rec,ncol = 2)
bysubRole <- annotate_figure(bysubRole,left = text_grob("Observed Mutation Frequency", rot = 90,face = "bold"), bottom = text_grob("Mutability",face = "bold"))

altogether <- plot_grid(bysubRole,domRec, ncol = 2,rel_widths = c(0.5,0.5),labels = c("","C"))
