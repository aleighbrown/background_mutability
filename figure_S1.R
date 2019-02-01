


#read in the census table
if(!exists("census_84_allnucleot")){
  census_84_allnucleot<- fread("/Users/browna6/mahmood/binom/analyzed\ data/nucleotide_mutation_table.csv")
  
}

all_census <- make_standardized_values_nuc_mut(census_84_allnucleot)
all_census[,mutability := mutability * 10^6]

#census_84_allnucleot[,Mut2 := paste0(substr(trinuc,1,1),"[",wildtype,">",mutant,"]",substr(trinuc,3,3))]
mutability_cdf <- ggplot(all_census,aes(x = mutability,color = Measure)) + stat_ecdf(size = 2) + 
  scale_color_manual(values = c("#3366ff","#e2940d")) + 
  theme(axis.text = element_text(size = 17)) + 
  theme(legend.position = c(0.5,0.5)) +
  labs(color = "Mutability Spectrum", x =expression(bold("-Log"[10], "Mutability")),y = bquote(bold("Cumulative Probability"))) + 
  guides(color = guide_legend(override.aes = list(size=10))) 




labelsAll <- unique(all_census[Measure == "Observed",c("mutability","Prop","PlotMutation")])


backtoback <- ggplot(all_census, aes(x = mutability)) + 
  geom_histogram(data = subset(all_census, Measure == "Observed"), aes(y = ..count../sum(..count..)), bins = 115, fill= "#e2940d",color = "black") + 
  geom_histogram(data = subset(all_census, Measure == "All Mutations"), aes(y = -..count../sum(..count..)), bins = 115, fill= "#3366ff",color = "black") +
  #scale_x_reverse() + 
  #coord_flip(ylim = c(-0.16,0.08)) +
  scale_y_continuous(breaks = c(-0.15,-0.075,0,0.075),labels = c("0.15","0.075","0","0.075")) +
  theme(axis.title = element_blank()) + 
  theme(plot.margin = unit(c(0,0,0,0.5), "cm"))

#annotate the mutations
backtoback <- backtoback + annotate("text", y = labelsAll[1:6,Prop] + 0.0017, x = labelsAll[1:6,mutability] - .18, label = labelsAll[1:6,PlotMutation])
#add the cdf
backtoback <- backtoback + annotation_custom(ggplotGrob(mutability_cdf),xmin = 5,xmax = 11, ymax = -0.019,ymin = -0.16)

casp8_muts <- make_standardized_values_nuc_mut(census_84_allnucleot[gene == "CASP8"])
casp8_muts[,mutability := mutability * 10^6]
tp53_muts <- make_standardized_values_nuc_mut(census_84_allnucleot[gene == "TP53"])
tp53_muts[,mutability := mutability * 10^6]

casp_hist <- ggplot(casp8_muts, aes(x = mutability)) + 
  geom_histogram(data = subset(casp8_muts, Measure == "Observed"), aes(y = ..count../sum(..count..)), bins = 115, fill= "#e2940d",color = "black") + 
  geom_histogram(data = subset(casp8_muts, Measure == "All Mutations"), aes(y = -..count../sum(..count..)), bins = 115, fill= "#3366ff",color = "black") +
  theme(axis.title = element_blank(),axis.text = element_text(size = 17)) +
  coord_cartesian(ylim = c(-.16,0.16)) +
  scale_y_continuous(breaks = c(-0.1,0,0.1),labels = c("0.1","0","0.1")) 
  

tp53_hist <- ggplot(tp53_muts, aes(x = mutability)) + 
  geom_histogram(data = subset(tp53_muts, Measure == "Observed"), aes(y = ..count../sum(..count..)), bins = 115, fill= "#e2940d",color = "black") + 
  geom_histogram(data = subset(tp53_muts, Measure == "All Mutations"), aes(y = -..count../sum(..count..)), bins = 115, fill= "#3366ff",color = "black") +
  theme(axis.title = element_blank(),axis.text = element_text(size = 17)) +
  coord_cartesian(ylim = c(-.16,0.16)) +
scale_y_continuous(breaks = c(-0.1,0,0.1),labels = c("0.1","0","0.1")) 
  
gene_hists <- plot_grid(casp_hist,tp53_hist,nrow = 2,labels = c("B","C")) + theme(plot.margin = unit(c(0,0.5,0,0), "cm")) 

final_observedexpected <- plot_grid(backtoback,gene_hists,ncol = 2, labels = "A", rel_widths = c(2,1))

final_observedexpected <- annotate_figure(final_observedexpected, bottom = text_grob(bquote(bold('Mutability *'~10^6))), 
                   left = text_grob("Proportion",rot = 90,face = "bold"))


