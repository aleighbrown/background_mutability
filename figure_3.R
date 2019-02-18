#read in the census table
if(!exists("rebuilt_amino")){
  rebuilt_amino <- fread("~/mahmood/binom/analyzed data/codon_substitution_table.csv")
}

rebuilt_amino_observed = rebuilt_amino[AminoCountCosmicPan != 0 & SubType %in% c("Missense","Nonsense", "Silent")]

#make the boxplot for all
all <- ggplot(subset(rebuilt_amino,SubType %in% c("Missense","Nonsense", "Silent")), aes(x=group, y=log10(AminoMutabilityPan),fill = SubType)) +
  geom_boxplot() +
  facet_grid(. ~ SubType) +
  theme(axis.title = element_blank(),legend.position = "none") +
  scale_fill_manual(values = c("#00BFC4","#F8766D","#7CAE00")) +
  #make the labes be positive
  scale_y_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) +
  #set the size of the axis label
  theme(axis.title = element_text(face = "bold")) +
  #et the size of the axis tick mark font
  theme(axis.title.x = element_blank()) +
  ylab(expression(bold(paste("-Log"[10], "Mutability"))))  




#make a cdf plot of all the possible mutations splitting it by the SubType either missense, nonsense or silent
type_cdf <- ggplot(rebuilt_amino_observed, aes(x = log10(AminoMutabilityPan), color = SubType)) + 
  #set the color manuall so that miss == Blue, Nonsense == Red, and silent == green
  scale_color_manual(values = c("#00BFC4","#F8766D","#7CAE00")) + 
  #set the size of the axis label
  theme(axis.title = element_text(face = "bold")) +
  #et the size of the axis tick mark font
  ylab("Cumulative Probability") +
  #xlab(expression(bold(paste("-Log"[10], "Mutability")))) +
  #put the legend in the bottom right corner and set the size of the text
  theme(
    legend.position = c(.95, .25),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.title = element_blank()) +
  #change the legend title
  #lims(x = c(0,max(log10(rebuilt_amino$mutability)))) +
  #make the ecdf and make the size of the lines a bit thicker
  stat_ecdf(size = 3) +
  scale_x_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) +
  #forces the legend samples to be much larger
  guides(colour = guide_legend(override.aes = list(size=5))) + theme(axis.title.x = element_blank())

#this is for where to draw the lines
y1 = -4.9
y2 = - 4.78

mbp <- ggplot(rebuilt_amino_observed, aes(y = log10(AminoMutabilityPan),x = SubType, fill = SubType)) +
  geom_violin() +
  scale_fill_manual(values = c("#00BFC4","#F8766D","#7CAE00")) +
  theme(legend.position="none",axis.text.x=element_blank(), axis.title.x = element_blank()) +
  #betwen silent and missense
  geom_segment(aes(x = 1, y = y2, xend = 3, yend = y2),size = 1) +
  annotate("text", x = 2, y = y2 + 0.01, label = "**", size = 5) + 
  #between missense and nonsense
  geom_segment(aes(x = 1.1, y = y1, xend = 1.9, yend = y1),size = 1) + 
  annotate("text", x = 1.5, y = y1 + 0.01, label = "**", size = 5) +
  #between nonsense and silent
  geom_segment(aes(x = 2.1, y = y1, xend = 2.9, yend = y1),size = 1) + 
  annotate("text", x = 2.5, y = y1 + 0.01, label = "**", size = 5) +
  
  ylab(expression(bold(paste("-Log"[10], "Mutability")))) +
  xlab(NULL) + 
  scale_y_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) + theme(aspect.ratio = 1)

  

#add the mutability boxplot as a custom grob
type_cdf <-  type_cdf + annotation_custom(ggplotGrob(mbp),xmin = -6.9, xmax = -5.8,ymin = 0.35)

aminoacid_plot <- plot_grid(type_cdf,all,align = "h",axis = 'b',labels = c("A","B"))

# ###for tests
#correlation between mutability and frequency
rebuilt_amino_observed$SubType <- as.factor(rebuilt_amino_observed$SubType)
kruskal.test(rebuilt_amino_observed$AminoMutabilityPan, rebuilt_amino_observed$SubType)
f <- dunn.test::dunn.test(rebuilt_amino_observed$AminoMutabilityPan, rebuilt_amino_observed$SubType)

#read in the census table
if(!exists("rebuilt_nuc")){
  rebuilt_nuc = fread("/Users/browna6/mahmood/binom/analyzed\ data/nucleotide_mutation_table.csv")

}

rebuilt_nuc_observed <- rebuilt_nuc[CosmicPanCount != 0 & SubType %in% c("Missense","Nonsense", "Silent")]


#make the boxplot for all
all_nuc <- ggplot(subset(rebuilt_nuc,SubType %in% c("Missense","Nonsense", "Silent")),aes(x=CosmicPanGroup, y=log10(mutabilityPan),fill = SubType)) +
  geom_boxplot() +
  facet_grid(. ~ SubType) +
  theme(axis.title = element_blank(),legend.position = "none") +
  scale_fill_manual(values = c("#00BFC4","#F8766D","#7CAE00")) +
  scale_y_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) +
  #set the size of the axis label
  theme(axis.title = element_text(face = "bold")) +
  #et the size of the axis tick mark font
  xlab("Observed Mutation Frequency") +
  ylab(expression(bold(paste("-Log"[10], "Mutability"))))  



#make a cdf plot of all the observed mutations splitting it by the SubType either missense, nonsense or silent
type_cdf_nuc <- ggplot(rebuilt_nuc_observed, aes(x = log10(mutabilityPan), color = SubType)) + 
  #set the color manuall so that miss == Blue, Nonsense == Red, and silent == green
  scale_color_manual(values = c("#00BFC4","#F8766D","#7CAE00")) + 
  #set the size of the axis label
  theme(axis.title = element_text(face = "bold")) +
  #et the size of the axis tick mark font
  ylab("Cumulative Probability") +
  xlab(expression(bold(paste("-Log"[10], "Mutability")))) +
  theme(legend.position = "none") +
  #lims(x = c(0,max(log10(rebuilt_amino$mutability)))) +
  #make the ecdf and make the size of the lines a bit thicker
  stat_ecdf(size = 3) +
  scale_x_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) 


y1 = -4.9
y2 = - 4.7

mbp_nuc <- ggplot(rebuilt_nuc_observed, aes(y = log10(mutabilityPan),x = SubType, fill = SubType)) +
  geom_violin() +
  scale_fill_manual(values = c("#00BFC4","#F8766D","#7CAE00")) +
  theme(legend.position="none",axis.text.x=element_blank(), axis.title.x = element_blank()) +
  #betwen silent and missense
  geom_segment(aes(x = 1, y = y2, xend = 3, yend = y2),size = 1) +
  annotate("text", x = 2, y = y2 + 0.01, label = "**", size = 5) + 
  #between missense and nonsense
  geom_segment(aes(x = 1.1, y = y1, xend = 1.9, yend = y1),size = 1) + 
  annotate("text", x = 1.5, y = y1 + 0.01, label = "**", size = 5) +
  #between nonsense and silent
  geom_segment(aes(x = 2.1, y = y1, xend = 2.9, yend = y1),size = 1) + 
  annotate("text", x = 2.5, y = y1 + 0.01, label = "**", size = 5) +
  
  ylab(expression(bold(paste("-Log"[10], "Mutability")))) +
  xlab(NULL) + 
  scale_y_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) + theme(aspect.ratio = 1)


#add the boxplot as a custom grob
type_cdf_nuc <-  type_cdf_nuc + annotation_custom(ggplotGrob(mbp_nuc),xmin = -6.9, xmax = -5.8,ymin = 0.35)

nucleotide_plot <- plot_grid(type_cdf_nuc,all_nuc,align = "h",axis = 'b',labels = c("C","D"))

ggarrange(aminoacid_plot,nucleotide_plot, nrow = 2)
# ###for tests

# rebuilt_nuc_observed$SubType <- as.factor(rebuilt_nuc_observed$SubType)
# kruskal.test(rebuilt_nuc_observed$mutabilityPan, rebuilt_nuc_observed$SubType)
# test_result = dunn.test::dunn.test(rebuilt_nuc_observed$mutabilityPan, rebuilt_nuc_observed$SubType)
