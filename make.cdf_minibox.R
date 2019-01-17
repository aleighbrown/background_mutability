#read in the census table
if(!exists("allamino_census")){
  allamino_census <- fread("~/mahmood/binom/analyzed data/cancer_census_all_exons_take5.csv")
}

allamino_census_observed = allamino_census[AminoCountCosmicPan != 0 & SubType %in% c("Missense","Nonsense", "Silent")]

#make the boxplot for all
all <- ggplot(allamino_census,aes(x=group, y=log10(AminoMutabilityPan),fill = SubType)) +
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
type_cdf <- ggplot(allamino_census_observed, aes(x = log10(AminoMutabilityPan), color = SubType)) + 
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
  #lims(x = c(0,max(log10(allamino_census$mutability)))) +
  #make the ecdf and make the size of the lines a bit thicker
  stat_ecdf(size = 3) +
  scale_x_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) +
  #forces the legend samples to be much larger
  guides(colour = guide_legend(override.aes = list(size=5))) + theme(axis.title.x = element_blank())

#this is for where to draw the lines
y1 = -4.9
y2 = - 4.78

mbp <- ggplot(allamino_census_observed, aes(y = log10(AminoMutabilityPan),x = SubType, fill = SubType)) +
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
allamino_census_observed$SubType <- as.factor(allamino_census_observed$SubType)
kruskal.test(allamino_census_observed$AminoMutabilityPan, allamino_census_observed$SubType)
f <- dunn.test::dunn.test(allamino_census_observed$AminoMutabilityPan, allamino_census_observed$SubType)
