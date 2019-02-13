# if(!exists("allamino_census_observed")){
#   #read in the census if I haven't yet
#   allamino_census_observed <- fread("~/mahmood/binom/analyzed data/allamino_census_observed.csv")
# }
# 
# if(!exists("corr_table")){
#   #read in the census if I haven't yet
#   corr_table <- fread("~/mahmood/binom/analyzed data/corr_table.csv")
# }
#read in the census table
if(!exists("allamino_census")){
  allamino_census <- fread( "/Users/browna6/mahmood/binom/analyzed\ data/cancer_census_all_exons_take3.csv")
}

allamino_census_observed = allamino_census[AminoCountCosmicPan !=0]


if(!exists("corr_table")){
  #read in the census if I haven't yet
  corr_table <- fread("~/mahmood/binom/analyzed data/corr_table_newcensus.csv")
}
# # 

alph = 0.1
vj = 1
fs = 12
#make a histogram of the correlations and color them differently by whether the correlation is significant or not
hstS <- ggplot(subset(corr_table,SubType == "Silent"),aes(x=`Spearman Rho`))+
  geom_histogram(data = subset(corr_table,SubType == "Silent" & Sig == "Y"), color = "black",fill = "#7CAE00") + 
  geom_histogram(data = subset(corr_table,SubType == "Silent" & Sig == "N"), color = "black",alpha = alph, fill = "#7CAE00") + 
  theme(legend.position = "none",axis.title = element_blank()) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) +
  coord_cartesian(xlim = c(-0.65, 0.825)) +
  theme(axis.text = element_text(face = "plain")) +
  annotation_custom(grobTree(textGrob("Silent  ", x=0.03,  y=0.95, hjust=0,vjust = vj,
                    gp=gpar(fontsize=fs))))


hstN <- ggplot(subset(corr_table,SubType == "Nonsense"),aes(x=`Spearman Rho`))+
  geom_histogram(data = subset(corr_table,SubType == "Nonsense" & Sig == "Y"), color = "black",fill = "#F8766D") + 
  geom_histogram(data = subset(corr_table,SubType == "Nonsense" & Sig == "N"), color = "black",alpha = alph, fill = "#F8766D") + 
  theme(legend.position = "none",axis.title = element_blank()) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) +
  coord_cartesian(xlim = c(-0.65, 0.825)) +
  theme(axis.text = element_text(face = "plain")) +
  scale_y_continuous(breaks = c(0,3,6,9)) +
  annotation_custom(grobTree(textGrob("Nonsense", x=0.03,  y=0.95, hjust=0,vjust = vj,
                            gp=gpar(fontsize=fs))))

#so I'm going to take the top 2 genes for each subtype and make a small scatter for each
hstM <- ggplot(subset(corr_table,SubType == "Missense"),aes(x=`Spearman Rho`))+
  geom_histogram(data = subset(corr_table,SubType == "Missense" & Sig == "Y"), color = "black",fill = "#00BFC4") + 
  geom_histogram(data = subset(corr_table,SubType == "Missense" & Sig == "N"), color = "black",alpha = alph, fill = "#00BFC4") + 
  theme(legend.position = "none",axis.title = element_blank()) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) +
  coord_cartesian(xlim = c(-0.65, 0.825))+
  theme(axis.text = element_text(face = "plain")) +
  annotation_custom(grobTree(textGrob("Missense  ", x=0.03,  y=0.95, hjust=0,vjust = vj,
                                      gp=gpar(fontsize=fs))))


# gfont = 20
# rfont = 20
# rvj = 0.1
corr_table[Sig == "Y",.SD[1:2], by=SubType]

#silent
#make_small_scatter <- function(gene = "none",subtype,fs = 20,Rvj = 0.1, RFont = 20,dtab,fontface = "bold.italic",frequencyType = "observed_wgs",measures)
sil_one <- make_small_scatter(gene = "PTPRT",subtype = "Silent",dtab = allamino_census_observed,measures = "l",frequencyType = "observed_wgs") +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
sil_two <- make_small_scatter(gene = "POLD1",subtype = "Silent",dtab = allamino_census_observed,measures = "l",frequencyType = "observed_wgs") +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
#non 
non_one <- make_small_scatter(gene = "TP53",subtype = "Nonsense",dtab = allamino_census_observed,fontface = "plain",measures = "l",frequencyType = "observed_wgs") +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
non_two<- make_small_scatter(gene = "CASP8",subtype = "Nonsense",dtab = allamino_census_observed,measures = "l",frequencyType = "observed_wgs") +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
#miss
miss_one <- make_small_scatter(gene = "RSPO2",subtype = "Missense",dtab = allamino_census_observed,measures = "l",frequencyType = "observed_wgs") +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
miss_two <- make_small_scatter(gene = "ARHGAP26",subtype = "Missense",dtab = allamino_census_observed,measures = "l",frequencyType = "observed_wgs") +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))


### make the scatters
sil <- plot_grid(sil_one, sil_two,ncol = 2, labels = c("F",""))
non <- plot_grid(non_one, non_two, ncol = 2, labels = c("E", ""))
mis <- plot_grid(miss_one, miss_two, ncol = 2, labels = c("D",""))


scat <- plot_grid(mis,non,sil, ncol = 1)

#scat <- annotate_figure(scat, left = text_grob("Observed Mutation Frequency", rot = 90,size = 20,face = "bold"), bottom = text_grob("Mutability",size = 20, face = "bold"))
scat <- annotate_figure(scat, left = text_grob("Observed Mutations", rot = 90,face = "bold"), bottom = text_grob("Mutability", face = "bold"))
scat <- scat + theme(plot.margin = unit(c(0,0,0,2), "cm"))

#make the combined plot for silent mutations
hists <- plot_grid(hstM,hstN,hstS, ncol = 1,align = "vh",labels = c("A","B","C")) 

hists <- annotate_figure(hists, left = text_grob("Number of Genes", rot = 90,face = "bold")) + theme(plot.margin = unit(c(0,0,0.1,0.1), "cm"))

#make a pretty table to draw on
plotintable <- corr_table[Sig == "Y",c("Gene","SubType","Spearman Rho","RoleinCancer","BinnedRole")]
#setorder(plotintable,-"Spearman Rho")
setnames(plotintable,"Spearman Rho","Rho")

plotintable[,GeneSub := paste(Gene,SubType)]
plotintable$GeneSub <- factor(plotintable$GeneSub, levels = plotintable$GeneSub[order(plotintable$Rho)])


plotintable <- plotintable[, head(.SD, 10), by=SubType]


small_miss <- plotintable[SubType == "Missense"][order(Rho)][,c("Gene","BinnedRole")]
miss_labs <- label_vector(small_miss)

miss_rho <- ggplot(subset(plotintable, SubType == "Missense"), aes(y = Rho, x = GeneSub)) + 
  geom_col( aes(fill = SubType)) + 
  coord_flip() + 
  scale_y_continuous(limits = c(0,0.802)) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) +
  scale_x_discrete(labels = miss_labs) + 
  scale_fill_manual(values = "#00BFC4") + 
  theme(legend.position = "none") +
  theme(axis.title = element_blank()) 

small_non <- plotintable[SubType == "Nonsense"][order(Rho)][,c("Gene","BinnedRole")]
non_labs <- label_vector(small_non)

non_rho <- ggplot(subset(plotintable, SubType == "Nonsense"), aes(y = Rho, x = GeneSub)) + 
  geom_col(aes(fill = SubType)) + 
  coord_flip() + 
  scale_y_continuous(limits = c(0,0.802)) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) +
  scale_x_discrete(labels = non_labs) + 
  scale_fill_manual(values = "#F8766D") +
  theme(legend.position = "none") +
  theme(axis.title = element_blank()) 

small_non <- plotintable[SubType == "Silent"][order(Rho)][,c("Gene","BinnedRole")]
sil_labs <- label_vector(small_non)

silent_rho <- ggplot(subset(plotintable, SubType == "Silent"), aes(y = Rho, x = GeneSub)) + 
  geom_col(aes(fill = SubType)) +
  coord_flip() + 
  scale_y_continuous(limits = c(0,0.802)) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) +
  scale_x_discrete(labels = sil_labs) + 
  scale_fill_manual(values = "#7CAE00") +
  theme(legend.position = "none") +
  theme(axis.title = element_blank()) 


rho_plot <- plot_grid(miss_rho,non_rho,silent_rho,nrow = 3) 


fig2 <- plot_grid(hists, rho_plot, nrow = 1,align = 'vh') 

fig2 <- annotate_figure(fig2, bottom = text_grob("Correlation Coefficient", face = "bold"))

plot_grid(fig2,scat,nrow = 1, ncol = 2, rel_widths = c(0.5,0.5)) 
