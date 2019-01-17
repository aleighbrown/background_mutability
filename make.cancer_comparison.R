cosmic_one = fread("~/mahmood/cosmic_85/cosmic_one_even_cleaner_cp.csv")
rebuilt = fread("~/mahmood/binom/analyzed data/rebuilt_census_table_withmutabilities.csv")

rebuilt[,GenomeMutation := paste0(chromosome_name, ":", genomic_coding,"-",genomic_coding, mutation)]

cosmic_one[, PanCount := .N , by = GenomeMutation]
cosmic_one[`Primary site` == "breast", BreastCount := .N , by = GenomeMutation]
cosmic_one[`Primary site` == "liver", LiverCount := .N , by = GenomeMutation]
cosmic_one[`Primary site` == "prostate", ProstateCount := .N , by = GenomeMutation]
cosmic_one[`Primary site` == "skin" & `Primary histology` == "malignant_melanoma", SkinMelCount := .N , by = GenomeMutation]
cosmic_one[`Primary site` == "lung" & `Primary histology` == "carcinoma", LungCarCount := .N , by = GenomeMutation]
cosmic_one[`Site subtype 1` == "colon" & `Primary histology` == "carcinoma", CoadCount := .N , by = GenomeMutation]



setDT(rebuilt)[cosmic_one, CosmicPanCount := PanCount, on = c(GenomeMutation = "GenomeMutation")]
rebuilt[is.na(CosmicPanCount),CosmicPanCount := 0]

setDT(rebuilt)[cosmic_one, CosmicBreastCount := BreastCount, on = c(GenomeMutation = "GenomeMutation")]
rebuilt[is.na(CosmicBreastCount),CosmicBreastCount := 0]

setDT(rebuilt)[cosmic_one, CosmicLiverCount := LiverCount, on = c(GenomeMutation = "GenomeMutation")]
rebuilt[is.na(CosmicLiverCount),CosmicLiverCount := 0]

setDT(rebuilt)[cosmic_one, CosmicProstateCount := ProstateCount, on = c(GenomeMutation = "GenomeMutation")]
rebuilt[is.na(CosmicProstateCount),CosmicProstateCount := 0]

setDT(rebuilt)[cosmic_one, CosmicSkinMelCount := SkinMelCount, on = c(GenomeMutation = "GenomeMutation")]
rebuilt[is.na(CosmicSkinMelCount),CosmicSkinMelCount := 0]

setDT(rebuilt)[cosmic_one, CosmicLungCarCount := LungCarCount, on = c(GenomeMutation = "GenomeMutation")]
rebuilt[is.na(CosmicLungCarCount),CosmicLungCarCount := 0]

setDT(rebuilt)[cosmic_one, CosmicCoadCount := CoadCount, on = c(GenomeMutation = "GenomeMutation")]
rebuilt[is.na(CosmicCoadCount),CosmicCoadCount := 0]

######amino counts
rebuilt[,AminoCountCosmicPan := sum(CosmicPanCount), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]
rebuilt[,AminoCountCosmicBreast := sum(CosmicBreastCount), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]
rebuilt[,AminoCountCosmicLiver := sum(CosmicLiverCount), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]
rebuilt[,AminoCountCosmicProstate := sum(CosmicProstateCount), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]
rebuilt[,AminoCountCosmicSkinMel := sum(CosmicSkinMelCount), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]
rebuilt[,AminoCountCosmicLungCar := sum(CosmicLungCarCount), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]

######amino mutabilites
rebuilt[,AminoMutabilityCosmicPan := sum(mutabilityPan), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]
rebuilt[,AminoMutabilityCosmicBreast := sum(mutabilityIDC), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]
rebuilt[,AminoMutabilityCosmicLiver := sum(mutabilityLICA), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]
rebuilt[,AminoMutabilityCosmicProstate := sum(mutabilityPRAD), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]
rebuilt[,AminoMutabilityCosmicSkinMel := sum(mutabilitySKCM), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]
rebuilt[,AminoMutabilityCosmicLungCar := sum(mutabilityLUAD), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]


cols = c("mutAA","gene","codon_pos","BinnedRole", names(rebuilt)[grep("Amino",names(rebuilt))])
amino_rebuilt = rebuilt[, ..cols]
setnames(amino_rebuilt, c("AminoAcid", "mutAA"), c("wildtype", "mutant"))
amino_rebuilt = unique(amino_rebuilt)

amino_rebuilt = assign.subtype(amino_rebuilt, "AminoCountCosmicPan", "GroupPan")
amino_rebuilt = assign.subtype(amino_rebuilt, "AminoCountCosmicBreast", "GroupBreast")
amino_rebuilt = assign.subtype(amino_rebuilt, "AminoCountCosmicLiver", "GroupLiver")
amino_rebuilt = assign.subtype(amino_rebuilt, "AminoCountCosmicSkinMel", "GroupSkin")
amino_rebuilt = assign.subtype(amino_rebuilt, "AminoCountCosmicLungCar", "GroupLung")
amino_rebuilt = assign.subtype(amino_rebuilt, "AminoCountCosmicProstate", "GroupProstate")


Breast <- ggplot(amino_rebuilt,aes(x=GroupBreast, y=log10(AminoMutabilityCosmicBreast),fill = SubType)) +
  geom_boxplot() +
  facet_grid(. ~ SubType) +
  theme(axis.title = element_blank(),legend.position = "none") +
  scale_fill_manual(values = c("#00BFC4","#F8766D","#7CAE00")) +
  scale_y_continuous(breaks=c(-7, -6, -5),labels = c("7", "6","5")) +
  #set the size of the axis label
  theme(axis.title = element_text(face = "bold"), axis.title.x = element_blank()) +
  #et the size of the axis tick mark font
  ylab(expression(bold(paste("-Log"[10], "Mutability"["Breast"]))))

Lung <- ggplot(amino_rebuilt,aes(x=GroupLung, y=log10(AminoMutabilityCosmicLungCar),fill = SubType)) +
  geom_boxplot() +
  facet_grid(. ~ SubType) +
  theme(axis.title = element_blank(),legend.position = "none") +
  scale_fill_manual(values = c("#00BFC4","#F8766D","#7CAE00")) +
  scale_y_continuous(breaks=c(-6.5,-5.5,-4.5),labels = c("6.5","5.5","4.5")) +
  #set the size of the axis label
  theme(axis.title = element_text(face = "bold"), axis.title.x = element_blank()) +
  #et the size of the axis tick mark font
  ylab(expression(bold(paste("-Log"[10], "Mutability"["Lung"]))))

Skin <- ggplot(amino_rebuilt,aes(x=GroupSkin, y=log10(AminoMutabilityCosmicSkinMel),fill = SubType)) +
  geom_boxplot() +
  facet_grid(. ~ SubType) +
  theme(axis.title = element_blank(),legend.position = "none") +
  scale_fill_manual(values = c("#00BFC4","#F8766D","#7CAE00")) +
  scale_y_continuous(breaks=c(-6,-5,-4),labels = c("6","5","4")) +
  #set the size of the axis label
  theme(axis.title = element_text(face = "bold"), axis.title.x = element_blank()) +
  #et the size of the axis tick mark font
  ylab(expression(bold(paste("-Log"[10], "Mutability"["Skin"]))))


Liver <- ggplot(amino_rebuilt,aes(x=GroupLiver, y=log10(AminoMutabilityCosmicLiver),fill = SubType)) +
  geom_boxplot() +
  facet_grid(. ~ SubType) +
  theme(axis.title = element_blank(),legend.position = "none") +
  scale_fill_manual(values = c("#00BFC4","#F8766D","#7CAE00")) +
  scale_y_continuous(breaks=c(-7,-6.5,-6,-5.5),labels = c("7","6.5","6","5.5")) +
  #set the size of the axis label
  theme(axis.title = element_text(face = "bold"), axis.title.x = element_blank()) +
  #et the size of the axis tick mark font
  ylab(expression(bold(paste("-Log"[10], "Mutability"["Liver"]))))

Prostate <- ggplot(amino_rebuilt,aes(x=GroupProstate, y=log10(AminoMutabilityCosmicProstate),fill = SubType)) +
  geom_boxplot() +
  facet_grid(. ~ SubType) +
  theme(axis.title = element_blank(),legend.position = "none") +
  scale_fill_manual(values = c("#00BFC4","#F8766D","#7CAE00")) +
  scale_y_continuous(breaks=c(-7,-6.5,-6,-5.5),labels = c("7","6.5","6","5.5")) +
  #set the size of the axis label
  theme(axis.title = element_text(face = "bold"), axis.title.x = element_blank()) +
  #et the size of the axis tick mark font
  ylab(expression(bold(paste("-Log"[10], "Mutability"["Prostate"]))))

full_plot = ggarrange(Breast, Lung, Skin, Liver, nrow = 2, ncol = 2, labels = c("A","B","C","D"))

quartz()
print(full_plot)
