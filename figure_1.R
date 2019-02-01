#msk_impact = fread("/Users/browna6/mahmood/binom/msk_impact_2017/data_mutations_with38.txt")

if(!exists("rebuilt_nuc")){
  rebuilt_nuc = fread("/Users/browna6/mahmood/binom/analyzed\ data/nucleotide_mutation_table.csv")
}

# reg = ":c.\\d*(.*)"
# 
# msk_impact[,mutation := str_match(HGVSc, reg)[2], by = 1:nrow(msk_impact)]
# rebuilt_nuc[,AminoMutabilityPan := sum(mutabilityPan), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]
# 
# rebuilt_nuc[,AminoCountCosmicPan := sum(CosmicPanCount), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]
# rebuilt_nuc[,AminoMSKCount := sum(MSKCount), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]

cols = c("AminoAcid","mutAA","gene","codon_pos","BinnedRole","AminoCountCosmicPan","AminoMutabilityPan","AminoMSKCount","InMSK")
amino_mutation_table = rebuilt_nuc[, ..cols]
setnames(amino_mutation_table, c("AminoAcid", "mutAA"), c("wildtype", "mutant"))
amino_mutation_table = unique(amino_mutation_table)

amino_mutation_table[AminoCountCosmicPan == 0, ObservationAmino := "Not Observed"]
amino_mutation_table[AminoCountCosmicPan != 0, ObservationAmino := "Observed"]

amino_mutation_table[AminoMSKCount == 0, ObservationAminoMSK := "Not Observed"]
amino_mutation_table[AminoMSKCount != 0, ObservationAminoMSK := "Observed"]


y1 = -4.9
y2 = - 4.7

aminocompareMSK <- ggplot(subset(amino_mutation_table, InMSK == "y"),aes(x = ObservationAminoMSK, y = log10(AminoMutabilityPan),fill = ObservationAminoMSK)) + 
  #geom_violin() +  
  geom_boxplot() +
  scale_fill_brewer(palette="GnBu") +
  xlab(NULL) +
  ylab(NULL) +
  scale_y_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) +
  theme(legend.position = "none") +
  geom_segment(aes(x = 1, y = y2, xend = 2, yend = y2),size = 1) +
  annotate("text", x = 1.5, y = y2 + 0.01, label = "**", size = 5) 

aminocompare <- ggplot(amino_mutation_table,aes(x = ObservationAmino, y = log10(AminoMutabilityPan),fill = ObservationAmino)) + 
  #geom_violin() +  
  geom_boxplot() +
  scale_fill_brewer(palette="Paired") +
  ylab(expression(bold(paste("-Log"[10], "Mutability")))) +
  xlab(NULL) + 
  scale_y_continuous(breaks=c(-6.5, -6,-5.5,-5),labels = c("6.5","6","5.5","5")) + theme(legend.position = "none") +
  geom_segment(aes(x = 1, y = y2, xend = 2, yend = y2),size = 1) +
  annotate("text", x = 1.5, y = y2 + 0.01, label = "**", size = 5) 



allplot <- plot_grid(aminocompare,aminocompareMSK, labels = c("A","B"))


# aM <- wilcox.test(AminoMutabilityPan ~ ObservationAminoMSK, amino_mutation_table)
# n <- wilcox.test(mutabilityPan ~ Observed, rebuilt)
# n <- wilcox.test(mutabilityPan ~ Observed, rebuilt)