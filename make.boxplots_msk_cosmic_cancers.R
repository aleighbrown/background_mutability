rebuilt = fread("/Users/browna6/mahmood/binom/analyzed\ data/rebuilt_census_table_withmutabilities.csv")

msk_impact = fread("/Users/browna6/mahmood/binom/msk_impact_2017/data_mutations_with38.txt")

make_small_melty_table <- function(dtab, group_columns = c("CosmicBreastGroup", "IDCGroup"), mut_column ="mutabilityIDC"){
  chosen_columns = c(group_columns, mut_column)
  temp = dtab[,..chosen_columns]
  temp[,rowid := .I]
  temp_m = melt(temp, id.vars = c("rowid",mut_column), value.var = group_columns)
  return(temp_m)
}

breast = make_small_melty_table(rebuilt[InMSK == "y"], group_columns = c("CosmicBreastGroup", "IDCGroup"), mut_column ="mutabilityIDC")
lung = make_small_melty_table(rebuilt[InMSK == "y"], group_columns = c("CosmicLungGroup", "LUADGroup"), mut_column ="mutabilityLUAD")
coad = make_small_melty_table(rebuilt[InMSK == "y"],  group_columns = c("CosmicCoadGroup", "COADGroup"), mut_column ="mutabilityCOAD")
skin = make_small_melty_table(rebuilt[InMSK == "y"],  group_columns = c("CosmicSkinGroup", "SKCMGroup"), mut_column ="mutabilitySKCM")

breast_plot = ggplot(breast, aes(x = value, y = log10(mutabilityIDC), fill = variable)) + geom_boxplot() +
  theme(axis.title.x = element_blank(),legend.position = "none") +
  scale_fill_manual(values = c("#3C77AF","#B3DBB8")) + 
  ylab(expression(bold(paste("-Log"[10], "Mutability"["Breast"]))))

lung_plot = ggplot(lung, aes(x = value, y = log10(mutabilityLUAD), fill = variable)) + geom_boxplot() +
  theme(axis.title.x = element_blank(),legend.position = "none") +
  scale_fill_manual(values = c("#3C77AF","#B3DBB8")) +
  ylab(expression(bold(paste("-Log"[10], "Mutability"["Lung"]))))

colon_plot = ggplot(coad, aes(x = value, y = log10(mutabilityCOAD), fill = variable)) + geom_boxplot() +
  theme(axis.title.x = element_blank(),legend.position = "none") + 
  scale_fill_manual(values = c("#3C77AF","#B3DBB8")) +
  ylab(expression(bold(paste("-Log"[10], "Mutability"["Colon"]))))

skin_plot = ggplot(skin, aes(x = value, y = log10(mutabilitySKCM), fill = variable)) + geom_boxplot() +
  theme(axis.title.x = element_blank(),legend.position = "none") +
  scale_fill_manual(values = c("#3C77AF","#B3DBB8")) +
  ylab(expression(bold(paste("-Log"[10], "Mutability"["Skin"]))))

ggarrange(breast_plot, lung_plot, colon_plot, skin_plot, nrow = 2, ncol = 2, labels = c("A","B","C","D"))
