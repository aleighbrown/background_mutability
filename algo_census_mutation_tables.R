#load the necessary functions
source('~/mahmood/binom/get_transcript_sequences.funcs.R')
#read in the transcripts
principal_transcripts = fread("~/mahmood/binom/analyzed data/principal_transcripts_census_genes.csv")
#get the exonIds 
exonIds = get_exonIds(principal_transcripts[,V4])
setDT(principal_transcripts)[exonIds, ChosenTranscript := ensembl_transcript_id, on = c(V2 = "ensembl_gene_id")]
#get the CDNA for the exons
exonswithCDNA = get_nucSequence(exonIds)
#get the mutation table
mutation_table = create_all_transcript_coding_tables(exonswithCDNA)
#add the hugo name back on 
setDT(mutation_table)[principal_transcripts, gene := i.hgnc, on = c(ensembl_transcript_id = "ChosenTranscript")]
setDT(mutation_table)[principal_transcripts, ensembl_gene_id := V2, on = c(ensembl_transcript_id = "ChosenTranscript")]
#add the codon positions
mutation_table[, codon_pos := return_codon_position_vector(pos), by = ensembl_transcript_id]
#add the codons
codons_table = create_full_cdna_table(exonswithCDNA)
setDT(mutation_table)[codons_table, codon := i.x, on = c(ensembl_transcript_id = "ensembl_transcript_id", codon_pos = "codon_pos")]
setDT(mutation_table)[codons_table, AminoAcid := i.amino, on = c(ensembl_transcript_id = "ensembl_transcript_id", codon_pos = "codon_pos")]

fwrite(mutation_table, "~/mahmood/binom/analyzed data/rebuilt_census_table.csv")
#get the signatures



mutation_table = append_mutability(mutation_table)
mutation_table = append_mutability(mutation_table, 9026, "LUAD")
mutation_table = append_mutability(mutation_table, 9032, "IDC")
mutation_table = append_mutability(mutation_table, 9033, "COAD")
mutation_table = append_mutability(mutation_table, 9037, "PRAD")
mutation_table = append_mutability(mutation_table, 9038, "LICA")
mutation_table = append_mutability(mutation_table, 9049, "BLCA")

mutation_table = append_mutability(mutation_table, 9013, "PACA")
mutation_table = append_mutability(mutation_table, 9025, "SKCM")


mutation_table = mutate.codons(mutation_table)

mutation_table[AminoAcid == mutAA & AminoAcid != "*", SubType := "Silent"]
mutation_table[AminoAcid != mutAA & mutAA != "*" & AminoAcid != "*", SubType := "Missense"]
mutation_table[mutAA == "*" & AminoAcid != "*", SubType := "Nonsense"]

census_genes = fread("~/mahmood/binom/census_genes_7_21_2018.csv")
setDT(mutation_table)[census_genes, RoleinCancer := `Role in Cancer`, on = c(gene = "SuggestedSymbol")]
setDT(mutation_table)[census_genes, MolecularGenetics := `Molecular Genetics`, on = c(gene = "SuggestedSymbol")]
setDT(mutation_table)[census_genes, BinnedRole := `BinnedRole`, on = c(gene = "SuggestedSymbol")]

mutation_table[,mutation := paste0(wildtype,">",mutant)]
mutation_table = mutation_table[gene != ""]
fwrite(mutation_table[!is.na(gene)], "~/mahmood/binom/analyzed data/rebuilt_census_table_withmutabilities.csv")

cosmic_one = fread("~/mahmood/cosmic_85/cosmic_one_even_cleaner_cp.csv")

#mutation_table = fread("~/mahmood/binom/analyzed data/rebuilt_census_table_withmutabilities.csv")

mutation_table[,GenomeMutation := paste0(chromosome_name, ":", genomic_coding,"-",genomic_coding, mutation)]

cosmic_one[, PanCount := .N , by = GenomeMutation]
cosmic_one[`Primary site` == "breast", BreastCount := .N , by = GenomeMutation]
cosmic_one[`Primary site` == "liver", LiverCount := .N , by = GenomeMutation]
cosmic_one[`Primary site` == "prostate", ProstateCount := .N , by = GenomeMutation]

setDT(mutation_table)[cosmic_one, CosmicPanCount := PanCount, on = c(GenomeMutation = "GenomeMutation")]
mutation_table[is.na(CosmicPanCount),CosmicPanCount := 0]
setDT(mutation_table)[cosmic_one, CosmicBreastCount := BreastCount, on = c(GenomeMutation = "GenomeMutation")]
mutation_table[is.na(CosmicBreastCount),CosmicBreastCount := 0]
setDT(mutation_table)[cosmic_one, CosmicLiverCount := LiverCount, on = c(GenomeMutation = "GenomeMutation")]
mutation_table[is.na(CosmicLiverCount),CosmicLiverCount := 0]
setDT(mutation_table)[cosmic_one, CosmicProstateCount := ProstateCount, on = c(GenomeMutation = "GenomeMutation")]
mutation_table[is.na(CosmicProstateCount),CosmicProstateCount := 0]

mutation_table[,AminoCountCosmicPan := sum(CosmicPanCount), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]
mutation_table[,AminoCountCosmicBreast := sum(CosmicBreastCount), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]
mutation_table[,AminoCountCosmicLiver := sum(CosmicLiverCount), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]
mutation_table[,AminoCountCosmicProstate := sum(CosmicProstateCount), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]

mutation_table[,AminoMutabilityLICA := sum(mutabilityLICA), by = c("gene", "AminoAcid", "mutAA", "codon_pos")]

cols = c("AminoAcid","mutAA","gene","codon_pos","BinnedRole","AminoCountPan","AminoMutabilityPan","AminoMutabilityIDC", "AminoMutabilityPRAD", "AminoMutabilityLICA", "AminoCountCosmicBreast", "AminoCountCosmicLiver", "AminoCountCosmicProstate","AminoCountCosmicPan")
amino_mutation_table = mutation_table[, ..cols]
setnames(amino_mutation_table, c("AminoAcid", "mutAA"), c("wildtype", "mutant"))
amino_mutation_table = unique(amino_mutation_table)

amino_mutation_table = assign.subtype(amino_mutation_table, "AminoCountPan", "GroupPanChang")
amino_mutation_table = assign.subtype(amino_mutation_table, "AminoCountCosmicPan", "GroupPan")
amino_mutation_table = assign.subtype(amino_mutation_table, "AminoCountCosmicBreast", "GroupBreast")
amino_mutation_table = assign.subtype(amino_mutation_table, "AminoCountCosmicLiver", "GroupLiver")
amino_mutation_table = assign.subtype(amino_mutation_table, "AminoCountCosmicProstate", "GroupProstate")

