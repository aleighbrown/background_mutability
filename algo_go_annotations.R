library(biomaRt)
library(GOfuncR)
library(data.table)

get_go_ids <- function(gene_list){
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  results <- as.data.table(getBM(attributes = c("hgnc_symbol", "go_id","name_1006","namespace_1003","ensembl_gene_id"), filters = "ensembl_gene_id", values = gene_list, mart = mart))
  #concate the name and the go id using pipe operator, since that's a nice one to delimit on imo
  results[,NameGo := paste(go_id,name_1006, sep = "|")]
  #aggregate the go ids for each name space
  full_results = as.data.table(aggregate(NameGo ~ hgnc_symbol + namespace_1003 + ensembl_gene_id, results, FUN=toString))
  full_results = full_results[namespace_1003 != ""]
  full_results = dcast(full_results,hgnc_symbol + ensembl_gene_id ~ namespace_1003, value.var = "NameGo")
  return(full_results)
}

rebuilt_nuc = fread("/Users/browna6/mahmood/binom/analyzed\ data/nucleotide_mutation_table.csv")
gene_list = unique(rebuilt_nuc[,c("gene","BinnedRole", "RoleinCancer","ensembl_gene_id")])

full_results = get_go_ids(gene_list[,ensembl_gene_id])

setDT(gene_list)[full_results, `Biological Process` := biological_process, on = c(ensembl_gene_id = "ensembl_gene_id")]
setDT(gene_list)[full_results, `Molecular Function` := molecular_function, on = c(ensembl_gene_id = "ensembl_gene_id")]
setDT(gene_list)[full_results, `Cellular Component` := cellular_component, on = c(ensembl_gene_id = "ensembl_gene_id")]


gene_list[,is_cand := 1]
tsg_list = gene_list[BinnedRole == "TSG", c("gene", "is_cand")]
onc_list = gene_list[BinnedRole == "Oncogene", c("gene", "is_cand")]

#doesn't find FAM46C
tsg_list[gene == "FAM46C", gene := "TENT5C"]
tsg_enrich = go_enrich(tsg_list)
tsg_stat = tsg_enrich[[1]]
by(tsg_stat, tsg_stat$ontology, head, n=3)

onc_enrich = go_enrich(onc_list)
onc_stat = onc_enrich[[1]]
by(onc_stat, onc_stat$ontology, head, n=3)

