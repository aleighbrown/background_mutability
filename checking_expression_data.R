rebuilt_amino <- fread("~/mahmood/binom/analyzed data/codon_substitution_table.csv")
rebuilt_amino[BScore <=5.462115e-05,Status85 := "Driver"]

#quickly rename some genes. never using Hugo symbols again after this project. Too many alias
rebuilt_amino[gene == "NSD3", gene := "WHSC1L1"]
rebuilt_amino[gene == "KNL1", gene := "CASC5"]
rebuilt_amino[gene == "AFDN", gene := "MLLT4"]

rebuilt_amino[gene == "WDCP", gene := "C2orf44"]
rebuilt_amino[gene == "NSD2", gene := "WHSC1"]

driver_genes = rebuilt_amino[Status85 == "Driver",unique(gene)]

 
ccle_data = fread("/Users/browna6/mahmood/binom/analyzed data/CCLE_RNAseq_genes_rpkm_20180929.gct.gz")

#melt the data
ccle_data_RK_melted = melt(ccle_data,measure.vars = names(ccle_data[,3:1021]), id.vars = c("Name","Description"))
setnames(ccle_data_RK_melted, c("Description","variable","value"),c("Gene","ccle_name","RKPM"))

#split up the tissue name into name and tissue
ccle_data_RK_melted[,tissue := str_split(str = ccle_name, pattern = "_", n = 2,simplify = T)[,2]]


tissue_gene_expression = ccle_data_RK_melted[,mean(RKPM), by = c("Gene","tissue")]



notthere = setdiff(rebuilt_amino[,unique(gene)],tissue_gene_expression[,unique(Gene)])

all_expression = tissue_gene_expression[Gene %in% rebuilt_amino[,unique(gene)]]
all_expression = all_expression[tissue != "" & tissue != "NS"]
driver_expression = tissue_gene_expression[Gene %in% driver_genes]

all_expression[V1 > 0.5 ]
driver_expression = tissue_gene_expression[Gene %in% driver_genes]
driver_expression[V1 > 0.5 ]
