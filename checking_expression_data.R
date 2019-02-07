rebuilt_amino = fread("/Users/browna6/mahmood/binom/analyzed\ data/codon_substitution_table.csv")
rebuilt_amino[,Status85:=(ifelse(BScore <=5.462115e-05,"Driver", ifelse(BScore<2.945374e-05,"Potential Driver","Undefined")))]
driver_genes = allamino_census[Status85 == "Driver",unique(Gene)]

 
ccle_data = fread("https://data.broadinstitute.org/ccle/CCLE_DepMap_18Q2_RNAseq_RPKM_20180502.gct")
 
#ccle_data_readcount = fread("https://data.broadinstitute.org/ccle/CCLE_DepMap_18Q2_RNAseq_reads_20180502.gct")
#this is the only gene which cannot be found in the RNA seq download, but on the website the gene can be found, so I downloaded straight from there
ccle_IRS4 = fread("/Users/browna6/mahmood/binom/analyzed\ data/IRS4.txt")

ccle_IRS4_melted = melt(ccle_IRS4,measure.vars = names(ccle_IRS4[,2:1458]), id.vars = "Gene")
setnames(ccle_IRS4_melted,c("variable","value"),c("ccle_name","RKPM"))
ccle_IRS4_melted[,tissue := str_split(str = ccle_name, pattern = "_", n = 2,simplify = T)[,2]]
ccle_IRS4_melted[is.na(RKPM),RKPM := 0]

#melt the data
ccle_data_RK_melted2 = melt(ccle_data,measure.vars = names(ccle_data[,3:1078]), id.vars = c("Name","Description"))
setnames(ccle_data_RK_melted2, c("Description","variable","value"),c("Gene","ccle_name","RKPM"))

#split up the tissue name into name and tissue
ccle_data_RK_melted2[,tissue := str_split(str = ccle_name, pattern = "_", n = 2,simplify = T)[,2]]
ccle_data_RK_melted = rbind(ccle_data_RK_melted,ccle_IRS4_melted,fill = TRUE)


tissue_gene_expression = ccle_data_RK_melted[,mean(RKPM), by = c("Gene","tissue")]


notthere = setdiff(allamino_census[Status85 == "Driver",unique(Gene)],ccle_data_RK_melted[,unique(Gene)])

driver_expression = tissue_gene_expression[Gene %in% driver_genes]