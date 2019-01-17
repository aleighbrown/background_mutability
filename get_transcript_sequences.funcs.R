
#given a list of CCDS  returns the exon ids, coding information, and genomic locations
get_exonIds <- function(transcript_list){
  if(!exists("mart")){
    require(biomaRt)
    mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  }
  exonIds = as.data.table(getBM(mart, attributes = c("external_gene_name","ensembl_exon_id","ensembl_gene_id","ensembl_transcript_id","cdna_coding_start","cdna_coding_end","exon_chrom_start","exon_chrom_end","genomic_coding_start","genomic_coding_end","chromosome_name","rank","strand","is_constitutive","cds_start","cds_end"), 
                                filters = "ccds", values = transcript_list))
  setorder(exonIds, rank, ensembl_transcript_id)
  return(exonIds)
}

#given a list of exonIds, uses the getSeq function to find the sequences, and the 1 bp flanks
get_nucSequence <- function(exonIds){
  #load the HSapiens object
  if(!exists("Hsapiens")){
    require(BSgenome.Hsapiens.NCBI.GRCh38)
  }
  #we don't need the exons which don't contain any coding DNA
  tempCDNA = exonIds[!is.na(cdna_coding_start)]
  setorder(tempCDNA, rank, ensembl_transcript_id)
  #normalize strand information
  tempCDNA[strand == -1, Strand := "-"]
  tempCDNA[strand == 1, Strand := "+"]
  tempCDNA[,strand := NULL]
  #get the sequence for all the exons with coding dna
  tempCDNA[,coding_seq := as.character(getSeq(Hsapiens, as.character(chromosome_name), start = genomic_coding_start, end = genomic_coding_end, strand = Strand))]
  #get the nucleotide right before the coding start
  tempCDNA[, UpStream := as.character(getSeq(Hsapiens, as.character(chromosome_name), start = genomic_coding_start - 1, end = genomic_coding_start - 1, strand = Strand))]
  #get the nucleotide right after the coding end
  tempCDNA[, DownStream := as.character(getSeq(Hsapiens, as.character(chromosome_name), start = genomic_coding_end + 1 , end = genomic_coding_end + 1, strand = Strand))]
  
  return(tempCDNA)
}

#given an exon's coding sequence, and the upstream and downstream nuc, breaks it up into contexts
sliding_context <- function(sequence){
  context_table = data.table()
  for(i in 1:(nchar(sequence) - 2)){
    context = substr(sequence,i,(i + 2))
    context_table = rbind(context_table,cbind(i,context))
  }
  
  context_table[,i := as.numeric(i) + 1]
  setnames(context_table, "i", "pos")
  return(context_table)
}
#given an exon's cdna context table, and the up stream and downstream, add those positions
add_up_down <- function(context_table, up, down){

  upstream = paste0(up,context_table[1,substr(context,1,2)])
  downstream = paste0(context_table[nrow(context_table),substr(context,2,3)],down)
  context_table = rbind(context_table,cbind(pos = 1, context = upstream, exon_boundary = "Y"), cbind(pos = (nrow(context_table) + 2), context = downstream, exon_boundary = "Y"), fill = T)
  context_table[,pos := as.numeric(pos)]
  setorder(context_table, pos)
  return(context_table)

}

#takes a coding table and widens it for the wildtype and possible mutants
mutate_coding_table <- function(coding_table){
  coding_table[,wildtype := substr(context, 2,2)]
  coding_table[wildtype == "A", `:=` (m1 = "C", m2 = "T", m3 = "G")]
  coding_table[wildtype == "C", `:=` (m1 = "A", m2 = "T", m3 = "G")]
  coding_table[wildtype == "T", `:=` (m1 = "C", m2 = "A", m3 = "G")]
  coding_table[wildtype == "G", `:=` (m1 = "C", m2 = "T", m3 = "A")]
  ids = setdiff(names(coding_table),c("m1","m2","m3"))
  temp = melt(coding_table, id.vars = ids, measure.vars = c("m1","m2","m3"))
  temp[,variable := NULL]
  setnames(temp, "value","mutant")
  setorder(temp, pos)
  return(temp)
}

#given a data table which contains all the exons for a single transcript, returns a data table which contains the exons glued together with their context
create_coding_table <- function(dt){
  coding = data.table()
  #for first exon rank with cdna to the max exon rank
  for(i in min(dt[,rank]):max(dt[,rank])){
    #split the exon into chunks 
    cxt_table = sliding_context(dt[rank == i, coding_seq])
    #add the start and end
    cxt_table = add_up_down(cxt_table, dt[rank == i, UpStream], dt[rank == i, DownStream])
    #include chromosome name
    cxt_table[,chromosome_name := dt[rank == i, chromosome_name]]
    #append the genomic coding, if a forward strand we start at the start
    if(dt[rank == i, Strand] == "+"){
      cxt_table[,genomic_coding := dt[rank == i,genomic_coding_start] : dt[rank == i, genomic_coding_end]]
      #if the reverse strand the coordin starts at the end
    }else if((dt[rank == i, Strand] == "-")){
      cxt_table[,genomic_coding := dt[rank == i,genomic_coding_end] : dt[rank == i, genomic_coding_start]]
    }
    #append which exon rank we're in 
    cxt_table[,ensembl_transcript_id := dt[rank == i,ensembl_transcript_id]]
    if(i != min(dt[,rank])){
      cxt_table[,pos := pos + (dt[rank == (i - 1), cds_end])]
    }
    cxt_table[,exon_rank := i]
    coding = rbind(coding, cxt_table, fill = T)
  }
    coding = mutate_coding_table(coding)
  
  return(coding)
}

#takes a table with transcripts, exon info, seq info, and makes a large table of all the genes
create_all_transcript_coding_tables <- function(transcript_table){
  tx_list = transcript_table[,unique(ensembl_transcript_id)]

  full_table = data.table()
  for(t in tx_list){
    ct = create_coding_table(transcript_table[ensembl_transcript_id == t])
    full_table = rbind(full_table,ct)

  }
  return(full_table)
}
#takes a vector of nucleotide positions arranged in mutation table format,where each nucleotide position is repeated 3 times, and returns of vector of the codon position
return_codon_position_vector <- function(pos){
  k = max(pos) / 3 
  vec = c()
  for(j in 1:k){
    vec = c(vec,rep(j, 9))
  }
  return(vec)
}
#takes the exons with CDNA table and returns the codons and there position for all the transcripts in that table
create_full_cdna_table <- function(tempCDNA){
  if(!exists("GENETABLE")){
    require(Biostrings)
    #create a temporary table to switch between the codons and the amino acids, this GENETIC_CODE comes from the biostrings package
    GENETABLE <- data.table(cbind(names(GENETIC_CODE),GENETIC_CODE))
    names(GENETABLE) <- c("codon","amino")
  }
  test = tempCDNA[,  list(paste(coding_seq, collapse="")), by = ensembl_transcript_id][,as.data.table(codons(DNAString(V1))), by = ensembl_transcript_id]
  test[,codon_pos := sequence(rle(as.character(test$ensembl_transcript_id))$lengths)]
  setDT(test)[GENETABLE, amino := i.amino, on = c(x = "codon")]
  return(test)
}

###take a string, position in the string, and a character, and substitutes that character at that string
subchar2 <- function(string, pos, char) { 
  substr(string, pos, pos) <- char
  return(string) 
}

mutate.codons = function(geneMutability){
  if(!exists("GENETABLE")){
    require(Biostrings)
    #create a temporary table to switch between the codons and the amino acids, this GENETIC_CODE comes from the biostrings package
    GENETABLE <- data.table(cbind(names(GENETIC_CODE),GENETIC_CODE))
    names(GENETABLE) <- c("codon","amino")
  }
  #add an index of the position of th nucleotide in the codon
  geneMutability[,codonposition := pos %% 3]
  geneMutability[codonposition == 0, codonposition:= 3]
  #return the mutant codon
  geneMutability[,mutcodon := subchar2(codon,codonposition,mutant), by = 1:nrow(geneMutability)]
  #translate the mutant codon
  setDT(geneMutability)[GENETABLE, mutAA := amino, on = c(mutcodon = "codon")]

  return(geneMutability)
}
