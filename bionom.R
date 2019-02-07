library(jsonlite)
library(data.table)
library(Biostrings)
library(stringr)
library(stringdist)
library(EnsDb.Hsapiens.v86)
library(rtracklayer)

#import the bed file of silent mutations
#MUTA_TRANSCRIPTS <- fread("/net/pan1/mutagene/analysis/Mutability/genes_transcripts.txt")

BED <- import("/net/pan1/mutagene/analysis/Mutability/mutagene_mutations.bed.gz",format = "BED")
BED#flip the start of the bed file
start(BED) <- end(BED)
#get the Hsapiens up
ENB = EnsDb.Hsapiens.v86

TX <- transcripts(ENB, columns=c("tx_id", "gene_id", "gene_name","tx_biotype","protein_sequence","protein_id"))
TX = keepStandardChromosomes(TX, pruning.mode = "tidy")
MAPPING <- as.data.table(cbind(tx_id=TX$tx_id, name=TX$gene_name, gene_id=TX$gene_id,sequence=TX$protein_sequence))

#get the lengths of all chromo
CHR.LENGTHS <- seqlengths(keepStandardChromosomes(genes(ENB), pruning.mode = "coarse"))


#read in all the SNPS and then make it into a GRange 
t <- fread("~/mahmood/snp_locals.gff",header = F)
SNPTable <- makeGRangesFromDataFrame(t,seqnames.field = "V1", start.field = "V2", end.field = "V2")
rm(t)

#rstudio runs the wrong verison of python on my local, so this just to make sure it gets the right one
Sys.setenv(PATH = paste("/Users/browna6/conda/bin",Sys.getenv("PATH"),sep=":"))

#for pancancer for the number of mutations total samples was 9450
TOTALSAMP = 12013

#create a temporary table to switch between the codons and the amino acids, this GENETIC_CODE comes from the biostrings package
GENETABLE <- data.table(cbind(names(GENETIC_CODE),GENETIC_CODE))
names(GENETABLE) <- c("codon","amino")


#function takes an accession and returns a count of the number of mutations per nucleotide for that model, default is the pancancer model
get.mutationPerNucleotide <- function(acc = 33999){
  #read the counts from the mutagene website for the given signature
  temp <- fread(paste0("https://dev.ncbi.nlm.nih.gov/research/mutagene/api/signature_data/count/",acc))
  #add a column for wildtype and mutant type
  temp[,wt := substr(V1,3,3)]
  temp[,mutant := substr(V1,5,5)]

  #create a table which sums the  number of mutations by the wildtype and the mutant
  mutationPerNucleotide <- temp[,sum(V2), by = c("wt", "mutant")]
  #since it's easiest to do things with table witht the full, and a T > C is equivalent to an A > G, we're going to double the table, then take the reverse complement of the doubled wt and mutatant
  mutationPerNucleotide <- rbind(mutationPerNucleotide,mutationPerNucleotide)
  mutationPerNucleotide[7:12,wt := as.character(reverseComplement(DNAStringSet(wt)))]
  mutationPerNucleotide[7:12,mutant := as.character(reverseComplement(DNAStringSet(mutant)))]
  mutationPerNucleotide[,prob := V1 / NUC_FREQ_EXOME[wt]]
  return(mutationPerNucleotide)
}



#function takes the gene and returns the sequence, uses mutagene api
get.Seq <- function(gene){
  #uses the mutagene API to download the sequence of the gene in question
  sequenceAPI <- "https://dev.ncbi.nlm.nih.gov/research/mutagene/api/sequence/"
  #upper case just in case 
  gene <- toupper(gene)
  seq <- fromJSON(paste0(sequenceAPI,gene))
  return(seq)
}

###take a string, position in the string, and a character, and substitutes that character at that string
subchar2 <- function(string, pos, char) { 
    substr(string, pos, pos) <- char
  return(string) 
}

subchar <- function(string, pos, char) { 
  for(i in pos) { 
    string <- gsub(paste("^(.{", i-1, "}).", sep=""), "\\1", string) 
  } 
  string 
}


##adds the condon that the nucleotide substition belongs to in the protein
bind.codon <- function(geneMutability, transcript){
  #make a data table of all the codons in the protein 
  c <-  as.data.table(codons(DNAString(transcript$sequence)))
  #  #give c a column with it's row number, this is the amino acid position, we'll use it later
  c <- c[,ind := .I]
  #repeat each codon three times, this will make it easier for me to index in a bind
  c <-  c[rep(seq_len(nrow(c)), each=3),]
  #give c a column with it's row number, used for matching later
  c <- c[,ind2 := .I]
  #now use the pos value of the DNA to match to it's respective codon
  setDT(geneMutability)[c, codon := i.x, on = c(pos = "ind2")]
  setDT(geneMutability)[c, aminoposition := i.ind, on = c(pos = "ind2")]

  
  return(geneMutability)
}


#function takes the gene and gets a table of its mutability and observed cancer mutations per nucleotide site. also provides a table with additional information used for later calculations default accession number is for pancancer, but others can be used
get.geneNucMutability <- function(gene, acc = 33999){
  #upper case just in case
  gene <- toupper(gene)

  #uses a pipe call to python function to read mutagene API, and returns this as a data.table
  #additionally, there's some things I kind of want to steal and smoosh together from two separate downloads so I'm calling it twice this is probably a bad lazy slow way to do this
  call1 <- paste("python ~/mahmood/binom/call_mutability_api.py", gene, as.character(acc),"5")
  tryCatch(geneMutability <- fread(call1),error=function(e) return(0))

  #good thing I have three different ways to test breakage....only one of which fails properly i'm sure...
  if ( !exists("geneMutability")){
    print("This gene not found.")
    print(gene)
    return(0)
  }else if(nrow(geneMutability) == 0){
    print("This gene not found.")
    print(gene)
    return(0)
  }
  

  #basically for calculation of probability for each nucleotide I want the trinuc sequence, there's probably a smart way to do this but I'm going to take the numeric
  #value out of the mutiation in the first data table using regex, add a column called position to that table, set it as the key, and then pull the trincule from stealme
  #and throw that sucka away

  geneMutability[,pos:= as.integer(str_extract(mutation,'[[:digit:]]+'))]
  
  #adding in an observed all col
  geneMutability[,observed_all := observed_wgs + observed_other]

  
  
  #I want to have easy access to the wildtype and mutant DNA base
  geneMutability[, wildtype := str_extract_all(geneMutability[,mutation],"[[:alpha:]]|[[:punct:]]",simplify = T)[,1]]
  
  geneMutability[, mutant := str_extract_all(geneMutability[,mutation],"[[:alpha:]]|[[:punct:]]",simplify = T)[,2]]
  
  #currently for scaling purposes mutagene puts a 10^6 in numerator, but I don't want or need that here, so pulling that out
  geneMutability[,mutability := mutability/10^6]

  return(geneMutability)
}
#takes the table of the mutability of nucleotide sites  + observed mutations and adds a column with the mutant codon -> mutant codon used later for calculating the probability of an amino acid substitution
mutate.codons <- function(geneMutability, transcript){
  
  #add the codons
  bind.codon(geneMutability, transcript)
  
  #add an index of the position of th nucleotide in the codon
  geneMutability[,codonposition := pos %% 3]
  geneMutability[codonposition == 0, codonposition:= 3]
  #return the mutant codon
  geneMutability[,mutcodon := subchar2(codon,codonposition,mutant), by = 1:nrow(geneMutability)]
  #translate the mutant codon
  setDT(geneMutability)[GENETABLE, mutAA := amino, on = c(mutcodon = "codon")]
  #translate the wt codon so I have that as well
  setDT(geneMutability)[GENETABLE, wildtypeAA := amino, on = c(codon = "codon")]
  return(geneMutability)
}



#this function returns a table of the mutability value for every amino acid subsitution possible through a SNP along a protein, similar to the nucleotide mutability, default is pancancer as 
#mutability, but others can be used
get.geneAminoMutability <- function(gene, acc = 33999){
  #this pulls the mutability by amino acid in the peptide
  #gene <- toupper(gene)
  #this uses the api to download the mutability
  call3 <- paste("python ~/mahmood/binom/call_mutability_api.py", gene, as.character(acc),"3")

  tryCatch(aminomuta <- fread(call3),error=function(e) return(0))

  #good thing I have three different ways to test breakage....only one of which fails properly i'm sure...
  if (!exists("aminomuta")){
    print("This gene not found.")
    print(gene)
    return(0)
  }else if(nrow(aminomuta) == 0){
    print("This gene not found.")
    print(gene)
    return(0)
  }
  
  
  #pull position information, for the amino acid
  aminomuta[,pos:= as.integer(str_extract(mutation,'[[:digit:]]+'))]
  #this make my life easier later on with taking wt and mutat aa so they're easier to futz with
  aminomuta[, wildtype := str_extract_all(aminomuta[,mutation],"[[:alpha:]]|[[:punct:]]",simplify = T)[,1]]
  aminomuta[, mutant := str_extract_all(aminomuta[,mutation],"[[:alpha:]]|[[:punct:]]",simplify = T)[,2]]
  

  #I don't want this scaled probability because that's not helpful for my purposes, so pull out that 10^6 factor
  aminomuta[,Gene := gene]
  return(aminomuta[,mutability:=mutability/10^6])
  
}

#reads in the frequency of trinucleotides in the human exome and writes it to a data.table, 
read.exomeFreq <- function(){
  exomefreq <- fread('~/mahmood/binom/exome_trinucleotides.txt')
  colnames(exomefreq) <-  c("trinucleotide","freq")
  #I'm appending the reverseve complements to the end with the same frequency because later I will want to run down this table
  #and take these frequenices to shove into geneMutatbility and running once down a 64 row table will be easiest
  exomefreq <- rbindlist(list(exomefreq, exomefreq[,list(as.character(reverseComplement(DNAStringSet(trinucleotide))),freq)]))
  return(exomefreq)
}


#calculates a weight of gene mutation by taking # unique  mutations / # possible sites
calc.geneWeight <- function(dtab, gene){
  #gene weight is the number of mutated nucleotide sites / number of nucleotides in the coding transcript, 
  geneMutability = copy(dtab[Gene == gene])
  geneWeight <- (geneMutability[CountNew != 0,length(unique(pos))]) / (geneMutability[,max(pos)])

  return(geneWeight)
}


#calulates the gene weight as the number of SNPS that happen within a window
calc.geneWeightSNPS <- function(gene, window = 50000){
    #
    #find the largest range that includes the gene, to do this call the GRange reduce function on the TX object
    gene_range <- GenomicRanges::reduce(TX[TX$gene_name == gene & TX$tx_biotype != "LRG_gene"])
    tot <- gene_range + window
    #make sure that if total range goes off the end, we're going to trim to one 
    tot <- trim(tot)
    #check that the window doesn't fall off the end of the gene, if the end of the expanded window is greater than the length of the chromosome it's on
    if(end(tot) > CHR.LENGTHS[as.character(seqnames(gene_range)@values)]){
      end(tot) <- CHR.LENGTHS[as.character(seqnames(gene_range)@values)]
      #if the start falls off it's just teh whole chromosome
      print("Window extends length of chromosome")
    }
    
    nSNPS <- countOverlaps(tot, SNPTable)
    denom <- width(tot)
    ov <- subsetByOverlaps(TX,tot)

    return(nSNPS / denom)
}

#calulates the gene weight as the number of synon mutations in cancer
calc.geneWeightSilent <- function(gene, window = 100000){
  #
  
  #find the isoform transcript from mutagene, it looks like mutagene takes the first isoform
  tscript <- MAPPING[name==gene,min(tx_id)]
  nsynmutations <- as.numeric(suppressWarnings(countOverlaps(TX[TX$tx_id == tscript] + window,BED)))
  return(nsynmutations / window)
}


#takes a transcript and returns a table with the trinuc frequency
calc.triNucTable <- function(transcript){
  #takes the transcript of a gene and uses a funciton out of Biostrings to compute 
  #table of frequencies, converts it to a data.table cause those as the bomb
  triNucTable <-  as.data.table(trinucleotideFrequency(DNAString(transcript$sequence)),keep.rownames = T)
  return(triNucTable)
}

#calculates the probability of observing a given Nucleotide substitution
calc.ProbMutinGene <- function(geneMutability,geneWeight, triNucTable){
  #add a column to geneMutability which contains the number of the trinuc in the gene
  setDT(geneMutability)[triNucTable, freqInGene := i.V2, on = c(trinuc = "V1")]

  #the probability in our data set to see a particular mutation in a gene is the probablity of that mutation to occur in general, reduced by
  #the geneweight if you factor in that any gene might be itself under selective pressure
  #geneMutability[,phtg := (pht * geneWeight)]
  geneMutability[,phtg := (pht * geneWeight)/freqInGene]

  return(geneMutability)
}

  
#calculates the probability of observing a given amino acid substitution given the mutability of the DNA, uses either the frequentist probabiliy of a nucleotide sub or a conditional prob or takes the column which is the neighborhood weight * probability
calc.ProbAminoAcidSubinGene <- function(protpos, geneMutability, mutantAmino, prob = "c"){
  if(prob == "f"){
    col <-  'phtg'
  }else if(prob == "c"){
    col <- 'condphtg'
  }else{
    col <- 'pGTH'
  }
  
  total <- (geneMutability[aminoposition == protpos & mutAA == mutantAmino, get(col)])
  
  probP <- 1 - prod(1 - total)
  probS <-sum(total)
  return(probP)
}


#takes a list of genes which neighbor a gene and calculates their mutation burden as the numer of silent mutations over the length
calc.neighborWeights <- function(neighbors){
  totalWeight <- 0
  for(n in neighbors){
    ts <- get.Seq(n)
    
    geneMut <- get.geneNucMutability(n,transcript = ts)
    #if the neighbor gene wasn't found, go to the next loop
    if (geneMut == 0){
      next
    }
      
    silentWeight <- calc.geneWeightSilent(geneMut)
    totalWeight <- totalWeight + silentWeight
    rm(ts,geneMut,silentWeight)
    
  }
  return(totalWeight)
}




#find all the transcripts in MutaGene for a gene or genes
find.transcripts <- function(gene){
  all_transcripts <- data.table()
  for(i in 1:length(gene)){
    all_transcripts <- rbind(all_transcripts,cbind(gene[i],MUTA_TRANSCRIPTS[like(V1,gene[i])]))
  }
  setnames(all_transcripts,c("Gene","Found","One","two"))
  #filter out the false positives
  
  all_transcripts <- all_transcripts[Gene == Found | !(Gene != str_extract(Found,".+?(?=_)"))]
  
  return(all_transcripts)
}

#given a list of Genes, return a data table of the count
get.fullAminoMutability <- function(all_genes, acc = 33999){
  #make an empty data table to be filled in a loop
  not_found <- character()
  allamino <- data.table()
  for (t in all_genes){
    print(t)
    #should write somethin to catch errors here...
    
    aminoMut <- get.geneAminoMutability(t, acc = acc)
    #didn't find the gene?
    if(aminoMut == 0){
      not_found <- c(t,not_found)

      next()
    }
    aminoMut[,Gene := t]
    
    allamino<- rbind(allamino,aminoMut)
    rm(aminoMut)
    
  }
  

  
  return(list(allamino,not_found))
}
#given a list of Genes, return a data table of the count
get.fullNucMutability <- function(all_genes, acc = 33999){
  #make an empty data table to be filled in a loop
  not_found <- character()
  allnucleo <- data.table()
  for (t in all_genes){
    print(t)
    
    nucMut <- get.geneNucMutability(t,acc = acc)
    #didn't find the gene?
    if(nucMut == 0){
      not_found <- c(t,not_found)
      next()
    }
    nucMut[,Gene := t]
    
    allnucleo<- rbind(allnucleo,nucMut)
    rm(nucMut)
    
  }
  
  
  
  return(list(allnucleo,not_found))
}

readrow <- function()
{ 
  n <- readline(prompt="Which row do you want? ")
  return(as.integer(n))
}

get_peptide_id <- function(gene){
  if(!exists("ensembl")){
    library(biomaRt)
    ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  }
  #use the mutagene transcript file to find the transcript id
  print(MUTA_TRANSCRIPTS[like(V1,gene)])
  rn <- readrow()
  ensembl_trans <- MUTA_TRANSCRIPTS[like(V1,gene)][rn,2]
  ans <- getBM(attributes=c('ensembl_transcript_id','hgnc_symbol',"ensembl_peptide_id"),filters = 'ensembl_transcript_id', values = ensembl_trans, mart = ensembl)
  print(ans)
  return(ans$ensembl_peptide_id)
}

#we don't have an api that returns the mutability values straight for a given model, so I'm going to write a work around
get_sig_table <- function(ACC = 34275){
  table_folder = "~/mahmood/binom/sig_tables/"
  current_tables <- list.files(table_folder)
  if(paste0(ACC,".csv") %in% current_tables){
    sig_table2 <- fread(paste0(table_folder,ACC,".csv"))
    return(sig_table2)
  } else{
  #to get the mutability for any given accession I query for brca1, then take the mean by context, wildtype, mutant, size BRCA1 has all the context, we get there
  br_nuc = get.geneNucMutability("BRCA1",acc = ACC)
  temp = br_nuc[,mean(mutability), by = c("context","wildtype","mutant")]
  setnames(temp, "V1", "mutability")
  fwrite(temp[mutability != 0],paste0(table_folder,as.character(ACC),".csv"))
  return(temp[mutability != 0])
  }
}

###this doesn't so much get the principal transcript as explore it
appr_anno= fread("http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/appris_data.principal.txt",header = F)
appr_anno[grep("PRINCIPAL",V5), Princ := as.integer(str_extract(V5,'[[:digit:]]+'))]

get_principal_transcript <- function(gene, ident = "hugo_symbol"){
  #read in the appris annotations on principal transcripts
  if(!exists("appr_anno")){
    appr_anno= fread("http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/appris_data.principal.txt",header = F)
    appr_anno[grep("PRINCIPAL",V5), Princ := as.integer(str_extract(V5,'[[:digit:]]+'))]
  }
  
  #if the identifer is hugo symbol, check to see if you can find it
  if(ident == "hugo_symbol"){
    if(!gene %in% appr_anno[,V1]){
      return("Hugo Symbol not found, try an alias or searching Ensembl gene")
    } else{
      print(appr_anno[V1 == gene])
    }
  } else if(ident == "ens"){
    if(!gene %in% appr_anno[,V2]){
      return("Ensembl gene not found")
    } else{
      print(appr_anno[V2 == gene])
    }
  }
}

append_mutability <- function(mt, acc = 33999, name = "Pan"){
  mut_table = get_sig_table(acc)
  n = paste0("mutability",name)
  setDT(mt)[mut_table, (n) := i.mutability, on = c(context = "context", wildtype = "wildtype", mutant = "mutant")]
  return(mt)
}

