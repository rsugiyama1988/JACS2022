# Sugiyama et al., JACS
# Process GenBank files
# By Ryosuke Sugiyama, PhD; zy105shsh@gmail.com
# Last update: June 5, 2022

rm(list = ls()) 
library(openxlsx)
library(pipeR)
library(rlist)
library(tidyverse)

library(rentrez)
library(genbankr)
library(seqinr)
library(Biostrings)
#set_entrez_key("xxx") # Activate API key

path1 <- 'xxx/JACS2022/Demo/' #set directory

##### Get nucleotide accessions containing a query rSAM #####

rSAM <- read.xlsx(paste(path1, "Key files/1-08_Dataset S1_source table_demo.xlsx", sep="")) %>>%
  filter(str_detect(rSAM.SSN.cluster, "(01|07)")==FALSE)
Uniprot <- read.xlsx(paste(path1, "Initial datasets/Megacluster-1-1-2_AS45_uniprot_demo.xlsx", sep="")) %>>%
  filter(Length>=300)

rSAM.nuc <- rSAM %>>%
  mutate(Nucleotide.reference = "",
         Nucleotides = "",
         No.of.nucleotides = 0)

pb <- txtProgressBar(min=1, max=dim(rSAM.nuc)[1], style=3)
for (i in 1:dim(rSAM.nuc)[1]) {
  if(rSAM.nuc$No.of.nucleotides[i]==0){
    if(str_detect(rSAM.nuc$Accession[i], "^WP_")){
      temp <- entrez_link(dbfrom = "protein", id = rSAM.nuc$Accession[i], db = "all")
      Sys.sleep(1)
      
      if(is.null(temp$links$protein_nuccore[1])==FALSE) {
        rSAM.nuc$Nucleotide.reference[i] <- temp$links$protein_nuccore[1]
        rSAM.nuc$Nucleotides[i] <- paste(temp$links$protein_nuccore, collapse = "; ")
        rSAM.nuc$No.of.nucleotides[i] <- length(temp$links$protein_nuccore)
      }
    } else {
      temp <- Uniprot %>>%
        filter(Sequence==rSAM.nuc$Sequence[i]) %>>%
        (str_split(.$EMBL, ";")) %>>%
        unlist %>>%
        (.[.!=""])
      sort
      
      if(length(temp)>0){
        rSAM.nuc$Nucleotide.reference[i] <- temp[1]
        rSAM.nuc$Nucleotides[i] <- paste(temp, collapse = "; ")
        rSAM.nuc$No.of.nucleotides[i] <- length(temp)
      }
      }
    }
  setTxtProgressBar(pb, i)
}

rbind(rSAM.nuc[rSAM.nuc$No.of.nucleotides>0,], rSAM.nuc[rSAM.nuc$No.of.nucleotides==0,]) %>>%
write.xlsx(paste(path1, "Intermediate files/GenBank processing/2-01_rSAM for precursor search.xlsx", sep=""))


##### Fetch GenBank files #####

rSAM.nuc <- read.xlsx(paste(path1, "Intermediate files/GenBank processing/2-01_rSAM for precursor search_fixed.xlsx", sep=""))
l <- list.files(paste(path1, "GenBank files/", sep=""))

temp <- rSAM.nuc$Nucleotide.reference %>>%
  (.[.!=""]) %>>%
  unique %>>%
  paste(".gb", sep="")
GB.query <- setdiff(temp,l) %>>%
  str_remove(".gb")

set <- 1:length(GB.query) %>>%
  split(ceiling(seq_along(.)/100))  

for (k in 1:length(set)){
  query <- GB.query[set[[k]]]
  print(paste("set", k, "ongoing"))
  m <- 1
  pb <- txtProgressBar(min=1, max=length(query), style=3)
  for (i in 1:length(query)){
    temp <- entrez_fetch(db="nuccore", query[i], rettype = "gbwithparts")
    write.table(temp, paste(path1, "GenBank files/", query[i], ".gb", sep=""),
                quote = F, row.names = F, col.names = F)
    rm(temp)
    Sys.sleep(1)
    setTxtProgressBar(pb, m)
    m <- m+1
  }
}


##### Extract genetic loci surrounding rSAM #####

f <-  function(filename){readGenBank(paste(path1, "GenBank files/", filename, sep=""), partial = T)}

query.all <- read.xlsx(paste(path1, "Intermediate files/GenBank processing/2-01_rSAM for precursor search_fixed.xlsx", sep=""))
query.all <- query.all %>>%
  mutate(GenBank = paste(Nucleotide.reference,".gb", sep=""))
set <- 1:dim(query.all)[1] %>%
  split(ceiling(seq_along(.)/200))  
load.error <- data.frame()

for (k in 1:length(set)){
  query <- query.all[set[[k]],]
  print(paste("set", k, "ongoing"))
  
  output <- list()
  DNA.sequences <- list()
  t <- 1
  pb <- txtProgressBar(min=0, max=dim(query)[1], style=3)
  
  for (i in 1:dim(query)[1]){
    gb <- try(f(query$GenBank[i]), silent = TRUE)
    
    if (class(gb) == "try-error"){
      load.error <- query[i,] %>>%
        mutate(error="GenBank unloaded") %>>%
        (rbind(load.error, .))
    } else if(is.null(cds(gb)$protein_id)) {
      load.error <- query[i,] %>>%
        mutate(error="No protein ID") %>>%
        (rbind(load.error, .))
    } else {
      strain <- sources(gb) %>>%
        attr("elementMetadata") %>>%
        (if (is.null(.$strain)){
          if (is.null(.$organism)){} else{.$organism}
        } else if (str_detect(.$organism, .$strain)){
          .$organism
        } else {
          paste(.$organism, .$strain)
        })
      Fw <- attr(gb, "sequence")  
      Rv <- attr(gb, "sequence") %>>% complement
      
      cds <- cds(gb)
      seq <- cds$translation %>>% as.vector
      rSAM.location <- paste("^", query$Sequence[i], "$", sep="") %>>% grep(seq)
      
      if(length(rSAM.location)==0)  {
        load.error <- query[i,] %>>%
          mutate(error="no query rSAM") %>>%
          (rbind(load.error, .))
      } else {
        
        dataset <- data.frame(Entry = query$Entry[i],
                              rSAM.SSN.cluster = query$rSAM.SSN.cluster[i],
                              rSAM.query = query$Accession[i],
                              Strain = strain,
                              Source = if(length(gb@accession)==1) {str_split(gb@accession, " ") %>>% unlist %>>% (.[1])} else {""},
                              Locus.tag = if(length(cds$locus_tag)>0) {cds$locus_tag} else {""},
                              Protein.ID = if(length(cds$protein_id)>0) {cds$protein_id} else {""},
                              Note = "",
                              Check = "",
                              Strand = attr(cds, "strand") %>>% (mapply(rep, attr(., "values"), attr(., "length"))) %>>% unlist %>>% as.factor,
                              Product = if(length(cds$product)>0) {cds$product} else {""},
                              Length = str_length(cds$translation),
                              Sequence = seq,
                              Start = attr(cds, "ranges") %>>% attr("start"),
                              End = attr(cds, "ranges") %>>% (attr(.,"start") + attr(.,"width")-1))
        
        region <- lapply(rSAM.location, function(x){max(x-10,1):min(x+10,length(seq))})
        
        for (m in 1:length(region)){
          sub.cds <- dataset[region[[m]],] %>>%
            mutate(Entry = paste(query$Entry[i], formatC(m, width = 2, flag = 0), sep="."))
          
          pseudo <- is.na(sub.cds$Protein.ID)
          sub.cds$Protein.ID[pseudo] <- "pseudo"
          
          sub.cds <- sub.cds %>>% mutate(FASTA = paste(">", Protein.ID,
                                                       " ", Product,
                                                       " [", Strain,
                                                       "]\n", Sequence, sep=""))
          sub.cds$FASTA[pseudo] <- ""
          for (n in 1:dim(sub.cds)[1]) {
            if(sub.cds$Strand[n]=="+"){
              sub.cds$Gene[n] <- str_sub(Fw, sub.cds$Start[n], sub.cds$End[n])
            } else if(sub.cds$Strand[n]=="-") {
              sub.cds$Gene[n] <- str_sub(Rv, sub.cds$Start[n], sub.cds$End[n]) %>>%
                str_split("") %>>% unlist %>>% rev %>>% str_c(collapse="")
            }
          }
          sub.cds$Gene[pseudo] <- ""
          
          L <- 1:dim(sub.cds)[1]
          rSAM.sublocation <- paste("^", query$Sequence[i], "$", sep="") %>>%
            grep(sub.cds$Sequence)
          if(length(rSAM.sublocation)==1){
            cl <- sapply(L, function(x){
              as.character(sub.cds$Strand[x:rSAM.sublocation]) %>>% (length(table(.)))
            })
            sub.cds$Note[cl==1] <- "clustered"
          }
          sub.cds$Note[rSAM.sublocation] <- "rSAM"
          
          
          Nucleotides <- c(sub.cds$Start, sub.cds$End) %>>% (str_sub(Fw, min(.), max(.))) %>>% list
          names(Nucleotides) <- sub.cds$Entry[1]
          
          DNA.sequences <- c(DNA.sequences, Nucleotides)
          
          sub.cds <- list(sub.cds)
          names(sub.cds) <- sub.cds$Entry[1]
          output <- c(output, sub.cds)
        }
      }
    }
    setTxtProgressBar(pb, t)
    t <- t+1
  }
  
  output.all <- list.rbind(output)
  
  save(output, file=paste(path1, "Intermediate files/GenBank processing/Partial data/Genetic loci_set ", k, ".Rdata", sep=""))
  save(DNA.sequences, file=paste(path1, "Intermediate files/GenBank processing/Partial data/DNA seq_set ", k, ".Rdata", sep=""))
  write.xlsx(output.all, paste(path1, "Intermediate files/GenBank processing/Partial data/Genetic loci_set ", k, ".xlsx", sep=""))
}

query.all$Obsolete <- ""
if(dim(load.error)[1]>0){
for (i in 1:dim(load.error)[1]){
  query.all$Obsolete[grep(load.error$Entry[i], query.all$Entry)] <- paste("load error", load.error$type[i])
}
}

query.all %>>% arrange(Obsolete, Entry) %>>%
  write.xlsx(paste(path1, "Intermediate files/GenBank processing/2-02_rSAM for precursor search_loaded.xlsx", sep=""))


files.xlsx <- list.files(paste(path1, "Intermediate files/GenBank processing/Partial data/", sep="")) %>>%
  (.[str_detect(.,"set.+xlsx")])
G.loci <- read.xlsx(paste(path1, "Intermediate files/GenBank processing/Partial data/", files.xlsx[1], sep=""))
if(length(files.xlsx)>1){
for (i in 2:length(files.xlsx)){
  G.loci <- read.xlsx(paste(path1, "Intermediate files/GenBank processing/Partial data/", files.xlsx[i], sep="")) %>>%
    (rbind(G.loci, .))
}
}
G.loci <- G.loci %>>% arrange(Entry, Start)
write.xlsx(G.loci, paste(path1, "Intermediate files/GenBank processing/2-03_Combined genetic loci_demo.xlsx", sep=""))




