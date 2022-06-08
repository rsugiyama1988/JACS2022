# Sugiyama et al., JACS
# Setup of rSAM queries for SSN
# By Ryosuke Sugiyama, PhD; zy105shsh@gmail.com
# Last update: June 8, 2022

rm(list = ls()) 
library(openxlsx)
library(pipeR)
library(rlist)
library(tidyverse)

library(seqinr)
library(taxize)
library(rentrez)
#set_entrez_key("xxx") # Activate API key

path1 <- 'xxx/JACS2022/Demo/' #set directory
path1 <- 'D:/OneDrive/Aica paper/JACS2022/Demo/' #set directory


##### Processing of UniProt accessions #####

Uniprot <- read.xlsx(paste(path1, "Initial datasets/Megacluster-1-1-2_AS45_uniprot_demo.xlsx", sep="")) %>>%
  filter(Length>=300)
Uniprot[is.na(Uniprot)] <- ""

# This table can be obtained 'Retrieve/ID mapping' function in Uniprot.org by
# submitting protein accessions listed in 'Megacluster-1-1-2_AS45_demo.txt'.

Uniprot.RefSeq <- str_split(Uniprot$RefSeq, ";") %>>%
  lapply(function(x){str_extract(x, "^WP_.+") %>>% (.[1])})
for(i in 1:length(Uniprot.RefSeq)){
  Uniprot$RefSeq[i] <- Uniprot.RefSeq[[i]]
}

Uniprot$ID <- Uniprot$Entry
Uniprot$ID[grep(FALSE, is.na(Uniprot$RefSeq))] <- Uniprot$RefSeq[grep(FALSE, is.na(Uniprot$RefSeq))]
Uniprot.2 <- Uniprot %>>%
  summarise(ID,
            Description=paste(">", ID, " ", Protein.name, " [", Organism, "]", sep=""),
            Sequence,
            Dataset="Uniprot",
            Protein.name,
            Organism)
Uniprot.2 <- str_detect(Uniprot.2$ID, "WP_") %>>%
  (Uniprot.2[c(grep(TRUE, .), grep(FALSE, .)),])

##### Processing of RefSeq accessions #####

RefSeq.queries <- paste(path1, "Initial datasets/", sep="") %>>%
  list.files %>>%
  (.[str_detect(.,"^query_blast")])
dataset <- c('Fxs', 'Grr', 'Xye', 'Dar')
RefSeq.all <- matrix(ncol=4)[-1,] %>>% data.frame
colnames(RefSeq.all) <- c("ID", "Description", 'Sequence', 'Dataset')

for (i in 1:length(RefSeq.queries)){
  fasta.q <- paste(path1, "Initial datasets/", RefSeq.queries[i], sep="") %>>%
    read.fasta(seqtype = "AA", as.string = T)
  
  l <- lapply(fasta.q, str_length) %>>% unlist
  fasta.q <- fasta.q[l>=300]
  
  
  RefSeq.all <- data.frame(ID = names(fasta.q),
                       Description = map(fasta.q, ~ {. %>>% attr("Annot")}) %>>% unlist,
                       Sequence = unlist(fasta.q),
                       Dataset=dataset[i]) %>>% (rbind(RefSeq.all,.))  
}

RefSeq.all <- RefSeq.all %>>%
  mutate(Protein.name = str_remove_all(Description, "(^>[^ ]+ | \\[.+\\]|MULTISPECIES: )"),
         Organism = str_extract(Description, " \\[.+\\]") %>>%
           str_sub(3,-2) %>>%
           str_replace_all("(aracter\\(0|\\[|\\])", ""))

rbind(RefSeq.all, Uniprot.2) %>>%
  write.xlsx(paste(path1, "Intermediate files/Data setup/1-01_all rSAMs with duplicates.xlsx", sep=""))


##### Generate a combined list of rSAMs #####

rSAM.all <- read.xlsx(paste(path1, "Intermediate files/Data setup/1-01_all rSAMs with duplicates.xlsx", sep=""))
unique.seq <- table(rSAM.all$Sequence) %>>% sort(decreasing = T) %>>% names

rSAM <- lapply(unique.seq, function(x){
  rSAM.all %>>%
    filter(Sequence==x) %>>%
    summarise(ID=ID[1],
              Description=Description[1],
              Sequence=Sequence[1],
              Dataset=unique(Dataset) %>>% paste(collapse = ";"),
              Protein.name=Protein.name[1],
              Organism=Organism[1])
}) %>>% list.rbind
write.xlsx(rSAM, paste(path1, "Intermediate files/Data setup/1-02_all rSAMs.xlsx", sep=""))

rSAM %>>%
  filter(str_detect(Dataset, "Uniprot")==FALSE) %>>%
  (write.table(.$ID, paste(path1, "Intermediate files/Data setup/1-03_RefSeq for domain retrieving.txt", sep=""),
            quote = F, row.names = F, col.names = F))


##### Annotate domain IDs #####

rSAM <- read.xlsx(paste(path1, "Intermediate files/Data setup/1-02_all rSAMs.xlsx", sep=""))

domain <- read.xlsx(paste(path1, "Initial datasets/Domain list of RefSeq proteins_demo.xlsx", sep=""))
domain[is.na(domain)] <- ""

domain <- Uniprot %>>%
  summarise(Query=ID,
            Sequence,
            TIGRFAM,
            InterPro) %>>%
  rbind(domain)

rSAM.domain <- lapply(rSAM$Sequence, function(x){
  temp <- domain %>>% filter(Sequence == x)
  apply(temp, 2, function(y){
    str_split(y, "[;,]") %>>%
      unlist %>>%
      unique %>>%
      sort %>>%
      paste(collapse = ";") %>>%
      str_remove("^;")
  }) %>>% t %>>% as.data.frame
})

rSAM.2 <- list.rbind(rSAM.domain) %>>% select(TIGRFAM, InterPro) %>>% (cbind(rSAM, .))


##### Annotate taxonomic lineage #####

Uniprot.taxa <- Uniprot %>>% summarise(Superkingdom, Phylum, Class, Order, Family, Genus, query= "", Sequence)

RefSeq <- rSAM %>>%
  filter(str_detect(Dataset, "Uniprot")==FALSE)
  
RefSeq.taxa.ref <- RefSeq$Organism %>>%
  str_extract_all("([A-Z][a-z]{4,} |[A-Z][a-z]{4,}$)") %>>%
  sapply(function(x){setdiff(x, "Candidatus ") %>>% str_remove(" ")})

taxa.list <- list()
load(file=paste(path1, "Intermediate files/Data setup/1-04_taxonomy reference.Rdata", sep="")) # load old 'taxa.list' to avoid duplicated search

taxa.query <- RefSeq.taxa.ref %>>%
  unlist %>>%
  unique %>>%
  sort %>>%
  setdiff(names(taxa.list))

taxa.list.add <- list()
if(length(taxa.query)>0){
for(i in 1:length(taxa.query)){
  temp <- classification(taxa.query[i],db="ncbi")[[1]]
  taxa.list.add <- c(taxa.list.add, list(temp))
  Sys.sleep(1)
}
names(taxa.list.add) <- taxa.query
}

taxa.list <- lapply(taxa.list.add, function(x){is.null(dim(x))==FALSE}) %>>%
  unlist %>>%
  (taxa.list.add[.]) %>>%
  (.[order(names(.))]) %>>%
  c(taxa.list)
save(taxa.list, file=paste(path1, "Intermediate files/Data setup/1-04_taxonomy reference.Rdata", sep=""))


rank.ref <- c("superkingdom", "phylum", "class", "order", "family", "genus")
temp <- matrix("", nrow=length(taxa.list), ncol=6) %>>%
  as.data.frame(row.names = names(taxa.list))
colnames(temp) <- rank.ref

for(i in 1:dim(temp)[1]){
  for(j in 1:length(rank.ref)){
    temp[i,j] <- taxa.list[[i]] %>>%
      filter(rank==rank.ref[j]) %>>%
      (if(dim(.)[1]>0){.$name} else {""})
  }
}

taxa.ref <- temp %>>% filter(superkingdom!="Eukaryota") %>>%
  arrange(superkingdom, phylum, class, order, family, genus) %>>%
  (mutate(., query=rownames(.)))
colnames(taxa.ref) <- c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "query")
write.xlsx(taxa.ref, paste(path1, "Intermediate files/Data setup/1-05_taxonomy reference.xlsx", sep=""))

taxa.ref <- read.xlsx(paste(path1, "Intermediate files/Data setup/1-05_taxonomy reference.xlsx", sep=""))

temp <- lapply(RefSeq.taxa.ref, function(x){intersect(x, taxa.ref$query)})
RefSeq.taxa <- lapply(temp, function(x){
  if(length(x)==1){
    taxa.ref %>>% filter(query==x)
  } else {
    data.frame(Superkingdom = "",
               Phylum = "",
               Class = "",
               Order = "",
               Family = "",
               Genus = "",
               query= "")
  }
}) %>>% list.rbind %>>% mutate(Sequence=RefSeq$Sequence)

combined.taxa <- rbind(Uniprot.taxa, RefSeq.taxa)

rSAM.taxa <- lapply(rSAM$Sequence, function(x){
  temp <- combined.taxa %>>%
    filter(Sequence == x)
  apply(temp, 2, function(y){
    str_split(y, "[;,]") %>>% unlist %>>% unique %>>% sort %>>% paste(collapse = ";") %>>% str_remove("^;")
  }) %>>%
    t %>>%
    as.data.frame
})

rSAM.3 <- list.rbind(rSAM.taxa) %>>%
  select(-query, -Sequence) %>>%
  (cbind(rSAM.2, .))
write.xlsx(rSAM.3, paste(path1, "Intermediate files/Data setup/1-06_rSAM with domain and taxonomy.xlsx", sep=""))

paste(rSAM.3$Description, "\n", rSAM.3$Sequence, sep="") %>>%
  write.table(paste(path1, "Key files/1-07_combined rSAM list for SSN.fasta", sep=""),
              quote = F, row.names = F, col.names = F)



##### Add clustering information from SSN #####

rSAM.3 <- read.xlsx(paste(path1, "Intermediate files/Data setup/1-06_rSAM with domain and taxonomy.xlsx", sep=""))
rSSN.CL <- read.xlsx(paste(path1, "Initial datasets/Node table_combined rSAM_demo.xlsx", sep=""))
#Node table exported from the sequence similarity network in cytoscape.

rSAM.3$rSAM.SSN.cluster <- ""
rSAM.3$rSAM.SSN.node <- ""

pb <- txtProgressBar(min=1, max=dim(rSAM.3)[1], style=3)
for (i in 1:dim(rSAM.3)[1]){
  temp <- grep(rSAM.3$ID[i], rSSN.CL$Description)
  rSAM.3$rSAM.SSN.cluster[i] <- rSSN.CL$Node.Count.Cluster.Number[temp][1] %>>%
    formatC(width = 2, flag = 0)
  rSAM.3$rSAM.SSN.node[i] <- rSSN.CL$shared.name[temp][1]
  setTxtProgressBar(pb, i)
}
rSAM.3$rSAM.SSN.cluster[rSAM.3$rSAM.SSN.cluster=="NA"] <- "singleton"


rSAM.4 <- rSAM.3 %>>%
  arrange(rSAM.SSN.cluster, Organism) %>>%
  mutate(Entry = paste("rSAM_demo_", formatC(1:dim(rSAM.3)[1], width = 5, flag = 0), sep="")) %>>%
  summarise(Entry,
            Accession = ID,
            rSAM.SSN.cluster,
            rSAM.SSN.node,
            Organism,
            Protein.name,
            Dataset,
            Length = str_length(Sequence),
            Sequence,
            TIGRFAM,
            InterPro,
            Superkingdom,
            Phylum,
            Class,
            Order,
            Family,
            Genus)

write.xlsx(rSAM.4, paste(path1, "Key files/1-08_Dataset S1_source table_demo.xlsx", sep=""))


##### Add domain and taxonomic information to SSN nodes #####

organise <- function(x){
  str_split(x, "[;,]") %>>% unlist %>>% unique %>>% sort %>>% paste(collapse = ";")
}

add <- list()
for(i in 1:dim(rSSN.CL)[1]){
  temp <- rSAM.4 %>>% filter(rSAM.SSN.node==rSSN.CL$shared.name[i]) 
  
  temp.2 <- temp %>>% summarise(Superkingdom = organise(Superkingdom),
                                Phylum = organise(Phylum),
                                Class = organise(Class),
                                Order = organise(Order),
                                Family = organise(Family),
                                Genus = organise(Genus),
                                TIGRFAM = organise(TIGRFAM),
                                InterPro = organise(InterPro),
                                Dataset = organise(Dataset))
  add <- c(add, list(temp.2))
}
rSSN.CL.2 <- list.rbind(add) %>>% (cbind(rSSN.CL, .))
rSSN.CL.2 %>>% select(-Description) %>>%
  write.csv(paste(path1, "Intermediate files/Data setup/1-09_Node table_combined rSAM_demo_modified.csv", sep=""))
