# Sugiyama et al., JACS
# Script-implemented precursor search
# By Ryosuke Sugiyama, PhD; zy105shsh@gmail.com
# Last update: June 8, 2022


rm(list = ls()) 
library(openxlsx)
library(pipeR)
library(rlist)
library(tidyverse)

library(patchwork)
library(staplr)
library(rentrez)
library(UniprotR)
library(seqinr)
library(taxize)
library(msa)
library(ggmsa)
library(ggseqlogo)
library(scales) # Utilities for scales and formatting
#set_entrez_key("xxx") # Activate API key


path1 <- 'xxx/JACS2022/Demo/' #set directory

prec <- read.xlsx(paste(path1, "Intermediate files/Manual analysis/4-03_Precursors after manual analysis_wo duplicates.xlsx", sep=""))
criteria <- read.xlsx(paste(path1, "Key files/3-01_precursor criteria.xlsx", sep=""))
taxa.ref.all <- read.xlsx(paste(path1, "Intermediate files/Data setup/1-05_taxonomy reference.xlsx", sep=""))
load(file=paste(path1, "Intermediate files/Data setup/1-04_taxonomy reference.Rdata", sep=""))

taxa.query <- prec$Strain %>>%
  str_extract_all("([A-Z][a-z]{4,} |[A-Z][a-z]{4,}$)") %>>%
  sapply(function(x){setdiff(x, "Candidatus ") %>>% str_remove(" ") %>>%
      intersect(taxa.ref.all$query)})

prec.no.taxa <- lapply(taxa.query, length) %>>% unlist %>>% (prec[.==0,])
taxa.add <- prec.no.taxa$Strain %>>%
  str_extract_all("([A-Z][a-z]{4,} |[A-Z][a-z]{4,}$)") %>>%
  sapply(function(x){setdiff(x, "Candidatus ") %>>% str_remove(" ")}) %>>% unlist %>>% unique

taxonomy.add <- list()
for(i in 1:length(taxa.add)){
  temp <- classification(taxa.add[i],db="ncbi")[[1]]
  taxonomy.add <- c(taxonomy.add, list(temp))
}
names(taxonomy.add) <- taxa.add
save(taxonomy.add, file=paste(path1, "Intermediate files/Data finalization/5-01_tax ref_add.Rdata", sep=""))
c(taxa.list, taxonomy.add) %>>%
save(file=paste(path1, "Intermediate files/Data finalization/5-02_tax ref_v2.Rdata", sep=""))


load(file=paste(path1, "Intermediate files/Data finalization/5-01_tax ref_add.Rdata", sep=""))
taxonomy.add <- taxonomy.add[is.na(taxonomy.add)==FALSE]
rank.ref <- c("superkingdom", "phylum", "class", "order", "family", "genus")
temp <- matrix("", nrow=length(taxonomy.add), ncol=6) %>>% as.data.frame(row.names = names(taxonomy.add))
colnames(temp) <- rank.ref

for(i in 1:dim(temp)[1]){
  for(j in 1:length(rank.ref)){
    temp[i,j] <- taxonomy.add[[i]] %>>% filter(rank==rank.ref[j]) %>>% (if(dim(.)[1]>0){.$name} else {""})
  }
}
taxa.ref <- temp %>>% filter(superkingdom!="Eukaryota") %>>%
  (mutate(., query=rownames(.)))
colnames(taxa.ref) <- c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "query")

taxa.ref.all <- rbind(taxa.ref.all, taxa.ref) %>>%
  arrange(Superkingdom, Phylum, Class, Order, Family, Genus) 
write.xlsx(taxa.ref.all, paste(path1, "Intermediate files/Data finalization/5-03_taxonomy reference_v2.xlsx", sep=""))



taxa.ref.all <- read.xlsx(paste(path1, "Intermediate files/Data finalization/5-03_taxonomy reference_v2.xlsx", sep=""))
taxa.query <- prec$Strain %>>%
  str_extract_all("([A-Z][a-z]{4,} |[A-Z][a-z]{4,}$)") %>>%
  sapply(function(x){setdiff(x, "Candidatus ") %>>% str_remove(" ") %>>% intersect(taxa.ref.all$query)})

prec.taxa <- lapply(taxa.query, function(x){
  if(length(x)==1){
    taxa.ref.all %>>% filter(query==x)
  } else {
    data.frame(Superkingdom = "",
               Phylum = "",
               Class = "",
               Order = "",
               Family = "",
               Genus = "",
               query= "")
  }
}) %>>% list.rbind

prec.2 <- cbind(prec, prec.taxa)
write.xlsx(prec.2, paste(path1, "Intermediate files/Data finalization/5-04_precursors with taxonomic lineage.xlsx", sep=""))

##### Add clustering information from precursor SSN #####

prec.2 <- read.xlsx(paste(path1, "Intermediate files/Data finalization/5-04_precursors with taxonomic lineage.xlsx", sep=""))
CL <- read.xlsx(paste(path1, "Additional datasets for finalization/Node table_combined precursor_demo.xlsx", sep=""))
rSAM <- read.xlsx(paste(path1, "Intermediate files/GenBank processing/2-04_rSAM without redundancy.xlsx", sep=""))

prec.2$prec.SSN.cluster <- ""
prec.2$prec.SSN.node <- ""

pb <- txtProgressBar(min=1, max=dim(prec.2)[1], style=3)
for (i in 1:dim(prec.2)[1]){
  temp <- grep(prec.2$Protein.ID[i], CL$Description)
  prec.2$prec.SSN.cluster[i] <- CL$Node.Count.Cluster.Number[temp][1] %>>% formatC(width = 3, flag = 0)
  prec.2$prec.SSN.node[i] <- CL$shared.name[temp][1]
  setTxtProgressBar(pb, i)
}
prec.2$prec.SSN.cluster[str_detect(prec.2$prec.SSN.cluster,"NA")] <- "singleton"

write.xlsx(prec.2, paste(path1, "Intermediate files/Data finalization/5-05_precursors with SSN information.xlsx", sep=""))


##### Add precursor information to SSN nodes #####

prec.3 <-  read.xlsx(paste(path1, "Intermediate files/Data finalization/5-05_precursors with SSN information.xlsx", sep=""))
CL <- read.xlsx(paste(path1, "Additional datasets for finalization/Node table_combined precursor_demo.xlsx", sep=""))

organise <- function(x){
  str_split(x, "[;,]") %>>% unlist %>>% unique %>>% sort %>>% paste(collapse = ";")
}

add <- list()
for(i in 1:dim(CL)[1]){
  temp <- prec.3 %>>% filter(prec.SSN.node==CL$shared.name[i]) 
  
  temp.2 <- temp %>>% summarise(Superkingdom = organise(Superkingdom),
                                Phylum = organise(Phylum),
                                Class = organise(Class),
                                Order = organise(Order),
                                Family = organise(Family),
                                Genus = organise(Genus),
                                prec.class = organise(Note),
                                rSAM.query = organise(rSAM.query),
                                rSAM.SSN.cluster = organise(rSAM.SSN.cluster))
  add <- c(add, list(temp.2))
}
CL.2 <- list.rbind(add) %>>% (cbind(CL, .))
CL.2 %>>% select(-Description) %>>%
  write.csv(paste(path1, "Intermediate files/Data finalization/5-06_Node table_combined precursor_demo_modified.csv", sep=""))


CL.r <- read.csv(paste(path1, "Intermediate files/Data setup/1-09_Node table_combined rSAM_demo_modified.csv", sep=""))
CL.r$prec <- ""
entry <- lapply(CL.r$shared.name, function(x){
  rSAM$Entry[rSAM$rSAM.SSN.node==x] %>>% unique
})

temp.2 <- lapply(entry, function(x){
  temp <- lapply(x, function(y){
    prec.3 %>>% filter(str_detect(Entry,y))
  }) %>>% list.rbind
  if(dim(temp)[1]>0){
    data.frame(prec.class=organise(temp$Note),
               prec.node=organise(temp$prec.SSN.node),
               prec.SSN.cluster=organise(temp$prec.SSN.cluster),
               prec.ID = organise(temp$Protein.ID))
  } else {
    data.frame(prec.class="",
               prec.node="",
               prec.SSN.cluster="",
               prec.ID="")
  }
}) %>>% list.rbind
cbind(CL.r, temp.2) %>>%
  write.csv(paste(path1, "Intermediate files/Data finalization/5-07_Node table_combined rSAM_demo_with precursor.csv", sep=""))






##### Finalize source data of Dataset S3 #####


prec <-  read.xlsx(paste(path1, "Intermediate files/Data finalization/5-05_precursors with SSN information.xlsx", sep=""))
prec$Note[str_detect(prec$Note, "putative")] <- "Putative"

prec$c <- ""
for(i in 1:dim(prec)[1]){
  prec$c[i] <- if(prec$rSAM.SSN.cluster[i]=="03"|prec$rSAM.SSN.cluster[i]=="08"){
    str_sub(prec$Sequence[i], -40, -1)
  } else {
    str_sub(prec$Sequence[i], -20, -1)
  }
}

prec$num <- 1:dim(prec)[1]

criteria <- data.frame(label = c("Fxs","Grr","Xye","Dar"),
                       motif.pattern = c("[F].[ANS].{2,3}$",
                                         "[WF].[ADRNQ]",
                                         "[WFY].[ARNDEQKS].[WFY].[ARNDEQKS][WFY].[ARNDEQKS]",
                                         "[W].[WFY]..")) 


prec$Note <- as.character(prec$Note)

prec$motif <- ""
for(j in 1:dim(criteria)[1]){
  temp <- prec %>>% filter(Note==criteria$label[j])
  if(dim(temp)[1]>0){
  temp2 <- str_locate_all(temp$c, criteria$motif.pattern[j])
  for(i in 1:dim(temp)[1]){
    t <- temp2[[i]] %>>% as.data.frame
    t2 <- sapply(1:dim(t)[1], function(x){
      str_sub(temp$c[i], start=t$start[x], end=t$end[x])
    }) %>>% unlist %>>% unique %>>% paste(collapse = " ")
    prec$motif[prec$num == temp$num[i]] <- t2
  }
}
}

prec2 <- prec %>>% filter(motif=="")

criteria2 <- data.frame(label = c("Yxd","Haa","Yhh"),
                        prec.clust = c("001","013","004"),
                        motif.pattern = c("Y.D",
                                          "HAA",
                                          "[Y]...[H]..[H].."))

for(j in 1:dim(criteria2)[1]){
  temp <- prec2 %>>% filter(str_detect(prec.SSN.cluster, criteria2$prec.clust[j]))
  if(dim(temp)[1]>0){
  temp2 <- str_locate_all(temp$c, criteria2$motif.pattern[j])
  for(i in 1:dim(temp)[1]){
    if(dim(temp2[[i]])[1]>0){
      t <- temp2[[i]] %>>% as.data.frame
      t2 <- sapply(1:dim(t)[1], function(x){
        str_sub(temp$c[i], start=t$start[x], end=t$end[x])
      }) %>>% unlist %>>% unique %>>% paste(collapse = " ")
      prec$motif[prec$num == temp$num[i]] <- t2
      prec$Note[prec$num == temp$num[i]] <- criteria2$label[j]
    }
  }
}
}

prec3 <- prec %>>% filter(motif=="")
criteria3 <- data.frame(label = "Putative",
                        cluster = ".{1,10}",
                        motif.pattern = "[WHFY].[ARNDEQKSL]")

for(j in 1:dim(criteria3)[1]){
  temp <- prec3 %>>% filter(str_detect(rSAM.SSN.cluster, criteria3$cluster[j]))
  if(dim(temp)[1]>0){
  temp2 <- str_locate_all(temp$c, criteria3$motif.pattern[j])
  for(i in 1:dim(temp)[1]){
    if(dim(temp2[[i]])[1]>0){
      t <- temp2[[i]] %>>% as.data.frame
      t2 <- sapply(1:dim(t)[1], function(x){
        str_sub(temp$c[i], start=t$start[x], end=t$end[x])
      }) %>>% unlist %>>% unique %>>% paste(collapse = " ")
      prec$motif[prec$num == temp$num[i]] <- t2
      prec$Note[prec$num == temp$num[i]] <- criteria3$label[j]
    }
  }
  }
}
prec$Note <- factor(prec$Note, levels=c("Fxs", "Grr", "Xye", "Dar", "Yxd", "Haa", "Yhh", "Putative"))

prec.4 <- prec %>>% 
  arrange(prec.SSN.cluster, Note, motif, prec.SSN.node, Strain, Protein.ID) %>>%
  summarise(Cp = prec.SSN.cluster,
            Cr = rSAM.SSN.cluster,
            Associated.rSAM = rSAM.query,
            Precursor.class = Note,
            Potential.motif = motif,
            Sequence,
            Length,
            Protein.ID,
            Product,
            Organism = Strain,
            Source,
            Locus.tag)

write.xlsx(prec.4, paste(path1, "Key files/5-08_Table S3_source table_demo.xlsx", sep=""))


##### Finalize source data of Dataset S2 #####

prec.4 <- read.xlsx(paste(path1, "Key files/5-08_Table S3_source table_demo.xlsx", sep=""))
rSAM <- read.xlsx(paste(path1, "Intermediate files/GenBank processing/2-04_rSAM without redundancy.xlsx", sep=""))

rSAM$prec <- ""
rSAM$prec.SSN <- ""
rSAM$prec.class <- ""

for(i in 1:dim(rSAM)[1]){
  temp <- prec.4$Associated.rSAM == rSAM$Accession[i]
  if(sum(temp)>0){
    rSAM$prec[i] <- organise(prec.4$Protein.ID[temp])
    rSAM$prec.SSN[i] <- organise(prec.4$Cp[temp])
    rSAM$prec.class[i] <- organise(prec.4$Precursor.class[temp])
  }
}

label <- c("Fxs", "Grr", "Xye", "Dar", "Yxd", "Haa", "Yhh") %>>%
  (data.frame(a=.,b=paste(.,";Putative", sep=""))) %>>%
  t %>>%
  as.vector %>>%
  c("Putative", "")
rSAM$prec.class <- factor(rSAM$prec.class, levels=label)
rSAM$prec.SSN <- factor(rSAM$prec.SSN, levels=unique(rSAM$prec.SSN) %>>% sort %>>% (c(.[-1], "")))
rSAM.2 <- rSAM %>>% arrange(Obsolete,rSAM.SSN.cluster,prec.SSN, prec.class, Organism)

rSAM.3 <-  rSAM.2 %>>% summarise(Entry,
                           Cr = rSAM.SSN.cluster,
                           Accession,
                           Cp = prec.SSN,
                           Precursor.ID = prec,
                           Precursor.class = prec.class,
                           Organism,
                           Protein.name,
                           Sequence,
                           Length,
                           Nucleotide.reference,
                           Nucleotides,
                           No.of.nucleotides,
                           Dataset,
                           TIGRFAM,
                           InterPro,
                           Superkingdom,
                           Phylum,
                           Class,
                           Order,
                           Family,
                           Genus,
                           Obsolete)
  
write.xlsx(rSAM.3, paste(path1, "Key files/5-09_Table S2_source table_demo.xlsx", sep=""))


