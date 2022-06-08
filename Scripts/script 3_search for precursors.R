# Sugiyama et al., JACS
# Script-implemented precursor search
# By Ryosuke Sugiyama, PhD; zy105shsh@gmail.com
# Last update: June 8, 2022


rm(list = ls()) 
library(openxlsx)
library(pipeR)
library(rlist)
library(tidyverse)

path1 <- 'xxx/JACS2022/Demo/' #set directory

df1.pre <- read.xlsx(paste(path1, "Intermediate files/GenBank processing/2-03_Combined genetic loci_demo.xlsx", sep=""))
df1.pre$Entry <- factor(df1.pre$Entry, levels=unique(df1.pre$Entry))
df1 <- df1.pre %>>% (split(., factor(.$Entry)))

criteria <- data.frame()

# criteria Fxs
criteria <- data.frame(label = "Fxs",
                       priority = 1,
                       dist.min = 1,
                       dist.max = 1,
                       len.min = 1,
                       len.max = 120,
                       motif.start = -10,
                       motif.end = -1,
                       motif.No.min = 1,
                       motif.No.max = 20,
                       motif.pattern = "[F].[ANS].{2,3}$",
                       strand = "any",
                       cluster = "02",
                       other = TRUE) %>>%
  (rbind(criteria, .))

# criteria Grr
criteria <- data.frame(label = "Grr",
                       priority = 1,
                       dist.min = 1,
                       dist.max = 5,
                       len.min = 100,
                       len.max = 300,
                       motif.start = -40,
                       motif.end = -1,
                       motif.No.min = 3,
                       motif.No.max = 20,
                       motif.pattern = "[WF].[ADRNQ]",
                       strand = "any",
                       cluster = "03",
                       other = "str_sub(df$seq, criteria$motif.start[j], end=criteria$motif.end[j]) %>>% replace_na(FALSE) %>>% str_count(pattern = \"GG\")>=2") %>>%
  (rbind(criteria, .))


# criteria Xye
criteria <- data.frame(label = "Xye",
                       priority = 1,
                       dist.min = 1,
                       dist.max = 1,
                       len.min = 1,
                       len.max = 100,
                       motif.start = -15,
                       motif.end = -1,
                       motif.No.min = 1,
                       motif.No.max = 1,
                       motif.pattern = "[WFY].[ARNDEQKS].[WFY].[ARNDEQKS][WFY].[ARNDEQKS]",
                       strand = "clustered",
                       cluster = "05",
                       other = TRUE) %>>%
  (rbind(criteria, .))

# criteria Dar
criteria <- data.frame(label = "Dar",
                       priority = 1,
                       dist.min = 1,
                       dist.max = 5,
                       len.min = 1,
                       len.max = 150,
                       motif.start = -20,
                       motif.end = -1,
                       motif.No.min = 1,
                       motif.No.max = 20,
                       motif.pattern = "[W].[WFY]..",
                       strand = "clustered",
                       cluster = "06",
                       other = TRUE) %>>%
  (rbind(criteria, .))


# criteria putative.1 (adjacent)
criteria <- data.frame(label = "putative-1",
                       priority = 2,
                       dist.min = 1,
                       dist.max = 1,
                       len.min = 1,
                       len.max = 120,
                       motif.start = -20,
                       motif.end = -1,
                       motif.No.min = 1,
                       motif.No.max = 20,
                       motif.pattern = "[WHFY].[ARNDEQKSL]",
                       strand = "any",
                       cluster = ".",
                       other = TRUE) %>>%
  (rbind(criteria, .))


# criteria putative.2 (clustered)
criteria <- data.frame(label = "putative-2",
                       priority = 3,
                       dist.min = 2,
                       dist.max = 5,
                       len.min = 1,
                       len.max = 120,
                       motif.start = -20,
                       motif.end = -1,
                       motif.No.min = 1,
                       motif.No.max = 20,
                       motif.pattern = "[WHFY].[ARNDEQKSL]",
                       strand = "clustered",
                       cluster = ".",
                       other = TRUE) %>>%
  (rbind(criteria, .))

# criteria putative.3 (wide)
criteria <- data.frame(label = "putative-3",
                       priority = 4,
                       dist.min = 2,
                       dist.max = 5,
                       len.min = 1,
                       len.max = 120,
                       motif.start = -20,
                       motif.end = -1,
                       motif.No.min = 1,
                       motif.No.max = 20,
                       motif.pattern = "[WHFY].[ARNDEQKSL]",
                       strand = "any",
                       cluster = ".",
                       other = TRUE) %>>%
  (rbind(criteria, .))

write.xlsx(criteria, paste(path1, "Key files/3-01_precursor criteria.xlsx", sep=""))


for (j in 1:dim(criteria)[1]){
  prec.list <- data.frame()
  
  t <- 1
  pb <- txtProgressBar(min=0, max=length(df1), style=3)
  for (i in 1:length(df1)){
    df <- df1[[i]]
    if(dim(df)[1]>1 & sum(str_detect(df$rSAM.SSN.cluster, criteria$cluster[j]))>0){
      n <- 1:dim(df)[1]
      rSAM <- grep("rSAM", df$Note)
      D <- sapply(rSAM, function(x){abs(n-x)}) %>>% apply(1,min)
      L <- df$Length
      M <- str_sub(df$Sequence, start = criteria$motif.start[j], end=criteria$motif.end[j]) %>>%
        str_extract_all(pattern = criteria$motif.pattern[j])
      df$Motifs <- lapply(M, function(x){paste(x, collapse = "; ")}) %>>% unlist
      M.no <- lapply(M, length) %>>% unlist
      S <- str_detect(df$Note, "[:alpha:]") %>>% replace_na(FALSE)
      Prec <- rbind(D>=criteria$dist.min[j],
                    D<=criteria$dist.max[j],
                    L>=criteria$len.min[j],
                    L<=criteria$len.max[j],
                    M.no>=criteria$motif.No.min[j],
                    M.no<=criteria$motif.No.max[j],
                    eval(parse(text=criteria$other[j])),
                    if(criteria$strand[j]=="any"){TRUE} else if(criteria$strand[j]=="clustered"){S}) %>>%
        (apply(.,2,sum)==dim(.)[1])
      df$Note[Prec] <- criteria$label[j]
      
      if(sum(Prec)>0){
        pr <- df[Prec,] %>>% summarise(Entry,
                                       rSAM.SSN.cluster,
                                       rSAM.query,
                                       rSAM.product = df$Product[rSAM],
                                       Note,
                                       Priority = criteria$priority[j],
                                       Product,
                                       Sequence,
                                       Length,
                                       Motifs,
                                       Strain,
                                       Source,
                                       Locus.tag,
                                       Protein.ID,
                                       FASTA,
                                       ID = paste(Locus.tag, Sequence, sep="_"))
        
        prec.list <- rbind(prec.list, pr)
      }
    }
    setTxtProgressBar(pb, t)
    t <- t+1
  }
  if(dim(prec.list)[1]>0){
    prec.list.2 <- prec.list %>>% arrange(rSAM.SSN.cluster, Motifs, Sequence)
    
    write.xlsx(prec.list.2, paste(path1, "Intermediate files/Precursor search/Partial data/precursor_type ", j, "_", criteria$label[j], ".xlsx", sep=""),
               overwrite = TRUE)
  }
}

##### Combine precursor lists #####

part <- list.files(paste(path1, "Intermediate files/Precursor search/Partial data/", sep="")) %>>% (.[str_detect(.,"^precursor_type")])

prec.list <- read.xlsx(paste(path1, "Intermediate files/Precursor search/Partial data/", part[1], sep=""))[1,]
for (j in 1:length(part)){
  xlsx <- read.xlsx(paste(path1, "Intermediate files/Precursor search/Partial data/", part[j], sep=""))
  pb <- txtProgressBar(min=0, max=dim(xlsx)[1], style=3)
  t <- 0
  for (i in 1:dim(xlsx)[1]){
    temp <- grep(paste("^", xlsx$ID[i], "$",sep=""), prec.list$ID)
    if(length(temp)==0){
      prec.list <- rbind(prec.list, xlsx[i,])
    }
    setTxtProgressBar(pb, t)
    t <- t+1
  }
}
prec.list <- prec.list %>>% arrange(rSAM.SSN.cluster, Strain, Entry, Priority, Locus.tag)
write.xlsx(prec.list, paste(path1, "Intermediate files/Precursor search/3-02_putative precursors_all.xlsx", sep=""))


##### Remove duplicated precursors associated with redundant rSAMs #####

prec.list <- read.xlsx(paste(path1, "Intermediate files/Precursor search/3-02_putative precursors_all.xlsx", sep=""))
prec.list <- prec.list %>>% mutate(f= paste(Strain, Sequence, sep="_"))
duplicates <- lapply(unique(prec.list$f), function(x){
  temp <- prec.list %>>% filter(f==x)
  if(dim(temp)[1]>1){
    temp <- temp %>>% mutate(RefSeq = str_detect(temp$Protein.ID, "^WP_"))
    temp %>>% arrange(desc(RefSeq))
  }
}) %>>% list.rbind
write.xlsx(duplicates, paste(path1, "Intermediate files/Precursor search/3-03_duplicated precursors.xlsx", sep=""))

# Check duplicates manually and save as '3-03_duplicated precursors_checked.xlsx'.

duplicates <- read.xlsx(paste(path1, "Intermediate files/Precursor search/3-03_duplicated precursors_checked.xlsx", sep=""),
                        sheet = 1) %>>% filter(RefSeq==FALSE)

prec.list.wo.duplicates <- sapply(duplicates$ID, function(x){
  grep(x, prec.list$ID)
}) %>>% unlist %>>% (prec.list[-.,])
write.xlsx(prec.list.wo.duplicates, paste(path1, "Intermediate files/Precursor search/3-04_putative precursors_wo duplicates.xlsx", sep=""))


redundant.rSAM <- setdiff(duplicates$rSAM.query, prec.list.wo.duplicates$rSAM.query)
rSAM <- read.xlsx(paste(path1, "Intermediate files/GenBank processing/2-02_rSAM for precursor search_loaded.xlsx", sep=""))
for(i in redundant.rSAM){
  rSAM$Obsolete[rSAM$Accession==i] <- "redundant"
  df1.pre$Check[df1.pre$rSAM.query==i] <- "redundant"
}
rSAM %>>% arrange(Obsolete, Entry) %>>%
  write.xlsx(paste(path1, "Intermediate files/GenBank processing/2-04_rSAM without redundancy.xlsx", sep=""))
df1.pre %>>% arrange(Check, Entry, Start) %>>%
  write.xlsx(paste(path1, "Intermediate files/GenBank processing/2-05_Combined genetic loci_demo_without redundancy.xlsx", sep=""))



##### Classify precursors based on priority #####

prec.list.wo.duplicates <- read.xlsx(paste(path1, "Intermediate files/Precursor search/3-04_putative precursors_wo duplicates.xlsx", sep="")) %>>%
  select(-f, -ID)

priority <- prec.list.wo.duplicates %>>%
  group_by(Entry) %>>%
  summarise(Class=min(Priority),
            Number=length(Priority),
            rSAM.SSN.cluster=unique(rSAM.SSN.cluster)[1]) %>>%
  arrange(Class, rSAM.SSN.cluster, desc(Number), Entry)


precursor.class <- list()
for (i in 1:4){
  temp <- priority %>>% filter(Class==i)
  temp.2 <- lapply(temp$Entry, function(x){
    prec.list.wo.duplicates %>>% filter(Entry==x) %>>% arrange(Priority, Locus.tag)
  }) %>>% list.rbind
  precursor.class <- c(precursor.class, list(temp.2))
}
names(precursor.class) <- paste("priority", 1:4)
write.xlsx(precursor.class, paste(path1, "Intermediate files/Precursor search/3-05_precursor classified by priority_demo.xlsx", sep=""))






