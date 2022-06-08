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

path1 <- 'xxx/JACS2022/Demo/' #set directory

df1.pre <- read.xlsx(paste(path1, "Intermediate files/GenBank processing/2-05_Combined genetic loci_demo_without redundancy.xlsx", sep=""))
remove <- df1.pre %>>% filter(Check=="redundant")
precursor <- read.xlsx(paste(path1, "Intermediate files/Precursor search/3-04_putative precursors_wo duplicates.xlsx", sep=""))
criteria <- read.xlsx(paste(path1, "Key files/3-01_precursor criteria.xlsx", sep=""))
priority <- lapply(1:5, function(x){
  if(x!=5){
    temp <- read.xlsx(paste(path1, "Intermediate files/Precursor search/3-05_precursor classified by priority_demo.xlsx", sep=""), sheet = x)
    unique(temp$Entry)
  } else {
    setdiff(df1.pre$Entry, precursor$Entry) %>>% setdiff(remove$Entry)
  }
})

##### modify genetic loci for manual analysis #####

note <- c("", "clustered", "rSAM", criteria$label, "pseudo")
priority.color <- c(rep("blue", 4), "green", "red", "magenta")
gene.fill <- c("grey90", "grey15",  "yellow", priority.color, "white")

df.mod <- list()
for(n in 1:5){
type <- paste("priority", n)

df.ref <- lapply(priority[[n]], function(x){
  df1.pre %>>% filter(Entry==x)
})
names(df.ref) <- priority[[n]]

pb <- txtProgressBar(min=0, max=length(df.ref), style=3)
df.ref.mod <- list()
df.g <- list()
i <- 1
for (i in 1:length(df.ref)){
  df <- df.ref[[i]]
  df$precursor <- ""
  rSAM <- grep("rSAM", df$Note)[1]
  if(n!=5){
  prec <- precursor %>>% filter(Entry==df$Entry[1]) %>>% select(Note, Priority, Protein.ID)
  for (j in rev(1:dim(prec)[1])){
    temp <- grep(prec$Protein.ID[j], df$Protein.ID)
    df$precursor[temp] <- prec$Note[j]
    df$Note[temp] <- prec$Note[j]
    df$Check[temp] <- "f"
  }
  if(n == 1){
    df$Check[str_detect(df$precursor, "(Fxs|Grr|Xye|Dar)")] <- "p"
  }
  if(n == 2){
    df$Check[df$precursor=="putative-1"] <- "p"
  }
  }
  
  if(df$Strand[rSAM]=="+"){
    df$Start.g <- df$Start-df$Start[rSAM]
    df$End.g <- df$End-df$Start[rSAM]
  } else if(df$Strand[rSAM]=="-"){
    df$Start.g <- -(df$End-df$End[rSAM])
    df$End.g <- -(df$Start-df$End[rSAM])
    df$Strand <- str_detect(df$Strand, pattern = "\\+") %>>%
      str_replace_all(pattern=c("TRUE"="-", "FALSE"="+"))
    df <- df[rev(rownames(df)),]
  }
  df.ref.mod <- select(df, -Start.g, -End.g, -precursor) %>>% (c(df.ref.mod, list(.)))
  
  df <- df %>>% select(Entry, Strain, Protein.ID, Note, precursor, Strand, Start.g, End.g, Product, Sequence)
  df.g.pre <- list()
  for(p in 1:dim(df)[1]){
    pre <- df[p,]
    if(pre$Strand=="+"){
      x <- c(pre$Start, max(pre$Start,pre$End-300),pre$End,max(pre$Start,pre$End-300),pre$Start)
      y <- c(0.25,0.25,1.25,2.25,2.25)
    } else {
      x <- c(pre$Start, min(pre$Start+300,pre$End),pre$End,pre$End, min(pre$Start+300,pre$End))
      y <- c(-1.25,-2.25,-2.25,-0.25,-0.25)
    }
    if(pre$Protein.ID=="pseudo"){
      pre$Protein.ID <- paste("pseudo", pre$Start.g, sep=".")
      pre$Note <- "pseudo"
    }
    df.g.pre <- data.frame(Entry=pre$Entry, ID=paste(pre$Protein.ID,p,sep="."), x, y, Note=pre$Note, Precursor=pre$precursor) %>>%
      (c(df.g.pre, list(.)))
  }
  df.g.pre <- list.rbind(df.g.pre)
  df.g.pre[is.na(df.g.pre)] <- ""
  df.g.pre$fill <- "grey90"
  df.g.pre$lty <- "solid"
  df.g.pre$lty[df.g.pre$Note=="pseudo"] <- "dotted"
  for (j in 1:length(note)){
  df.g.pre$fill[grep(note[j], df.g.pre$Note)] <- gene.fill[j]
  df.g.pre$fill[grep(note[j], df.g.pre$Precursor)] <- gene.fill[j]
  }
  df.g <- c(df.g, list(df.g.pre))
  setTxtProgressBar(pb, i)
}
names(df.g) <- names(df.ref)
df.ref.mod <- list.rbind(df.ref.mod) %>>% mutate(Note=replace_na(Note, ""))
df.mod <- c(df.mod, list(df.ref.mod))
save(df.g, file=paste(path1, "Intermediate files/Manual analysis/Visualized genetic loci/partial/Genetic loci for visualization_", type, ".Rdata", sep=""))
paste(path1, "Intermediate files/Manual analysis/Genetic loci_", type, "_manual analysis.xlsx", sep="") %>>%
  (write.xlsx(df.ref.mod, ., overwrite = TRUE))
}



##### visualize genetic loci for manual analysis #####
n <- 1
for(n in 1:5){
  type <- paste("priority", n)
  load(file=paste(path1, "Intermediate files/Manual analysis/Visualized genetic loci/partial/Genetic loci for visualization_", type, ".Rdata", sep=""))
  
  title.a <- vector()
  title.b <- vector()
  
  if(n!=5){
    precursor <- read.xlsx(paste(path1, "Intermediate files/Precursor search/3-05_precursor classified by priority_demo.xlsx", sep=""), sheet = n)
      
    seq <- str_sub(precursor$Sequence, -20, -1)
    motif <- precursor$Motifs %>>% str_split("; ") %>>% lapply(unique)
    pb <- txtProgressBar(min=0, max=dim(precursor)[1], style=3)
    
    for(x in 1:dim(precursor)[1]){
      a <- seq[x]
      for(i in motif[[x]]){
        if(!is.na(i))
        a <- str_replace_all(a, i, paste("\"*underline(", i, ")*\"", sep=""))%>>%  str_replace_all("\\*\"\"","")
      }
      a <- paste(precursor$Protein.ID[x], ": ...", a," (", precursor$Length[x], " aa)", sep="")
      str_sub(a, 1,1) <- paste("paste(\"", str_sub(a, 1,1), sep="")
      str_sub(a, -1,-1) <- paste(str_sub(a, -1,-1), "\", sep = \"\")", sep="")
      b <- paste(precursor$Entry[x], " (Cr-", precursor$rSAM.SSN.cluster[x], ")\n",
                 precursor$Strain[x], sep="")
      str_sub(b, 1,1) <- paste("paste(\"", str_sub(b, 1,1), sep="")
      str_sub(b, -1,-1) <- paste(str_sub(b, -1,-1), "\", sep = \"\")", sep="")
      title.a <- c(title.a, a)
      title.b <- c(title.b, b)
      title.c <- rep("black", dim(precursor)[1])
      if(n == 1){
      title.c[str_detect(precursor$Note, "(Fxs|Grr|Xye|Dar)")] <- "blue"
      }
      if(n == 2){
        title.c[precursor$Note == "putative-1"] <- "darkgreen"
      }
      setTxtProgressBar(pb, x)
    }
  } else {
    for(x in 1:length(df.g)){
      b <- paste(df.g[[x]]$Entry[1], "\n",
                 df.g[[x]]$Strain[1], sep="")
      str_sub(b, 1,1) <- paste("paste(\"", str_sub(b, 1,1), sep="")
      str_sub(b, -1,-1) <- paste(str_sub(b, -1,-1), "\", sep = \"\")", sep="")
      title.a <- c(title.a, "")
      title.b <- c(title.b, b)
      title.c <- rep("black", length(df.g))
      setTxtProgressBar(pb, x)
    }
  }
  
  pack <- 1:length(title.b) %>% split(ceiling(seq_along(.)/20))
  pb <- txtProgressBar(min=0, max=length(pack), style=3)
  k <- 1
  x <- 1
  i <- 1
  r <- c(-14000, 16000)
  for (k in 1:length(pack)){
    if(n==5){
    graph.e <- lapply(1:length(pack[[k]]), function(x){
      df.g[names(df.g)==names(df.g)[pack[[k]]][x]][[1]] %>>% data.frame(num=x)
    }) %>>% list.rbind
    } else {
    graph.e <- lapply(1:length(pack[[k]]), function(x){
      df.g[names(df.g)==precursor$Entry[pack[[k]]][x]][[1]] %>>% data.frame(num=x)
    }) %>>% list.rbind
    }
    graph.e[is.na(graph.e)] <- ""
    graph.e$num <- factor(graph.e$num, levels = unique(graph.e$num))
    range <- tapply(graph.e$x, graph.e$num, range) %>>% list.rbind %>>%
      (data.frame(ta=title.a[pack[[k]]],
                  tb=title.b[pack[[k]]],
                  m=.[,1],
                  M=.[,2],
                  num=levels(graph.e$num),
                  x=0,
                  y=0,
                  col=title.c[pack[[k]]]))
    range$num <- factor(range$num, levels = range$num)
    
    for(i in 1:length(pack[[k]])){
      temp <- graph.e %>>% filter(num==i & str_detect(ID, precursor$Protein.ID[pack[[k]]][i]))
      range$x[i] <- ifelse(temp$y[1]>0, range(temp$x)[1], range(temp$x)[2])
      range$y[i] <- ifelse(temp$y[1]>0, 2.6, -2.6)
      range$hjust[i] <- ifelse(temp$y[1]>0, 0.3, 0.7)
    }
    
    gd <- ggplot()+
      geom_segment(data=graph.e, aes(x=r[1],xend=r[2],y=0,yend=0), color="grey80", size=0.10, linetype= "dashed") +
      geom_segment(data=range, aes(x=m,xend=M,y=0,yend=0), color="grey20", size=0.2) +
      geom_polygon(data=graph.e, aes(x=x, y=y, group=ID), fill=graph.e$fill, color="black", size=0.15, linetype="solid") +
      geom_text(data=range, aes(x=x, y=y, hjust=hjust), color="black", label="*",
                vjust=0.8, size=(10*5/14)) + 
      facet_grid(num~.)+
      scale_y_continuous(limits = c(-3,3), expand = c(0,0)) +
      scale_x_continuous(limits = r, expand = c(0.01,0.01)) +
      theme_classic() +
      labs(x=NULL,y=NULL)
    gd <- gd + theme(axis.line = element_blank(), 
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     strip.text = element_blank(),
                     legend.position = "none",
                     panel.background = element_rect(fill="transparent", colour=NA),
                     plot.background=element_rect(fill="transparent", colour=NA),
                     strip.background=element_rect(fill="transparent", colour=NA))
    
    
    gt <- ggplot()+
      geom_text(data=range, aes(label=ta),x=0, y=-1, parse = TRUE,
                hjust=0,vjust=1, size=(5.5*5/14), color=range$col) + 
      geom_text(data=range, aes(label=tb), x=0, y=-0.5, parse = TRUE,
                hjust=0,vjust=0, size=(5.5*5/14), color="black") + 
      facet_grid(num~.)+
      scale_y_continuous(limits = c(-3,3), expand = c(0,0)) +
      scale_x_continuous(expand = c(0.01,0.01)) +
      theme_classic() +
      labs(x=NULL,y=NULL) +
      theme(axis.line = element_blank(), 
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            strip.text = element_blank(),
            legend.position = "none",
            panel.background = element_rect(fill="transparent", colour=NA),
            plot.background=element_rect(fill="transparent", colour=NA),
            strip.background=element_rect(fill="transparent", colour=NA))
    g <- (gd|gt)  + plot_layout(widths = c(4, 3))
    ggsave(paste(path1, "Intermediate files/Manual analysis/Visualized genetic loci/partial/Genetic loci_", type, "_manual analysis_", formatC(k, width = 3, flag=0), ".pdf", sep=""), g, width=15.6, height=1.2*length(pack[[k]])+0.6, units="cm", dpi = 600)
    setTxtProgressBar(pb, k)
  }
  pdf.list <- list.files(paste(path1, "Intermediate files/Manual analysis/Visualized genetic loci/partial/", sep="")) %>>% 
    (.[str_detect(., paste("^priority ", n, "_", sep=""))])
  staple_pdf(input_files = paste(path1, "Intermediate files/Manual analysis/Visualized genetic loci/partial/", pdf.list, sep=""),
             output_filepath = paste(path1, "Intermediate files/Manual analysis/Visualized genetic loci/", type, "for manual analysis.pdf", sep=""),
             overwrite = T)
}
  
##### Combine precursor lists after manual check #####
# Manually check the genetic loci to correct precursor annotations. Save as '...manual analysis_checked.xlsx'

files.xlsx <- list.files(paste(path1, "Intermediate files/Manual analysis/", sep="")) %>>%
  (.[str_detect(.,"Genetic loci.+checked.xlsx")])
prec.checked <- read.xlsx(paste(path1, "Intermediate files/Manual analysis/", files.xlsx[1], sep="")) %>>%
  filter(Check=="p")
for (i in 2:length(files.xlsx)){
  prec.checked <- read.xlsx(paste(path1, "Intermediate files/Manual analysis/", files.xlsx[i], sep="")) %>>%
    filter(Check=="p") %>>%
    (rbind(prec.checked, .))
}
prec.checked <- prec.checked %>>% arrange(Entry) %>>% mutate(f= paste(Sequence, Strain, sep="_"),
                                                             ID = paste(Sequence, Entry, sep="_"))
write.xlsx(prec.checked, paste(path1, "Intermediate files/Manual analysis/4-01_Precursors after manual analysis_all.xlsx", sep=""))


  prec.checked <- read.xlsx(paste(path1, "Intermediate files/Manual analysis/4-01_Precursors after manual analysis_all.xlsx", sep=""))
  duplicates <- lapply(unique(prec.checked$f), function(x){
    temp <- prec.checked %>>% filter(f==x)
    if(dim(temp)[1]>1){
      temp <- temp %>>% mutate(RefSeq = str_detect(temp$Protein.ID, "^WP_"))
      temp %>>% arrange(desc(RefSeq))
    }
  }) %>>% list.rbind
  write.xlsx(duplicates, paste(path1, "Intermediate files/Manual analysis/4-02_Precursors after manual analysis_duplicates.xlsx", sep=""))
  
  ##### Remove duplicates generated during manual analysis #####
  
  duplicates <- read.xlsx(paste(path1, "Intermediate files/Manual analysis/4-02_Precursors after manual analysis_duplicates_checked.xlsx", sep="")) %>>%
    filter(RefSeq==FALSE)
  
  if(dim(duplicates)[1]>0){
  prec.checked.wo.duplicates <- sapply(duplicates$ID, function(x){
    grep(x, prec.checked$ID)
  }) %>>% unlist %>>% (prec.checked[-.,])
  } else {
    prec.checked.wo.duplicates <- prec.checked
  }
  prec.checked.wo.duplicates %>>% arrange(rSAM.SSN.cluster, Strain) %>>% select(-f, -ID) %>>%
  write.xlsx(paste(path1, "Intermediate files/Manual analysis/4-03_Precursors after manual analysis_wo duplicates.xlsx", sep=""))
  
  
##### Generate FASTA for precursor SSN #####
  
  prec.checked.wo.duplicates$FASTA %>>%
    write.table(paste(path1, "Key files/4-04_combined precursor list for SSN.fasta", sep=""),
                quote = F, row.names = F, col.names = F)
  


    
