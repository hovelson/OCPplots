library(ggplot2)
library(reshape2)
library(plyr)
require(gridExtra)
require(scales)
setwd("E:/UMBox/Box Sync/Tomlins/")

tmp <- read.table("IE/geneComb/combined.CNcalls.20141105.txt",header=TRUE,sep="\t")
genes <- unique(tmp$Gene)

# OCP2 genes by genome order
#geneO <- read.table("E:/UMBox/Box Sync/Tomlins/ref/OCP2.20131028.designed.genes_genomeOrder.txt",header=FALSE,sep="\t")
geneO <- read.table("E:/UMBox/Box Sync/Tomlins/ref/OCP.20130724.designed.noTrack.GC.genes_genomeOrder.txt",header=FALSE,sep="\t")
geneO <- geneO[geneO$V1 %in% genes,]
geneO$order <- seq(1,nrow(geneO),1)
#geneO$order <- row.names(geneO)
colnames(geneO) <- c("variable","chr","order")
tmp2 <- merge(tmp,geneO,by.x=c("Gene"),by.y=c("variable"),all.x=TRUE)
tmp2$order2 <- as.numeric(tmp2$order)
tmp2$Sample <- gsub("_","",tmp2$Sample)
#tmp2$Sample <- gsub("PN00","PN",tmp2$Sample)
#tmp2$Sample <- gsub("PN0","PN",tmp2$Sample)

# Order gene factor levels by genome position
tmp2$Gene <- factor(tmp2$Gene,levels=unique(tmp2[order(tmp2$order2),1]))
tmp2$Sample <- factor(tmp2$Sample,unique(as.character(tmp2$Sample)))
tmp2$Chr <- factor(tmp2$chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7",
                                     "chr8","chr9","chr10","chr11","chr12","chr13","chr14",
                                     "chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                     "chr22","chrX"))

# Vertical chr separators & chr labels
geneO$order2 <- as.numeric(geneO$order)
tmp3 <- ddply(geneO, .(chr), summarise, ordmax=max(order2),ordmin=min(order2))
tmp3$order2 <- tmp3$ordmax+0.5
tmp3$clab <- (tmp3$ordmax + tmp3$ordmin)*0.5
lines <- tmp3$order2
lines <- lines[1:22]
clabs <- tmp3$clab
cnames <- tmp3$chr

# error_bar limits
limits <- aes(ymax=log2(CopyNumberRatio + ProbeError),ymin=log2(CopyNumberRatio - ProbeError))
#limits <- aes(ymax=min(log2(CopyNumberRatio + ProbeError),3),ymin=max(log2(CopyNumberRatio - ProbeError),-3))

# grab gene labels
labels <- tmp2$Gene

# grab gene labels
tmp2$Gene2 <- ""
tmp2$colrs <- "black"
tmp2$Cat <- ""
for (i in 1:nrow(tmp2)) {
  t <- log2(tmp2$CopyNumberRatio[i])
  if (t > 0.8 | t < (-0.8)) {
    tmp2$Gene2[i] <- as.character(tmp2$Gene[i])
    if (t > 0.8) {
      tmp2$colrs[i] <- "red"
      tmp2$Cat[i] <- "Gain"
    }
    else {
      tmp2$colrs[i] <- "blue"
      tmp2$Cat[i] <- "Loss"
    }
  }
}
brks <- tmp2$order2
names(brks) <- as.character(tmp2$Gene2)


# Plot function
plotHistFunc <- function(x, panel=NULL, na.rm = TRUE, ...) {
  nm <- unique(x$Sample)
  p <- list()
  for (i in 1:length(nm)) {
    sample <- nm[i]
    s <- as.vector(nm[i])
    tp <- tmp2[tmp2$Sample %in% s,]
    brks <- tp$order2
    names(brks) <- as.character(tp$order2)
    #p <- ggplot(x,aes_string(x = nm[i])) + geom_histogram(alpha = .5,fill = "dodgerblue")
    p[[i]] <- ggplot(tp,aes(x=order2,y=log2(CopyNumberRatio))) + 
      theme_bw() +
      geom_vline(xintercept=lines,linetype="dashed",color="gray30",size=0.25) +
      geom_hline(yintercept=0,linetype="solid",color="gray60",size=0.25) +
      geom_pointrange(limits,size=0.25) +
      geom_point(aes(color=Cat),alpha=0.75)  +
      annotate("text",x=clabs,y=-4,label=cnames,size=1.5) +
      scale_color_manual(breaks=c("","Gain","Loss"),values=c("gray80","red","blue")) +
      scale_y_continuous("log2CN",limits=c(-4,4),breaks=seq(-4,4,1)) +
      scale_x_discrete("",breaks=brks,labels=tp$Gene2) +
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,face="bold",size=7,color="black"),
            panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
            axis.ticks.x=element_blank(),
            plot.margin=unit(c(0,0.5,-0.5,0),"cm"),
            plot.title=element_text(size=10),
            legend.position="none") +
      ggtitle(paste(paste("Copy Number Analysis: ",panel,sep=""),sample,sep="\n"))
    #ggsave(plots,filename=paste("myplot",nm[i],".png",sep=""))
  } 
  return(p)
}

#tlist <- plotHistFunc(tmp2[!(tmp2$Sample %in% c("BL193","BL193A","BL193B","HMG","PRMNOldPool")),],panel="OCPv2") ## execute function

#tlist <- plotHistFunc(tmp2[tmp2$Sample %in% c("PN30","PN60","PN54","PN71","PN75","PN51"),],panel="OCPv2") ## execute function
#tlist <- plotHistFunc(tmp2[tmp2$Sample %in% c("PH3","PH5","PH8","PH13","PH16","PH20","PH22"),],panel="OCPv2") ## execute function
tlist <- plotHistFunc(tmp2,panel="OCPv2") ## execute function
n <- length(tlist)
nRow <- n
isEven <- nRow%%3
do.call("grid.arrange",c(tlist,nrow=nRow))

# Filter high-level, concordant alterations
tsonco <- read.table("lifetech/manuscript/allgenes.TS_onco.20141120.txt",header=TRUE,sep="\t")
tsonco2 <- unique(tsonco[,c("Gene","Category")])  

tmp_hi <- tmp2[tmp2$NumProbes > 3,]
#tmp_hi <- tmp_hi[tmp_hi$NumProbes > 3,]
tmp_hi <- merge(tmp_hi,tsonco2,by=c("Gene"))
tmp_hi <- tmp_hi[log2(tmp_hi$CopyNumberRatio) > 0.8 | log2(tmp_hi$CopyNumberRatio) < -1.0,]
tmp_hi <- tmp_hi[(tmp_hi$Cat %in% c("Loss") & tmp_hi$Category %in% c("Tumor Suppressor") ) | (tmp_hi$Cat %in% c("Gain") & tmp_hi$Category %in% c("Oncogene") ) ,]


# Combined
tmp2$Gene2 <- as.character(tmp2$Gene)
tmp2[!tmp2$Gene %in% c("CDKN2A","PTEN","ARID1B","TP53","NCOR1","FAT1"),]$Gene2 <- "Other"
tmp2$Gene2 <- factor(tmp2$Gene2)
#plottmp <- tmp2[tmp2$Sample %in% c("PN30","PN60","PN54","PN71","PN75","PN51","PN22","PN28"),]
plottmp <- tmp2
brks <- plottmp$order2
names(brks) <- as.character(plottmp$order2)

# AC check
#tmp2$log2CN <- log2(tmp2$CopyNumberRatio)
#acr <- tmp2[(tmp2$log2CN < -1 | tmp2$log2CN > 0.6) & !tmp2$Chr %in% c("chrX"),]
#acr2 <- tmp2[(tmp2$log2CN < -1 | tmp2$log2CN > 0.8) & !tmp2$Chr %in% c("chrX"),]
#dim(acr2)
plottmp$Sample2 <- factor(as.character(plottmp$Sample),levels=sort(unique(as.character(plottmp$Sample))))
plottmp$log2CN <- log2(plottmp$CopyNumberRatio)
plottmp[plottmp$log2CN < -2.0,]$log2CN <- -2.0
plottmp2 <- plottmp[!plottmp$Sample %in% c("BL193A","BL193B","HMG","PRMNOldPool"),]
plottmp2$Cat <- "Other"
plottmp2[which(regexpr("^UMUC",plottmp2$Sample,perl=TRUE) > 0),]$Cat <- "UMUC"
plottmp2$Cat <- factor(plottmp2$Cat,levels=c("UMUC","Other"))

ggplot(plottmp2[plottmp2$NumProbes > 3,],aes(x=reorder(Gene,order),y=log2CN,color=Sample2),shape=21) + 
  theme_bw() +
  #  geom_vline(xintercept=lines,linetype="dashed",color="gray30",size=0.25) +
  geom_hline(yintercept=0,linetype="solid",color="gray60",size=0.25) +
  #geom_pointrange(alpha=0.75,limits,size=0.5) +
  #geom_errorbar(limits,color="black",width=0) +
  geom_point(size=2.5,alpha=0.5)  +
  #annotate("text",x=clabs,y=-4,label=cnames,size=1.5) +
  #scale_color_brewer(palette="Set1") +
  scale_y_continuous("log2CN",limits=c(-2,5.5),breaks=seq(-2,5.0,1.0)) +
  #scale_x_discrete("",breaks=brks,labels=plottmp$Gene) +
  scale_x_discrete() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,face="bold",size=6,color="black"),
        axis.text.y=element_text(size=6,color="black"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=unit(c(0,0.5,-0.5,0),"cm"),
        plot.title=element_text(size=10),
        legend.position="bottom",
        legend.text=element_text(size=8),
        strip.text.y=element_text(angle=0),
        strip.text.x=element_text(size=8),
        strip.background=element_blank(),
        legend.key=element_blank()) +
  #facet_grid(Sample~Chr,scales="free_x",space="free") +
  #facet_grid(Sample~.) +
  facet_grid(Cat~Chr,scales="free_x",space="free_x") +
  guides(col = guide_legend(nrow = 3)) +
  ggtitle("BL Cell Lines CN Analysis: OCPv2\n")

# Ovarian STIC CN plot
pairmx <- as.data.frame(matrix(NA,8,2))
colnames(pairmx) <- c("Sample","Patient") 
pairmx[,1] <- factor(unique(sort(as.character(plottmp$Sample))))
pairmx[,2] <- c(rep(seq(1,4,1)))
pairmx$Patient <- factor(paste("Patient",pairmx$Patient))
plottmp$Sample <- factor(plottmp$Sample,levels=c("OV1","OV2","OV3","OV4","OV5","OV6","OV7","OV8"))
plottmp2 <- merge(plottmp,pairmx,by=c("Sample"))

ggplot(plottmp2[plottmp2$NumProbes > 3,],aes(x=reorder(Gene,order),y=log2(CopyNumberRatio),color=Sample)) + 
  theme_bw() +
  #geom_vline(xintercept=lines,linetype="dashed",color="gray30",size=0.25) +
  geom_hline(yintercept=0,linetype="solid",color="gray60",size=0.25) +
  geom_pointrange(alpha=0.75,limits,size=0.5) +
  #geom_errorbar(limits,color="black",width=0) +
  #geom_point(shape=21,size=3,alpha=0.75)  +
  #annotate("text",x=clabs,y=-4,label=cnames,size=1.5) +
  scale_color_brewer(palette="Set1") +
  scale_y_continuous("log2CN",limits=c(-2,2),breaks=seq(-2,2,0.5)) +
  #scale_x_discrete("",breaks=brks,labels=plottmp$Gene) +
  scale_x_discrete() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,face="bold",size=6,color="black"),
        axis.text.y=element_text(size=6,color="black"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=unit(c(0,0.5,-0.5,0),"cm"),
        plot.title=element_text(size=10),
        legend.position="bottom",
        strip.text.y=element_text(angle=0),
        strip.background=element_blank(),
        legend.key=element_blank()) +
  #facet_grid(Sample~Chr,scales="free_x",space="free") +
  #facet_grid(Sample~.) +
  facet_grid(Patient~.) +
  ggtitle("OV STIC Copy Number Analysis: OCPv2\n")


# Test for CN concordance
plottmp3 <- plottmp2[plottmp2$NumProbes > 3,]
ovconc <- dcast(plottmp3,Gene~Sample,value.var="CopyNumberRatio")

pdf("ovarian_STIC/CN_concordance.20150416.pdf",onefile=TRUE,width=7,height=5)
ggplot(ovconc,aes(OV1,OV5)) + geom_point() +
  scale_y_continuous(limits=c(0,2.5)) +
  scale_x_continuous(limits=c(0,2.5)) +
  geom_smooth(method="lm") +
  theme_bw() +
  ggtitle("Patient 1\n")

ggplot(ovconc,aes(OV2,OV6)) + geom_point() +
  scale_y_continuous(limits=c(0,2.5)) +
  scale_x_continuous(limits=c(0,2.5)) +
  geom_smooth(method="lm") +
  theme_bw()+
  ggtitle("Patient 2\n")

ggplot(ovconc,aes(OV3,OV7)) + geom_point() +
  scale_y_continuous(limits=c(0,2.5)) +
  scale_x_continuous(limits=c(0,2.5)) +
  geom_smooth(method="lm") +
  theme_bw()+
  ggtitle("Patient 3\n")

ggplot(ovconc,aes(OV4,OV8)) + geom_point() +
  scale_y_continuous(limits=c(0,2.5)) +
  scale_x_continuous(limits=c(0,2.5)) +
  geom_smooth(method="lm") +
  theme_bw()+
  ggtitle("Patient 4\n")

dev.off()

# IE samples
ieann <- read.table("IE/geneComb/qry_TISSUE_by_IE.txt",header=TRUE,sep="\t")
ia1 <- ieann[,c("INTERNAL_ACCESSION","TISSUE_COMMENT")]
colnames(ia1) <- c("Sample_","Type")
ia1$Sample <- gsub("ST-14-","",ia1$Sample_)
ia1$Sample <- gsub("IE-0+","IE",ia1$Sample)
ptmp2 <- merge(plottmp,ia1,by=c("Sample"))
#ptmp2 <- ptmp2[!ptmp2$Gene %in% c("MET"),]
ptmp2[ptmp2$CopyNumberRatio < -2,]$CopyNumberRatio <- -2
ggplot(ptmp2[!ptmp2$Sample %in% c("IE23","IE26","IE31","IE18"),]) + 
  theme_bw() +
  geom_vline(xintercept=lines,linetype="dashed",color="gray30",size=0.25) +
  geom_hline(yintercept=0,linetype="solid",color="gray60",size=0.25) +
  #  geom_pointrange(limits,size=0.25) +
  geom_point(aes(x=order2,y=log2(CopyNumberRatio),color=Type),alpha=0.75)  +
  #  geom_boxplot() +
  annotate("text",x=clabs,y=-4.0,label=cnames,size=1.5) +
  #  scale_color_manual(breaks=c("","Gain"),values=c("gray80","purple")) +
  scale_color_brewer(palette="Set1") +
  #  scale_color_manual(breaks=c("ARID1B","CDKN2A","FAT1","NCOR1","PTEN","TP53","Other"),
  #                              values=c("purple","red","blue","darkgreen","gray80","brown","orange")) +
  scale_y_continuous("log2CN",limits=c(-4.0,2),breaks=seq(-4.0,2,0.5)) +
  scale_x_discrete("",breaks=brks,labels=plottmp$Gene) +
  #  scale_x_discrete() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,face="bold",size=6,color="black"),
        axis.text.y=element_text(size=6,color="black"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=unit(c(0,0.5,-0.5,0),"cm"),
        plot.title=element_text(size=10),
        legend.position="bottom",
        strip.text.y=element_text(angle=0),
        strip.background=element_blank()) +
  #facet_grid(Sample~Chr,scales="free_x",space="free") +
  #  facet_grid(.~Chr,scale="free_x",space="free_x") +
  ggtitle("Copy Number Analysis: OCP\n")

###########
# midpt dist matrix
if (FALSE) {
  middist <- read.table("CNtesting/instability/OCP2.20131028.designed.midpoint_dist.20141103.txt",header=TRUE,sep="\t")
  
  # Test out instability metric
  g <- sort(unique(tmp$Gene))
  s <- sort(unique(tmp$Sample))
  
  df <- as.data.frame(matrix(NA,length(s),2))
  
  csum <- ddply(tmp,.(Sample),summarize,cns = sum(CopyNumberRatio))
  for (i in 1:length(s)) {
    s <- 0
    for (j in 1:length(g)) {
      tg <- c(g[j])
      td <- middist[,colnames(middist) %in% tg]
      tc <- 
        s <- s + 
        
        
    }
  }
  
}
