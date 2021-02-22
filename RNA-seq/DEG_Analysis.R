#Code is used for differential gene expression analysis of RNA-seq data for COVID-19 2020 project
#Code developed by Hardip Patel and Alexander Capraro, August 2020

library(tidyverse)
library(edgeR)
library(plotly)
library(gprofiler2)
library(gplots)
library(RUVSeq)
library(Hmisc)
library(pvclust)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(viridis)
library(EnhancedVolcano)
library(rmulti)
library(rstatix)
library(dplyr)
library(ggpubr)

#### Read input files ####
sampleinfo <- read.table("metadata.table", header = T, sep = "\t")
##dashes replaced with dots and start digits removed
sampleinfo <- sampleinfo %>% mutate(sname = gsub("-",".",filebase))
sampleinfo <- sampleinfo %>% dplyr::mutate(individualid = case_when(is.na(individualid) ~ "NA", TRUE ~ as.character(individualid)))
sampleinfo <- sampleinfo %>% unite(col=sid, treatment,agegroup,tissue,individualid,experiment,technicalrepid, remove = F) 

countsummary <- read.table("covseqcounts.txt.summary", header = T, sep = "\t", row.names = 1)
colnames(countsummary) <- gsub(".sorted.bam", "", colnames(countsummary))
colnames(countsummary) <- gsub("^B2.", "", colnames(countsummary))
colnames(countsummary) <- gsub(".xo$", "", colnames(countsummary))
colnames(countsummary) <- gsub("^X", "", colnames(countsummary))
countsummary <- countsummary[rowSums(countsummary)>0,]
countsummary <- data.frame(t(countsummary))
countsummary$sname <- row.names(countsummary)
sampleinfo <- left_join(sampleinfo,countsummary)

##human genome mapped is assigned - virus mapped
sampleinfo <- mutate(sampleinfo, hmapped = Assigned - virusmapped)

##figure for proportion of reads (un)assigned to the virus and the human genes
##sampleinfo %>% dplyr::select(sname,virusmapped, hmapped, Unassigned_Unmapped, Unassigned_NoFeatures) %>% gather("feature", "counts", -sname) %>% mutate(sid = as.factor(sid), feature = as.factor(feature)) %>% left_join(sampleinfo[,c("sid","readcounts")], by = "sid") %>% mutate(proportion = counts/readcounts) %>% dplyr::select(sid,feature,proportion) %>% plot_ly(y=~sid, x=~proportion, color = ~feature, type="bar") %>% layout(barmode = "stack")
csinfo <- select(sampleinfo,sid,readcounts,virusmapped,virus.pos,virus.neg,Assigned,Unassigned_Unmapped,Unassigned_NoFeatures,hmapped) %>% 
  dplyr::group_by(sid) %>% 
  dplyr::summarise_all(sum)

csinfo <- csinfo %>% separate(sid, into = c("treatment","agegroup","tissue","individual","batch","replicate"), remove = F)

select(csinfo, sid, virusmapped, hmapped, Unassigned_Unmapped, Unassigned_NoFeatures) %>% gather("feature", "counts", -sid) %>% 
  mutate(sid = as.factor(sid), feature = as.factor(feature)) %>% 
  left_join(csinfo[,c("sid","readcounts")], by = "sid") %>% 
  mutate(proportion = counts/readcounts) %>% 
  dplyr::select(sid,feature,proportion) %>% 
  plot_ly(y=~sid, x=~proportion, color = ~feature, type="bar") %>% 
  layout(barmode = "stack")


##read counts file
counts <- read.table("covseqcounts.txt", header = T, sep = "\t")
colnames(counts) <- gsub(".sorted.bam$", "", colnames(counts), perl = T)
colnames(counts) <- gsub("^X", "", colnames(counts), perl = T)
colnames(counts) <- gsub("^B2.", "", colnames(counts), perl = T)
colnames(counts) <- gsub(".xo$", "", colnames(counts), perl = T)

##remove columns for one of the samples
cols2retain <- c(colnames(counts)[1:6],sampleinfo$sname)
counts <- counts[,cols2retain]

##merge lane counts
counts <- counts %>% 
  dplyr::select(!(Chr:Length)) %>% 
  gather(sname,value,-Geneid) %>%
  dplyr::left_join(sampleinfo[,c("sname","sid")]) %>%
  dplyr::select(-sname) %>%
  dplyr::group_by(sid,Geneid) %>%
  dplyr::summarise_all(sum) %>%
  spread(sid,value)

lcounts <- as.matrix(select(counts, -Geneid))
row.names(lcounts) <- counts$Geneid
##arrange csinfo table to match the order of colnames in lcounts, this will ensure that you can use the table as is
csinfo <- csinfo[match(colnames(lcounts), csinfo$sid),]

#### Differential gene expression analysis GLM ####

performdiffx <- function(rawcounts, sampleinfo, topthreshold = 0.995, rmvirus = T, plotlabel = "all") {
  
  ##hard filter to retain genes expressed in topthreshold% in any one sample
  seltop <- function (x,threshold=0.99){slist<-cumsum(x[order(x,decreasing=TRUE)])/sum(x)<=threshold; bottom<-labels(slist)[slist==FALSE]; x[bottom]<-NA; return(x)}
  filt <- apply(rawcounts, 2, seltop, threshold=0.995)
  filt <- filt[rowSums(filt,na.rm=TRUE)>0,]
  fcounts <- rawcounts[row.names(filt),]
  
  ##remove virus genes if required
  if (rmvirus) {
    virusgenes <- paste("GU280_gp", str_pad(1:11,2,pad=0),sep = "")
    fcounts <- fcounts[! row.names(fcounts) %in% virusgenes,]
  }
  
  dge <- DGEList(counts = as.matrix(fcounts), lib.size = sampleinfo$readcounts, group = sampleinfo$treatment)
  dge <- calcNormFactors(dge)
  tpm <- cpm.DGEList(dge, log = T)
  
  design <- model.matrix(~sampleinfo$individual+sampleinfo$treatment)
  y <- estimateDisp(dge,design)
  fit <- glmFit(y, design)
  qlf <- glmLRT(fit)
  diffex <- topTags(qlf, n = nrow(fcounts))
  dx <- diffex$table
  dx$gid <- row.names(dx)
  
  dxgenes <- tpm[dx %>% filter(FDR <= 0.05) %>% .$gid,]
  write.table(dx, file = paste("diffex_glm",plotlabel, topthreshold, "txt", sep ="."), col.names = T, row.names = T, quote = F, sep = "\t")
  heatmap.2(as.matrix(tpm[dx %>% filter(FDR <= 0.05) %>% .$gid,]), trace = 'none', density.info = 'none', col = viridis)
} 

sampleinfo <- csinfo
rawcounts <- lcounts
performdiffx(rawcounts=rawcounts, sampleinfo=sampleinfo, topthreshold = 0.995, rmvirus = T, plotlabel = "all")

sampleinfo <- csinfo %>% filter(! (individual=="BC" & batch=="B1"))
rawcounts <- lcounts[,sampleinfo$sid]
performdiffx(rawcounts=rawcounts, sampleinfo=sampleinfo, topthreshold = 0.995, rmvirus = T, plotlabel = "allminusBC")

sampleinfo <- csinfo %>% filter(! (individual=="BC" & batch=="B1") & tissue == "B")
rawcounts <- lcounts[,sampleinfo$sid]
performdiffx(rawcounts=rawcounts, sampleinfo=sampleinfo, topthreshold = 0.995, rmvirus = T, plotlabel = "minusBC.tissueB")

sampleinfo <- csinfo %>% filter(! (individual=="BC" & batch=="B1") & tissue == "N")
rawcounts <- lcounts[,sampleinfo$sid]
performdiffx(rawcounts=rawcounts, sampleinfo=sampleinfo, topthreshold = 0.995, rmvirus = T, plotlabel = "minusBC.tissueN")

sampleinfo <- csinfo %>% filter(! (individual=="BC" & batch=="B1") & agegroup == "OA")
rawcounts <- lcounts[,sampleinfo$sid]
performdiffx(rawcounts=rawcounts, sampleinfo=sampleinfo, topthreshold = 0.995, rmvirus = T, plotlabel = "minusBC.OA")

sampleinfo <- csinfo %>% filter(! (individual=="BC" & batch=="B1") & agegroup == "YA")
rawcounts <- lcounts[,sampleinfo$sid]
performdiffx(rawcounts=rawcounts, sampleinfo=sampleinfo, topthreshold = 0.995, rmvirus = T, plotlabel = "minusBC.YA")

sampleinfo <- csinfo %>% filter(! (individual=="BC" & batch=="B1") & agegroup == "CH" & tissue == "N")
rawcounts <- lcounts[,sampleinfo$sid]
performdiffx(rawcounts=rawcounts, sampleinfo=sampleinfo, topthreshold = 0.995, rmvirus = T, plotlabel = "minusBC.CH.tissueN")

sampleinfo <- csinfo %>% filter(! (individual=="BC" & batch=="B1") & agegroup != "CH")
rawcounts <- lcounts[,sampleinfo$sid]
performdiffx(rawcounts=rawcounts, sampleinfo=sampleinfo, topthreshold = 0.995, rmvirus = T, plotlabel = "minusBC.Adults")




