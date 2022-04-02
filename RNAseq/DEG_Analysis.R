# Code developed for QC and differential expression analysis of RNA-seq data (3 batches) for COVID-19 Age Study
# Code developed by Hardip Patel & Alexander Capraro

library(tidyverse)
library(edgeR)
library(plotly)
library(gprofiler2)
library(RUVSeq)
library(Hmisc)
library(pvclust)
library(ggpubr)
library(ggrepel)
library(gplots)
library(viridis)
library(EnhancedVolcano)
options(scipen = 100)

# Load files & functions ####
sinfo <- read_delim("sinfo.txt", delim="\t")
cinfo <- read_delim("cinfo.txt", delim="\t")

diffx <- function(rawcounts, sdetails, rmvirus, plotlabel) {
  cmat <- as.matrix(select(rawcounts, -Geneid))
  row.names(cmat) <- rawcounts$Geneid
  filt <- cmat[rowSums(cmat,na.rm=TRUE)>0,]
  fcounts <- cmat[row.names(filt),]
  
  ##remove virus genes if required
  if (rmvirus) {
    virusgenes <- paste("GU280_gp", str_pad(1:11,2,pad=0),sep = "")
    fcounts <- fcounts[! row.names(fcounts) %in% virusgenes,]
  }
  
  design <- model.matrix(~sdetails$individualid + sdetails$treatment)
  x <- DGEList(counts=as.matrix(fcounts), lib.size = sdetails$readcount, group = sdetails$treatment)
  x <- calcNormFactors(x, method = "upperquartile")
  x <- estimateGLMCommonDisp(x, design)
  x <- estimateGLMTagwiseDisp(x, design)
  
  fit <- glmFit(x, design)
  res <- residuals(fit, type="deviance")
  normset <- RUVr(fcounts, row.names(fcounts), k=1, res)
  normcounts <- normset$normalizedCounts
  
  design <- model.matrix(~sdetails$individualid + sdetails$treatment)
  dge <- DGEList(counts=normcounts, lib.size = sdetails$readcount, group = sdetails$treatment)
  y <- calcNormFactors(dge)
  #tpm <- cpm.DGEList(y, log=T)
  y <- estimateDisp(y,design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit)
  diffex <- topTags(lrt, n = nrow(fcounts))
  dx <- diffex$table %>%
    mutate(gid = row.names(diffex$table)) %>%
    select(gid, logFC,logCPM,LR,PValue,FDR)
  dxgenes <- dx[dx %>% filter(FDR <= 0.05, logFC > 0.585 | logFC < -0.585) %>% .$gid,]
  write_delim(dx, file = paste("./DGE/diffex_glm",plotlabel, "txt", sep ="."), delim="\t")
} 

# Counts Plots ####

#Calculate % viral reads

viral_percent <- sinfo %>%
  mutate(percentage = (virusmapped/readcount)*100) %>%
  select(sname,treatment,individualid,technicalrepid,sex,agegroup,tissue,percentage) %>%
  filter(treatment != "Mock")

write.table(viral_percent,file="./viral_percent.txt",sep="\t",row.names=F,quote=F)


figure1 <- sinfo %>% 
  select(sname,Unassigned_Unmapped,virusmapped,humanmapped) %>%
  pivot_longer(!sname, names_to = "category", values_to = "counts") %>%
  group_by(sname) %>% mutate(proportions = counts/sum(counts)) %>%
  plot_ly(x = ~sname, y = ~counts, type = 'bar',color = ~category,text = ~counts, hoverinfo = 'text') %>% 
  layout(yaxis = list(title = 'Read counts'),barmode = "stack")


figure2 <- sinfo %>% 
  select(sname,Unassigned_Unmapped,virusmapped,humanmapped) %>%
  pivot_longer(!sname, names_to = "category", values_to = "counts") %>%
  group_by(sname) %>% mutate(proportions = counts/sum(counts)) %>%
  plot_ly(x = ~sname, y = ~proportions*100, type = 'bar',color = ~category,text = ~round(proportions*100, digits = 1), hoverinfo = 'text') %>% 
  layout(yaxis = list(title = 'Percentage (%)'),barmode = "stack")


figure1

figure2

# MDS Plots ####

drawmds_all <- function(rawcounts, sdetails, rmvirus = T, plotlabel = "all") {
  cmat <- as.matrix(select(rawcounts, -Geneid))
  row.names(cmat) <- rawcounts$Geneid
  filt <- apply(cmat, 2, seltop, threshold=topthreshold)
  filt <- filt[rowSums(filt,na.rm=TRUE)>0,]
  fcounts <- cmat[row.names(filt),]
  
  ##remove virus genes if required
  if (rmvirus) {
    virusgenes <- paste("GU280_gp", str_pad(1:11,2,pad=0),sep = "")
    fcounts <- fcounts[! row.names(fcounts) %in% virusgenes,]
  }
  
  design <- model.matrix(~ sdetails$individualid + sdetails$treatment)
  y <- DGEList(counts=as.matrix(fcounts), lib.size = sdetails$readcount, group = sdetails$treatment)
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  
  fit <- glmFit(y, design)
  res <- residuals(fit, type="deviance")
  normset <- RUVr(fcounts, row.names(fcounts), k=2, res)
  z <- plotMDS(log(normcounts),gene.selection = "common", top = 20000,plot = F)$cmdscale.out
  mdsval <- as.tibble(z) %>% mutate(sname = row.names(z))
  colnames(mdsval) <- c("leading logFC dim1", "leading logFC dim2", "sname")
  mdsval <- left_join(mdsval, sdetails, by = "sname")
  mdsval %>% ggplot() + 
    geom_point(aes(x=`leading logFC dim1`, y = `leading logFC dim2`, shape = agegroup, colour = treatment), size = 3) + 
    geom_label_repel(aes(x=`leading logFC dim1`, y = `leading logFC dim2`,label = individualid, fill=experiment), alpha = 0.5, size = 2, segment.color = 'black') +
    theme_classic()
} 

sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1") & !(individualid == "HK" & treatment == "Mock") & !(individualid == "CH" & treatment == "Virus") & !(individualid == "RG" & treatment == "Virus"), tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
drawmds_all(rawcounts = rawcounts, sdetails = sdetails)


# DGE Mock vs. Virus 0.2 ####


#Age-independent Nasal DGE (0.2 MOI)
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1") &
                     !(treatment == "Virus6") &
                     !(treatment == "Nafa") &
                     !(treatment == "5T"), 
                   tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "all_minusOutliers")

#Child Nasal DGE

sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1") & 
                     !(treatment == "Virus6") &
                     !(treatment == "Nafa")&
                     !(treatment == "5T"), 
                   tissue == "N" & agegroup == "CH")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "CH_HNE")

#Young Adult DGE

sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1") & 
                     !(treatment == "Virus6") &
                     !(treatment == "Nafa")&
                     !(treatment == "5T"), 
                   tissue == "N" & agegroup == "YA")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "YA")

#Older Adult DGE

sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1") & 
                     !(treatment == "Virus6") &
                     !(treatment == "Nafa")&
                     !(treatment == "5T"), 
                   tissue == "N" & agegroup == "OA")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "OA")

#Bronchial Child DGE

sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1") & 
                     !(treatment == "Virus6") &
                     !(treatment == "Nafa")&
                     !(treatment == "5T"), 
                   tissue == "B" & agegroup == "CH")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "CH_HBE")

# DGE Mock vs. Virus 0.6 ####

#All

sdetails <- filter(sinfo, !(treatment == "Virus") &
                     !(treatment == "Nafa") &
                     !(treatment == "B"), 
                   individualid == "BM" | individualid == "HK" | individualid == "HL" |
                     individualid == "021CsH" | individualid == "MB" | individualid == "SLW" |
                     individualid == "AT" | individualid == "PT" | individualid == "RG") #Only pairwise comparison

rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "All_0.6vMock")


#Child

sdetails <- filter(sinfo, !(treatment == "Virus") &
                     !(treatment == "Nafa") &
                     !(treatment == "B"), 
                   individualid == "BM" | individualid == "HK" | individualid == "HL") #Only pairwise comparison

rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "CH_0.6vMock")

#Young Adult

sdetails <- filter(sinfo, !(treatment == "Virus") &
                     !(treatment == "Nafa") &
                     !(treatment == "B"), 
                   individualid == "021CsH" | individualid == "MB" | individualid == "SLW") 
rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "YA_0.6vMock")

#Older Adult

sdetails <- filter(sinfo, !(treatment == "Virus") &
                     !(treatment == "Nafa") &
                     !(treatment == "B"), 
                   individualid == "AT" | individualid == "PT" | individualid == "RG" | individualid == "RB" & treatment == "Mock") 

rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "OA_0.6vMock")

# DGE Virus 0.2 vs. Virus 0.6 ####

#All

sdetails <- filter(sinfo, !(treatment == "Mock") &
                     !(treatment == "Nafa") &
                     !(treatment == "B"), 
                   individualid == "BM" | individualid == "HK" | individualid == "HL" |
                     individualid == "021CsH" | individualid == "MB" | individualid == "SLW" |
                     individualid == "AT" | individualid == "PT" | individualid == "RG") #Only pairwise comparison

rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "All_0.6v0.2")


#Child

sdetails <- filter(sinfo, !(treatment == "Mock") &
                     !(treatment == "Nafa") &
                     !(treatment == "B"), 
                   individualid == "BM" | individualid == "HK" | individualid == "HL") #Only pairwise comparison

rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "CH_0.6v0.2")

#Young Adult

sdetails <- filter(sinfo, !(treatment == "Mock") &
                     !(treatment == "Nafa") &
                     !(treatment == "B"), 
                   individualid == "021CsH" | individualid == "MB" | individualid == "SLW") 

rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "YA_0.6v0.2")

#Older Adult

sdetails <- filter(sinfo, !(treatment == "Mock") &
                     !(treatment == "Nafa") &
                     !(treatment == "B"), 
                   individualid == "AT" | individualid == "PT" | individualid == "RG" | 
                     individualid == "RB" & treatment == "Virus" | 
                     individualid == "WB" & treatment == "Virus" & technicalrepid == "t1")

rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "OA_0.6v0.2")

# DGE Mock vs. Nafa ####

#All

sdetails <- filter(sinfo, !(treatment == "Virus") &
                     !(treatment == "Virus6") &
                     !(treatment == "B"), 
                   individualid == "BM" | individualid == "HK" | individualid == "LT" |
                     individualid == "021CsH" | individualid == "MB" | individualid == "SLW" |
                     individualid == "RB" | individualid == "PT" | individualid == "RG") #Only pairwise comparison

rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "All_NafavMock")


#Child

sdetails <- filter(sinfo, !(treatment == "Virus") &
                     !(treatment == "Virus6") &
                     !(treatment == "B"), 
                   individualid == "BM" | individualid == "HK" | individualid == "LT") #Only pairwise comparison

rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "CH_NafavMock")

#Young Adult

sdetails <- filter(sinfo, !(treatment == "Virus") &
                     !(treatment == "Virus6") &
                     !(treatment == "B"), 
                   individualid == "021CsH" | individualid == "MB" | individualid == "SLW") 

rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "YA_NafavMock")

#Older Adult

sdetails <- filter(sinfo, !(treatment == "Virus") &
                     !(treatment == "Virus6") &
                     !(treatment == "B"), 
                   individualid == "RB" | individualid == "PT" | individualid == "RG")

rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "OA_NafavMock")

# DGE Virus 0.2 vs. Nafa ####

#All

sdetails <- filter(sinfo, !(treatment == "Mock") &
                     !(treatment == "Virus6") &
                     !(treatment == "B"), 
                   individualid == "BM" | individualid == "HK" | individualid == "LT" |
                     individualid == "021CsH" | individualid == "MB" | individualid == "SLW" |
                     individualid == "RB" | individualid == "PT" | individualid == "RG") #Only pairwise comparison

rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "All_Nafav0.2")


#Child

sdetails <- filter(sinfo, !(treatment == "Mock") &
                     !(treatment == "Virus6") &
                     !(treatment == "B"), 
                   individualid == "BM" | individualid == "HK" | individualid == "LT") #Only pairwise comparison

rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "CH_Nafav0.2")

#Young Adult

sdetails <- filter(sinfo, !(treatment == "Mock") &
                     !(treatment == "Virus6") &
                     !(treatment == "B"), 
                   individualid == "021CsH" | individualid == "MB" | individualid == "SLW") 

rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "YA_Nafav0.2")

#Older Adult

sdetails <- filter(sinfo, !(treatment == "Mock") &
                     !(treatment == "Virus6") &
                     !(treatment == "B"), 
                   individualid == "RB" | individualid == "PT" | individualid == "RG" |
                     individualid == "WB" & technicalrepid == "t1")

rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "OA_Nafav0.2")

# DGE Age Comparison ####

diffx_age <- function(rawcounts, sdetails, rmvirus, plotlabel) {
  cmat <- as.matrix(select(rawcounts, -Geneid))
  row.names(cmat) <- rawcounts$Geneid
  filt <- cmat[rowSums(cmat,na.rm=TRUE)>0,]
  fcounts <- cmat[row.names(filt),]
  
  ##remove virus genes if required
  if (rmvirus) {
    virusgenes <- paste("GU280_gp", str_pad(1:11,2,pad=0),sep = "")
    fcounts <- fcounts[! row.names(fcounts) %in% virusgenes,]
  }
  
  design <- model.matrix(~groups$agegroup)
  x <- DGEList(counts=as.matrix(fcounts), lib.size = sdetails$readcount, group = groups$agegroup)
  x <- calcNormFactors(x, method = "upperquartile")
  x <- estimateGLMCommonDisp(x, design)
  x <- estimateGLMTagwiseDisp(x, design)
  
  fit <- glmFit(x, design)
  res <- residuals(fit, type="deviance")
  normset <- RUVr(fcounts, row.names(fcounts), k=1, res)
  normcounts <- normset$normalizedCounts
  
  dge <- DGEList(counts=normcounts, lib.size = sdetails$readcount, group = groups$agegroup)
  y <- calcNormFactors(dge)
  y <- estimateDisp(y,design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit)
  diffex <- topTags(lrt, n = nrow(fcounts))
  dx <- diffex$table %>%
    mutate(gid = row.names(diffex$table)) %>%
    select(gid, logFC,logCPM,LR,PValue,FDR) %>%
    filter(FDR <= 0.05 & logFC < -0.585 | logFC > 0.585)
  write_delim(dx, file = paste("D:/OneDrive - UNSW/miCF/SARS-CoV2_Study/RNAseq_Batch3/DGE/diffex_glm",plotlabel,"txt", sep ="."), delim="\t")
}

#Child vs. Adult - Mock
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1"),
                   treatment == "Mock" & tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
groups <- sdetails %>%
  select(., agegroup) %>%
  mutate(agegroup = ifelse(grepl("CH",.$agegroup),"CH","Adult"))

diffx_age(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "ChildvAdult_Mock")

#Young Adult vs. All - Mock
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1"),
                   treatment == "Mock" & tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
groups <- sdetails %>%
  select(., agegroup) %>%
  mutate(agegroup = ifelse(grepl("YA",.$agegroup),"YA","All"))

diffx_age(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "YoungAdultvsAll_Mock")

#Older Adult vs. All - Mock
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1"),
                   treatment == "Mock" & tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
groups <- sdetails %>%
  select(., agegroup) %>%
  mutate(agegroup = ifelse(grepl("OA",.$agegroup),"YA","All"))

diffx_age(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "OlderAdultvsAll_Mock")


#Child vs. Adult - Virus
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1"),
                   treatment == "Virus" & tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
groups <- sdetails %>%
  select(., agegroup) %>%
  mutate(agegroup = ifelse(grepl("CH",.$agegroup),"CH","Adult"))

diffx_age(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "ChildvAdult_Virus")

#Young Adult vs. All - Virus
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1"),
                   treatment == "Virus" & tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
groups <- sdetails %>%
  select(., agegroup) %>%
  mutate(agegroup = ifelse(grepl("YA",.$agegroup),"YA","All"))

diffx_age(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "YoungAdultvsAll_Virus")

#Older Adult vs. All - Virus
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1"),
                   treatment == "Virus" & tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
groups <- sdetails %>%
  select(., agegroup) %>%
  mutate(agegroup = ifelse(grepl("OA",.$agegroup),"YA","All"))

diffx_age(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "OlderAdultvsAll_Virus")

# Mock vs. Virus Volcano plots ####

dx <- read.table(file = "./DGE/diffex_glm.CH_HNE.txt",sep="\t",header=T)

{
  dx$PValue <- as.numeric(as.character(dx$PValue))
  
  keyvals <- ifelse(
  dx$logFC <= -0.585 & dx$FDR < 0.05, '#E8495C',
  ifelse(dx$logFC >= 0.585 & dx$FDR < 0.05, '#38ACE2',
         'grey'))
  names(keyvals)[keyvals == '#38ACE2'] <- 'Upregulated'
  names(keyvals)[keyvals == 'grey'] <- 'Not significant'
  names(keyvals)[keyvals == '#E8495C'] <- 'Downregulated'
  
  EnhancedVolcano(dx,
                  lab = row.names(dx),
                  selectLab = c('SPIKE'),
                  x='logFC',
                  y='FDR',
                  pCutoff = 0.05,
                  FCcutoff = 0.585,
                  colCustom = keyvals,
                  colAlpha = 1, 
                  legendPosition = 'none',
                  xlim=c(-10,10), ylim=c(0,200)) #export pdf as 5.5x6
}


dx <- read.table(file = "./DGE/diffex_glm.YA.txt",sep="\t",header=T)

{
  dx$PValue <- as.numeric(as.character(dx$PValue))
  
  keyvals <- ifelse(
    dx$logFC <= -0.585 & dx$FDR < 0.05, '#E8495C',
    ifelse(dx$logFC >= 0.585 & dx$FDR < 0.05, '#38ACE2',
           'grey'))
  names(keyvals)[keyvals == '#38ACE2'] <- 'Upregulated'
  names(keyvals)[keyvals == 'grey'] <- 'Not significant'
  names(keyvals)[keyvals == '#E8495C'] <- 'Downregulated'
  
  EnhancedVolcano(dx,
                  lab = row.names(dx),
                  selectLab = c('SPIKE'),
                  x='logFC',
                  y='FDR',
                  pCutoff = 0.05,
                  FCcutoff = 0.585,
                  colCustom = keyvals,
                  colAlpha = 1, 
                  legendPosition = 'none',
                  xlim=c(-10,10), ylim=c(0,55)) #export pdf as 5.5x6
}

dx <- read.table(file = "./DGE/diffex_glm.OA.txt",sep="\t",header=T)

{
  dx$PValue <- as.numeric(as.character(dx$PValue))
  
  keyvals <- ifelse(
    dx$logFC <= -0.585 & dx$FDR < 0.05, '#E8495C',
    ifelse(dx$logFC >= 0.585 & dx$FDR < 0.05, '#38ACE2',
           'grey'))
  names(keyvals)[keyvals == '#38ACE2'] <- 'Upregulated'
  names(keyvals)[keyvals == 'grey'] <- 'Not significant'
  names(keyvals)[keyvals == '#E8495C'] <- 'Downregulated'
  
  EnhancedVolcano(dx,
                  lab = row.names(dx),
                  selectLab = c('SPIKE'),
                  x='logFC',
                  y='FDR',
                  pCutoff = 0.05,
                  FCcutoff = 0.585,
                  colCustom = keyvals,
                  colAlpha = 1, 
                  legendPosition = 'none',
                  xlim=c(-10,10), ylim=c(0,25)) #export pdf as 5.5x6
}
                
# DGE Age Comparison - Mock & Nafa ####

diffx_age <- function(rawcounts, sdetails, rmvirus, plotlabel) {
  cmat <- as.matrix(select(rawcounts, -Geneid))
  row.names(cmat) <- rawcounts$Geneid
  filt <- cmat[rowSums(cmat,na.rm=TRUE)>0,]
  fcounts <- cmat[row.names(filt),]
  
  ##remove virus genes if required
  if (rmvirus) {
    virusgenes <- paste("GU280_gp", str_pad(1:11,2,pad=0),sep = "")
    fcounts <- fcounts[! row.names(fcounts) %in% virusgenes,]
  }
  
  design <- model.matrix(~sdetails$agegroup)
  x <- DGEList(counts=as.matrix(fcounts), lib.size = sdetails$readcount, group = sdetails$agegroup)
  x <- calcNormFactors(x, method = "upperquartile")
  x <- estimateGLMCommonDisp(x, design)
  x <- estimateGLMTagwiseDisp(x, design)
  
  fit <- glmFit(x, design)
  res <- residuals(fit, type="deviance")
  normset <- RUVr(fcounts, row.names(fcounts), k=1, res)
  normcounts <- normset$normalizedCounts
  
  dge <- DGEList(counts=normcounts, lib.size = sdetails$readcount, group = sdetails$agegroup)
  y <- calcNormFactors(dge)
  y <- estimateDisp(y,design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit)
  diffex <- topTags(lrt, n = nrow(fcounts))
  dx <- diffex$table %>%
    mutate(gid = row.names(diffex$table)) %>%
    select(gid, logFC,logCPM,LR,PValue,FDR) %>%
    filter(FDR <= 0.05) %>%
    filter(logFC < -0.585 | logFC > 0.585)
  write_delim(dx, file = paste("./DGE/diffex_glm",plotlabel,"txt", sep ="."), delim="\t")
}

#### Mock

#Child vs. Young Adult - Mock
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1") &
                     !(agegroup == "OA"),
                   treatment == "Mock" & tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx_age(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "CHvsYA_Mock")

#Child vs. Older Adult - Mock
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1") &
                     !(agegroup == "YA"),
                   treatment == "Mock" & tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx_age(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "CHvsOA_Mock")

#Young Adult vs. Older Adult - Mock
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1") &
                     !(agegroup == "CH"),
                   treatment == "Mock" & tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx_age(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "YAvsOA_Mock")

#### Virus 0.2 MOI

#Child vs. Young Adult - Virus 0.2
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1") &
                     !(agegroup == "OA"),
                   treatment == "Virus" & tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx_age(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "CHvsYA_V2")

#Child vs. Older Adult - Virus 0.2
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1") &
                     !(agegroup == "YA"),
                   treatment == "Virus" & tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx_age(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "CHvsOA_V2")

#Young Adult vs. Older Adult - Virus 0.2
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1") &
                     !(agegroup == "CH"),
                   treatment == "Virus" & tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx_age(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "YAvsOA_V2")

#### Nafamostat

#Child vs. Young Adult - Nafa
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1") &
                     !(agegroup == "OA"),
                   treatment == "Nafa" & tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx_age(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "CHvsYA_Nafa")

#Child vs. Older Adult - Nafa
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1") &
                     !(agegroup == "YA"),
                   treatment == "Nafa" & tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx_age(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "CHvsOA_Nafa")

#Young Adult vs. Older Adult - Nafa
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1") &
                     !(agegroup == "CH"),
                   treatment == "Nafa" & tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
diffx_age(rawcounts=rawcounts, sdetails=sdetails, rmvirus = T, plotlabel = "YAvsOA_Nafa")
