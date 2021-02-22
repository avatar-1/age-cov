#Code is used for age-specific differential gene expression analysis of RNA-seq data for COVID-19 2020 project
#Also contains code for the analysis of age-associated genes gathered from Chow et al 2020
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
library(dplyr)
library(ggpubr)
library(rstatix)

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

###remove JM sample as it is to be excluded
sampleinfo <- filter(sampleinfo, individualid != "JM")

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


#### Pairwise age group vs. the rest MOCK ####

seltop <- function (x,threshold=0.99){slist<-cumsum(x[order(x,decreasing=TRUE)])/sum(x)<=threshold; bottom<-labels(slist)[slist==FALSE]; x[bottom]<-NA; return(x)}
sampleinfo <- csinfo %>% filter(! (individual=="BC" & batch=="B1") & treatment == "Mock" & tissue == "N")
rawcounts <- lcounts[,sampleinfo$sid]
filt <- apply(rawcounts, 2, seltop, threshold=0.995)
filt <- filt[rowSums(filt,na.rm=TRUE)>0,]
fcounts <- rawcounts[row.names(filt),]

# Child vs. All MOCK #
#create DGE list and normalise
dge <- DGEList(counts = as.matrix(fcounts), lib.size = sampleinfo$readcounts, group = c("CH","CH","CH","CH","CH","Adult","Adult","Adult","Adult","Adult","Adult","Adult","Adult","Adult"))
dge <- calcNormFactors(dge)
tpm <- cpm.DGEList(dge, log = T)

#create experimental design and glm test
design <- model.matrix(~0+group, data = dge$samples)
y <- estimateDisp(dge,design)
fit <- glmFit(y, design)

qlf.CHvsALL <- glmLRT(fit, contrast = c(-1,1))
topTags(qlf.CHvsALL)
diffex.CHvsALL <- topTags(qlf.CHvsALL, n = nrow(fcounts))
dx.CHvsALL <- diffex.CHvsALL$table
dx.CHvsALL <- dx.CHvsALL[dx.CHvsALL$FDR < 0.05, ]
write.table(dx.CHvsALL, file = "diffex.CHvsALL1.txt", sep = "\t", quote = F)

# YA vs. All MOCK #
#create DGE list and normalise
dge <- DGEList(counts = as.matrix(fcounts), lib.size = sampleinfo$readcounts, group = c("All","All","All","All","All","All","All","All","All","All","YA","YA","YA","YA"))
dge <- calcNormFactors(dge)
tpm <- cpm.DGEList(dge, log = T)

#create experimental design and glm test
design <- model.matrix(~0+group, data = dge$samples)
y <- estimateDisp(dge,design)
fit <- glmFit(y, design)

qlf.YAvsALL <- glmLRT(fit, contrast = c(-1,1))
topTags(qlf.YAvsALL)
diffex.YAvsALL <- topTags(qlf.YAvsALL, n = nrow(fcounts))
dx.YAvsALL <- diffex.YAvsALL$table
dx.YAvsALL <- dx.YAvsALL[dx.YAvsALL$FDR < 0.05, ]
write.table(dx.YAvsALL, file = "diffex.YAvsALL.txt", sep = "\t", quote = F)

# OA vs. All MOCK #
#create DGE list and normalise
dge <- DGEList(counts = as.matrix(fcounts), lib.size = sampleinfo$readcounts, group = c("All","All","All","All","All","OA","OA","OA","OA","OA","All","All","All","All"))
dge <- calcNormFactors(dge)
tpm <- cpm.DGEList(dge, log = T)

#create experimental design and glm test
design <- model.matrix(~0+group, data = dge$samples)
y <- estimateDisp(dge,design)
fit <- glmFit(y, design)

qlf.OAvsALL <- glmLRT(fit, contrast = c(-1,1))
topTags(qlf.OAvsALL)
diffex.OAvsALL <- topTags(qlf.OAvsALL, n = nrow(fcounts))
dx.OAvsALL <- diffex.OAvsALL$table
dx.OAvsALL <- dx.OAvsALL[dx.OAvsALL$FDR < 0.05, ]
write.table(dx.OAvsALL, file = "diffex.OAvsALL.txt", sep = "\t", quote = F)

#### Pairwise age group vs. the rest VIRUS ####

seltop <- function (x,threshold=0.99){slist<-cumsum(x[order(x,decreasing=TRUE)])/sum(x)<=threshold; bottom<-labels(slist)[slist==FALSE]; x[bottom]<-NA; return(x)}
sampleinfo <- csinfo %>% filter(! (individual=="BC" & batch=="B1") & treatment == "Virus" & tissue == "N")
rawcounts <- lcounts[,sampleinfo$sid]
filt <- apply(rawcounts, 2, seltop, threshold=0.995)
filt <- filt[rowSums(filt,na.rm=TRUE)>0,]
fcounts <- rawcounts[row.names(filt),]

# CH vs. All Virus #

seltop <- function (x,threshold=0.99){slist<-cumsum(x[order(x,decreasing=TRUE)])/sum(x)<=threshold; bottom<-labels(slist)[slist==FALSE]; x[bottom]<-NA; return(x)}
sampleinfo <- csinfo %>% filter(! (individual=="BC" & batch=="B1") & treatment == "Virus" & tissue == "N")
rawcounts <- lcounts[,sampleinfo$sid]
filt <- apply(rawcounts, 2, seltop, threshold=0.995)
filt <- filt[rowSums(filt,na.rm=TRUE)>0,]
fcounts <- rawcounts[row.names(filt),]

#create DGE list and normalise
dge <- DGEList(counts = as.matrix(fcounts), lib.size = sampleinfo$readcounts, group = c("CH","CH","CH","CH","CH","CH","CH","CH","Adult","Adult","Adult","Adult","Adult","Adult","Adult","Adult","Adult","Adult"))
dge <- calcNormFactors(dge)
tpm <- cpm.DGEList(dge, log = T)

#create experimental design and glm test
design <- model.matrix(~0+group, data = dge$samples)
y <- estimateDisp(dge,design)
fit <- glmFit(y, design)

qlf.CHvsALL <- glmLRT(fit, contrast = c(-1,1))
topTags(qlf.CHvsALL)
diffex.CHvsALL <- topTags(qlf.CHvsALL, n = nrow(fcounts))
dx.CHvsALL <- diffex.CHvsALL$table
dx.CHvsALL <- dx.CHvsALL[dx.CHvsALL$FDR < 0.05, ]
write.table(dx.CHvsALL, file = "diffex.CHvsALL_virus.txt", sep = "\t", quote = F)

# YA vs. All Virus #
#create DGE list and normalise
dge <- DGEList(counts = as.matrix(fcounts), lib.size = sampleinfo$readcounts, group = c("All","All","All","All","All","All","All","All","All","All","All","All","All","All","YA","YA","YA","YA"))
dge <- calcNormFactors(dge)
tpm <- cpm.DGEList(dge, log = T)

#create experimental design and glm test
design <- model.matrix(~0+group, data = dge$samples)
y <- estimateDisp(dge,design)
fit <- glmFit(y, design)

qlf.YAvsALL <- glmLRT(fit, contrast = c(-1,1))
topTags(qlf.YAvsALL)
diffex.YAvsALL <- topTags(qlf.YAvsALL, n = nrow(fcounts))
dx.YAvsALL <- diffex.YAvsALL$table
dx.YAvsALL <- dx.YAvsALL[dx.YAvsALL$FDR < 0.05, ]
write.table(dx.YAvsALL, file = "diffex.YAvsALL_virus.txt", sep = "\t", quote = F)

# OA vs. All Virus #
#create DGE list and normalise
dge <- DGEList(counts = as.matrix(fcounts), lib.size = sampleinfo$readcounts, group = c("All","All","All","All","All","All","All","All","OA","OA","OA","OA","OA","OA","All","All","All","All"))
dge <- calcNormFactors(dge)
tpm <- cpm.DGEList(dge, log = T)

#create experimental design and glm test
design <- model.matrix(~0+group, data = dge$samples)
y <- estimateDisp(dge,design)
fit <- glmFit(y, design)

qlf.OAvsALL <- glmLRT(fit, contrast = c(-1,1))
topTags(qlf.OAvsALL)
diffex.OAvsALL <- topTags(qlf.OAvsALL, n = nrow(fcounts))
dx.OAvsALL <- diffex.OAvsALL$table
dx.OAvsALL <- dx.OAvsALL[dx.OAvsALL$FDR < 0.05, ]
write.table(dx.OAvsALL, file = "diffex.OAvsALL_virus.txt", sep = "\t", quote = F)

#### Age-associated gene correlation Nasal + Bronchial (USED IN PAPER) ####

#Mock
tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm[,c(1:16,18:40,42:44)] #Remove BC_B1
tpm_plots <- tpm[,grepl('Mock', colnames(tpm))] #only keep mock

#lung_age <- read.table(file = "/Users/z3416833/Documents/miCF/SARS-CoV2_Study/RNAseq_Complete/Age-related/Chow\ et\ al\ 2020/aging-genes-lung.orig.txt", sep = "\t", header = T)
lung_age <- read.table(file = "aging-genes-lung.revised.txt", sep = "\t", header = T)
agegroups <- read.table(file = "mock_agegroups.txt", sep = "\t", header = T)

# Positively associated genes
lung_age_up <- lung_age[lung_age$cluster == "AgeUp",]
tpm_up <- tpm_plots[lung_age_up$genes, ]
tpm_up <- tpm_up[complete.cases(tpm_up),] #Remove NAs
zscore_up <- data.frame(scale(t(tpm_up)))
zscore_up <- data.frame(t(zscore_up))
zscore_up$geneID <- row.names(zscore_up)
zscore_up_long <- zscore_up %>% gather(ID, zscore, Mock_CH_B_376AC_B2_t1:Mock_YA_N_SB_B2_t1)
zscore_up_long <- merge(zscore_up_long, agegroups, by = "ID")
zscore_up_long$AgeGroup <- factor(zscore_up_long$AgeGroup, levels = c("B","CH","YA","OA"))

stat_test <- zscore_up_long %>% 
  wilcox_test(zscore ~ AgeGroup, p.adjust.method = "bonferroni")
stat_test

ggplot(zscore_up_long, aes(x = AgeGroup, y = zscore, fill = AgeGroup)) + 
  geom_boxplot() +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  ylim(-3,4) #export 4.5x6

# Negatively associated genes

lung_age_down <- lung_age[lung_age$cluster == "AgeDown",]
tpm_down <- tpm_plots[lung_age_down$genes, ]
tpm_down <- tpm_down[complete.cases(tpm_down),]
zscore_down <- data.frame(scale(t(tpm_down)))
zscore_down <- data.frame(t(zscore_down))
zscore_down$geneID <- row.names(zscore_down)
zscore_down_long <- zscore_down %>% gather(ID, zscore, Mock_CH_B_376AC_B2_t1:Mock_YA_N_SB_B2_t1)
zscore_down_long <- merge(zscore_down_long, agegroups, by = "ID")
zscore_down_long$AgeGroup <- factor(zscore_down_long$AgeGroup, levels = c("B","CH","YA","OA"))

stat_test <- zscore_down_long %>% 
  wilcox_test(zscore ~ AgeGroup, p.adjust.method = "bonferroni", paired = F)
stat_test

ggplot(zscore_down_long, aes(x = AgeGroup, y = zscore, fill = AgeGroup)) + 
  geom_boxplot() +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  ylim(-3,4)

compare_means(AgeGroup ~ zscore, data = zscore_down_long, paired = TRUE)

# Positively associated genes Individual ages
lung_age_up <- lung_age[lung_age$cluster == "AgeUp",]
tpm_up <- tpm_plots[lung_age_up$genes, ]
tpm_up <- tpm_up[complete.cases(tpm_up),] #Remove NAs
tpm_up$geneID <- row.names(tpm_up)
tpm_up_long <- tpm_up %>% gather(ID, zscore, Mock_CH_B_376AC_B2_t1:Mock_YA_N_SB_B2_t1)
tpm_up_long <- merge(tpm_up_long, agegroups, by = "ID")

ggplot(tpm_up_long, aes(x = Age, y = zscore)) +
  geom_boxplot(aes(group = Age)) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  ylim(-3,4) #export 4.5x6

# Negatively associated genes Individual ages
lung_age_down <- lung_age[lung_age$cluster == "AgeDown",]
tpm_down <- tpm_plots[lung_age_down$genes, ]
tpm_down <- tpm_down[complete.cases(tpm_down),] #Remove NAs
tpm_down$geneID <- row.names(tpm_down)
tpm_down_long <- tpm_down %>% gather(ID, zscore, Mock_CH_B_376AC_B2_t1:Mock_YA_N_SB_B2_t1)
tpm_down_long <- merge(tpm_down_long, agegroups, by = "ID")

ggplot(tpm_down_long, aes(x = Age, y = zscore)) +
  geom_boxplot(aes(group = Age)) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  ylim(-3,4) #export 4.5x6

