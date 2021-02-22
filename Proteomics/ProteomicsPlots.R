library(tidyverse)
library(forcats)
library(ggplot2)
library(gplots)
library(viridis)
library(ggthemes)
library(rstatix)
library(ggpubr)

#### RNA-Protein Correlation plots ####

Mock_corr <- read.table(file = "Mock_All_RNA+Protein.txt", sep = "\t", header = T) #imports log2FC for nasal CH and bronchial. Only contains genes expressed in both cell types

corr <- cor.test(x=Mock_corr$Protein, y=Mock_corr$RNA, method = 'spearman')
corr
ggplot(Mock_corr, aes(x = Protein, y = RNA)) +
  geom_point() +
  theme_minimal() +
  geom_smooth(method=lm, color='#2C3E50') + 
  ylim(0,15)

Virus_corr <- read.table(file = "Virus_All_RNA+Protein.txt", sep = "\t", header = T) #imports log2FC for nasal CH and bronchial. Only contains genes expressed in both cell types

corr <- cor.test(x=Virus_corr$Protein, y=Virus_corr$RNA, method = 'spearman')
corr
ggplot(Virus_corr, aes(x = Protein, y = RNA)) +
  geom_point() +
  theme_minimal() +
  geom_smooth(method=lm, color='#2C3E50') + 
  ylim(0,15)

#### Age associated genes  ####

prot <- read.table(file = "lfq_matrix.txt", sep = "\t", header = T)
prot <- prot[,c(1:11)] #Mock nasal

lung_age <- read.table(file = "aging-genes-lung.revised.txt", sep = "\t", header = T)
agegroups <- read.table(file = "mock_agegroups.txt", sep = "\t", header = T)

# Positively associated genes
lung_age_up <- lung_age[lung_age$cluster == "AgeUp",]
prot_up <- prot[lung_age_up$genes, ]
prot_up <- prot_up[complete.cases(prot_up),] #Remove NAs
zscore_up <- data.frame(scale(t(prot_up)))
zscore_up <- data.frame(t(zscore_up))
zscore_up$geneID <- row.names(zscore_up)
zscore_up_long <- zscore_up %>% gather(ID, zscore, M_B_376AC:M_YA_SB)
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
prot_down <- prot[lung_age_down$genes, ]
prot_down <- prot_down[complete.cases(prot_down),]
zscore_down <- data.frame(scale(t(prot_down)))
zscore_down <- data.frame(t(zscore_down))
zscore_down$geneID <- row.names(zscore_down)
zscore_down_long <- zscore_down %>% gather(ID, zscore, M_B_376AC:M_YA_SB)
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

#### Mock Differentially abundant proteins HEATMAP ####

prot <- read.table(file = "lfq_matrix.txt", sep = "\t", header = T)
prot <- prot[,c(3:11)] #Mock nasal

heatmap_average <- function(genes){
  heatmap_plot <- data.frame((subset(prot, rownames(prot) %in% genes))) #Subset only genes in list
  average_plot <- data.frame(row.names = row.names(heatmap_plot))
  average_plot$Mock_CH <- rowMeans(heatmap_plot[,1:2])
  average_plot$Mock_YA <- rowMeans(heatmap_plot[,6:9])
  average_plot$Mock_OA <- rowMeans(heatmap_plot[,3:5])
  heatmap.2(as.matrix(average_plot),scale = "row", trace = 'none', density.info = 'none', col = viridis)
} #Heatmap of average for each condition + age group 

heatmap_average(genes = c("ITGB4","CFAP73","MYO18A","SURF4","CLIP1","NAMPT","NFKB2","SERPINB1","ALDH4A1","S100A7","SLC25A12","GCHFR","CUL4A","PLCB3","MMP10","TMED2","DDX47","CTSC","ADH1C","HGS","FNTA","CLIC3","PTPN23","CYB5R2","VPS52","TIGAR","RANGAP1","ESRP1","PPIG","DDX21","CD109","ABCC1","PMM2","PTK7","AK3","S100A9","ACOT9","PDZD11","RPS15","SSR1","HNRNPH3","NUP160","PYGL","AKR1C3","TXNRD2","HAGH","NXN","C3","SULT1A1","C1orf194","RPL26","NDRG1","FTO","KDSR","TNFAIP8","ARHGAP1","HNRNPH2","ECHDC1","ARF3","HMGCL","NMT1","CALM3","DYNLT1","SH3GLB1","RPE","CSNK2B","RAE1","BID","DYNLL2","MACROH2A2","MBP","TOM1","SRRM1","PCYOX1","SKIV2L","HMBS","PRDX2","ADD1","RANBP2","NUCKS1","PRMT1","H3C13","RANBP1","LMAN2","RPS18","SCARB2","OAS2","CAB39","UFM1","CSTF2","TSGA10","DNAJC10","NUDT21","NOL3","CPNE1","DDOST","PDLIM1","MX1","CRYBG1","BTF3","SRSF7","LAMTOR1","LRRC59","ALDH6A1","OTUB1","PCMT1","MTPN","PDAP1","TNPO2","SH3BGRL","DLAT","DDI2","SCGB1A1","EPS15","RPL35A","NDUFS5","PCK2","UNC93B1","IVL","SEC14L2","CRNN","ZYX","RSPH1","GPX1","LYPD3","BPGM","MB","MAP1A","DHRS9","BPIFA1","IL1RN","ASS1","CLIC4","ST13","CORO1C","RSPH4A","RAD50","DNAH5","NME1","TWF2","FMO3","IFT172","SAMHD1","DNAH9","CAVIN1","KIF21A","SCAMP1","BLVRA","COX7C","ANXA8L1","CLCA4","CCAR1","MAT2A","SELENBP1","DYNC2H1","PLXNB2","NUP62","GSTA1","ASCC3","RPL35","TTC25","TUBB6","AK7","GLG1","COPS4","PCYT2","SUGT1","RAN","GLIPR2","MFF","MYOF","CORO7","COX5A","LLGL2","PHGDH","MRTO4","GPRC5A","SYPL1","PRPF40A","VCL","RUVBL2","MAGOHB","ADH7","SPAG6","ABLIM1","TPP2","TP53BP1","WASL","SUMF2","CAPS","RSPH9","GCLM","EEF1D","APPL2","CYRIB","PPA1","NUCB2","PARP1","VCPIP1","PSPC1","FAM83H","TGM2","CROCC","LY6D","AKR1B10","LAMB3","TCEA1","PEA15","BOLA2","ANXA8","CAPNS1","LRRC46","RAB10","MRPL16","EDC4","ZFPL1","SPINT1","TPR","ARHGDIB","RPL31","EFHC1","MYO1D","CETN2","DAP","FLNA","SNRPB2","GAPVD1","USP7","SRRT","ITPR3","NUDC","SERPINB13","UBE2K","CPD","ALDH3B1","NCLN","TMPRSS11D","TOM1L2","ILF2","ABHD6","PLS3","SERPINH1","DYNLRB2","PREP","ACAA1","HDGFL2","SEPTIN11","FAM114A1","MAPRE3","TFAM","ACP1","IL18","CFAP52","PYGB","GPC1","CRK","API5","CDC5L","PPM1G","DNALI1","NDUFA9","TMEM43","SLC27A2","NUMA1","UGDH","DARS1","OSCP1","MACF1","DDX1","B2M","PAFAH1B3","ERLEC1","PSAP","ABRACL","LAD1","CSNK1A1","CALD1","CNN2","CFAP20","HSD17B12","GLTP","PDXDC1","SEC24B","SCAMP2","IPO4","CDV3","COPS7B","PTRH2","CHMP4B","ERGIC1","PLRG1","PLEKHA5","G6PD","PAPSS1","YWHAB","TPD52L1","HM13","GOLGB1","ERICH3","WASHC2A","RRBP1","PDCD6","ERH","USP10","OPA1","S100A2","PTPN1","SPATA18","TIMM9","KPNA6","NBAS","ATP5PF","AGR3","GRB2","DDX19B","HSP90AA1","PPL","FUBP1","SUPT16H","GMPR2","TPPP3","PLEKHS1","SRSF1","KIAA0513","FSCN1","MCCC1","CAT","C4B","LGALS3","RPL23A","ST14","SND1","ETFDH","FSTL1","ATL3","LGALS1","EIF4A3","HNRNPA0","HYPK","BCAP29","LZTFL1","NDUFB9","RPL19","NT5C2","RBMX","RPA2","EPS8","SNCA","DNAI1","TOMM34","RPL22","H1-10","SQOR","TNKS1BP1","CD9","MRPL53","F3","GTF3C3","FER1L6","MYO1B","RPS8","DNAJB13","CUL2","HNRNPM","STEAP4","SSH3","REXO2","PCBP2","TUBB4B","HLA-A","GLB1","SUCLG1"))
  
#### Mock IPA ####

canonical <- read.table(file = "CanonicalMockPlot.txt", sep = "\t", header = TRUE, row.names = 1)
heatmap.2(as.matrix(canonical),dendrogram='none', Rowv = NA, Colv = NA, trace = 'none', density.info = 'none', col=brewer.pal(9,"YlOrRd"), cexCol = 1, cexRow = 1)

#### Viral proteins ####

virus <- read.table(file = "ViralProteinLFQ.txt", sep = "\t", header = T)

virus_input <- reshape2::melt(virus)
virus_input <- virus_input[order(virus_input$Gene.names),] #order by gene
virus_input$Age <- ifelse(grepl("CH", virus_input$variable), "Child", 
                        ifelse(grepl("YA", virus_input$variable), "Young Adult", 
                               ifelse(grepl("OA", virus_input$variable), "Older Adult","Bronchial")))
virus_input$Age <- factor(x = virus_input$Age, levels = c("Older Adult","Young Adult","Child", "Bronchial"))
ggplot(virus_input[which(virus_input$value > 0), ], aes (x = Age, y = value, fill = Age)) + 
  geom_boxplot(alpha = 1) +
  geom_point() +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")+ 
  theme(legend.position = "none") +
  facet_wrap(~Gene.names, ncol = 7,
             scales = "fixed") +
  ylab("log2 LFQ intensity")

#horizontal
ggplot(virus_input[which(virus_input$value > 0), ], aes (x = Age, y = value, fill = Age)) + 
  geom_boxplot(alpha = 1) +
  geom_point()+
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")+ 
  theme(legend.position = "none") +
  ylab("log2 LFQ intensity")  +
  coord_flip() +
  facet_wrap(~Gene.names, ncol = 1,
             scales = "fixed") + #facet graph to get multiple groups of boxplots/lines/points
  #change theme elements so graph does not appear facetted
  theme(panel.border = element_blank(), #remove borders on facets
        panel.spacing = unit(0, "lines"), #remove spacing btween panels
        strip.background = element_rect(color = "white"))

stat_test <- virus_input %>%
  group_by(Gene.names) %>%
  t_test(value ~ Age, p.adjust.method = "bonferroni", paired = F)
stat_test

#### Mock vs Virus IPA Comparison ####

canonical <- read.table(file = "Canonical_MockvVirus.txt", sep = "\t", header = TRUE, row.names = 1)
heatmap.2(as.matrix(canonical),dendrogram='none', Rowv = NA, Colv = NA, trace = 'none', density.info = 'none', col=brewer.pal(9,"YlOrRd"), cexCol = 1, cexRow = 1)

#### ISG heatmap ####
prot <- read.table(file = "log2_lfq_matrix.txt", sep = "\t", header = T, row.names = 1)

heatmap_average <- function(genes){
  heatmap_plot <- data.frame((subset(prot, rownames(prot) %in% genes))) #Subset only genes in list
  average_plot <- data.frame(row.names = row.names(heatmap_plot))
  average_plot$Mock_OA <- rowMeans(heatmap_plot[,5:7])
  average_plot$Mock_YA <- rowMeans(heatmap_plot[,8:11])
  average_plot$Mock_CH <- rowMeans(heatmap_plot[,3:4])
  average_plot$Mock_B <- rowMeans(heatmap_plot[,1:2])
  average_plot$Virus_OA <- rowMeans(heatmap_plot[,18:20])
  average_plot$Virus_YA <- rowMeans(heatmap_plot[,21:24])
  average_plot$Virus_CH <- rowMeans(heatmap_plot[,15:17])
  average_plot$Virus_B <- rowMeans(heatmap_plot[,12:14])
  Mean <- summarise_all(average_plot, mean)
  average_plot <- rbind(average_plot, Mean)
  rownames(average_plot)[rownames(average_plot) == "1"] <- "Mean"
  heatmap.2(as.matrix(average_plot),scale = "row", Rowv = F, Colv = F, breaks=seq(-1.5,1.5,0.2), trace = 'none', density.info = 'none', col = viridis)
} #Heatmap of average for each condition + age group 
boxplot_zscore <- function(genes, scales, cols){
  agegroups <- read.table(file = "agegroups.txt", sep = "\t", header = T)
  input_prot <- data.frame((subset(prot, rownames(prot) %in% genes))) #Subset only genes in list
  input_prot <- input_prot[complete.cases(input_prot),]
  zscore <- data.frame(scale(t(input_prot)))
  zscore <- data.frame(t(zscore))
  zscore$geneID <- row.names(zscore)
  zscore_long <- zscore %>% gather(ID, zscore, M_B_376AC:V_YA_SB)
  zscore_long <- merge(zscore_long, agegroups, by = "ID")
  zscore_long$AgeGroup <- factor(zscore_long$AgeGroup, levels = c("OA_Mock","YA_Mock","CH_Mock","B_Mock","OA_Virus","YA_Virus","CH_Virus","B_Virus"))
  
  .GlobalEnv$stat_test <- zscore_long %>% 
    wilcox_test(zscore ~ AgeGroup, p.adjust.method = "bonferroni", paired = F)
  print(stat_test)
  
  plot<-  ggplot(zscore_long, aes(x = AgeGroup, y = zscore, fill = AgeGroup)) + 
    geom_boxplot() +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    ylim(-3,4)
  #export 4.5x6
  return(list(plot=plot, stat_test = stat_test))
}

heatmap_average(genes = c("BLVRA","PARP4","PPA1",
                          "ASCC3","ATP6V1H","CASP3","NME7","NT5C2","RANBP2","VCPIP1",
                          "ACOT9","B2M","BAG1","MX1","OAS2","PARP1","SAMHD1","SCARB2","UNC93B1","ZC3HAV1","GBP6",
                          "PML",
                          "CLIC4","IL1RN","MYOF","NAMPT","RRBP1")) # ALL ISGs

heatmap_average(genes = c("PARP4","BLVRA","CAMK2D","GCLC","CAST","ANXA1","SLC25A24","CAPZA2","KPNB1","PPA1","SNX6","TOR1AIP1","KIF13B","CAPN2","ERP44",
                          "LRBA","SSB","SNX2","ABI1","DECR1","DNAJA1","CUL1","RAB3GAP1","NME7","GPD2","APPL1","SUB1","UBA6","UTRN","CASP3","RANBP2","NT5C2","XRN2","ATP6V1H","CD2AP",
                          "DTX3L","GBP6","ZC3HAV1","ACOT9","STAT2","STAT3","BAG1","UNC93B1","CNP","PARP1","MX1","CASP1","ADAR","LAP3","LGALS3BP","SLFN5","SAMHD1","PSME1","SPATS2L","AGRN","ANKFY1","STAT1","CHMP5","PSME2","SCARB2","TRIM25",
                          "RABGGTA","COASY","OGFR","ABCC1","PML","NAPA","ATP13A1",
                          "LMNB1","MYOF","MUC1","EHD4","RALB","RRBP1","ACSL1","GLRX","CKAP4","CASP7","PGD","KYNU","NAGK","STOM","MARCKS","MVP","CLIC4","IL1RN","GSTO1","PRKCD","LYN")) # ALL ISGs

heatmap_average(genes = c("PARP4","BLVRA","CAMK2D","GCLC","CAST","ANXA1","SLC25A24","CAPZA2","KPNB1","PPA1","SNX6","TOR1AIP1","KIF13B","CAPN2","ERP44",
                          "LRBA","SSB","SNX2","ABI1","DECR1","DNAJA1","CUL1","RAB3GAP1","NME7","GPD2","APPL1","SUB1","UBA6","UTRN","CASP3","RANBP2","NT5C2","XRN2","ATP6V1H","CD2AP")) #C1 & C2
boxplot_zscore(genes = c("PARP4","BLVRA","CAMK2D","GCLC","CAST","ANXA1","SLC25A24","CAPZA2","KPNB1","PPA1","SNX6","TOR1AIP1","KIF13B","CAPN2","ERP44",
                         "LRBA","SSB","SNX2","ABI1","DECR1","DNAJA1","CUL1","RAB3GAP1","NME7","GPD2","APPL1","SUB1","UBA6","UTRN","CASP3","RANBP2","NT5C2","XRN2","ATP6V1H","CD2AP"))

heatmap_average(genes = c("DTX3L","GBP6","ZC3HAV1","ACOT9","STAT2","STAT3","BAG1","UNC93B1","CNP","PARP1","MX1","CASP1","ADAR","LAP3","LGALS3BP","SLFN5","SAMHD1",
                          "PSME1","SPATS2L","AGRN","ANKFY1","STAT1","CHMP5","PSME2","SCARB2","TRIM25")) # C3
boxplot_zscore(genes = c("DTX3L","GBP6","ZC3HAV1","ACOT9","STAT2","STAT3","BAG1","UNC93B1","CNP","PARP1","MX1","CASP1","ADAR","LAP3","LGALS3BP","SLFN5","SAMHD1",
                         "PSME1","SPATS2L","AGRN","ANKFY1","STAT1","CHMP5","PSME2","SCARB2","TRIM25"))

heatmap_average(genes = c("RABGGTA","COASY","OGFR","ABCC1","PML","NAPA","ATP13A1")) #C4
boxplot_zscore(genes = c("RABGGTA","COASY","OGFR","ABCC1","PML","NAPA","ATP13A1")) #C4

heatmap_average(genes = c("LMNB1","MYOF","MUC1","EHD4","RALB","RRBP1","ACSL1","GLRX","CKAP4","CASP7","PGD","KYNU","NAGK","STOM","MARCKS","MVP","CLIC4","IL1RN","GSTO1","PRKCD","LYN")) #C5
boxplot_zscore(genes = c("LMNB1","MYOF","MUC1","EHD4","RALB","RRBP1","ACSL1","GLRX","CKAP4","CASP7","PGD","KYNU","NAGK","STOM","MARCKS","MVP","CLIC4","IL1RN","GSTO1","PRKCD","LYN")) #C5


#### Canonical pathway virus age groups vs. ALL ####

# Child
child <- read.table("CHvAll_Virus_Canonical.txt", sep ="\t", header = TRUE) #Import IPA child canonical pathway
child$Ingenuity.Canonical.Pathways <- factor(child$Ingenuity.Canonical.Pathways, levels = child$Ingenuity.Canonical.Pathways[order(child$X.log.p.value.)])
ggplot(child, aes(y=Ingenuity.Canonical.Pathways, x=X.log.p.value., fill=zScore)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradient2_tableau(palette = "Orange-Blue Diverging",limits = c(-3,3)) #limits sets the scale bar for fill colour

#YA
YA <- read.table("YAvAll_Virus_Canonical.txt", sep ="\t", header = TRUE) #Import IPA child canonical pathway
YA$Ingenuity.Canonical.Pathways <- factor(YA$Ingenuity.Canonical.Pathways, levels = YA$Ingenuity.Canonical.Pathways[order(YA$X.log.p.value.)])
ggplot(YA, aes(y=Ingenuity.Canonical.Pathways, x=X.log.p.value., fill=zScore)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradient2_tableau(palette = "Orange-Blue Diverging",limits = c(-3,3)) #limits sets the scale bar for fill colour

#OA
OA <- read.table("OAvAll_Virus_Canonical.txt", sep ="\t", header = TRUE) #Import IPA child canonical pathway
OA$Ingenuity.Canonical.Pathways <- factor(OA$Ingenuity.Canonical.Pathways, levels = OA$Ingenuity.Canonical.Pathways[order(OA$X.log.p.value.)])
ggplot(OA, aes(y=Ingenuity.Canonical.Pathways, x=X.log.p.value., fill=zScore)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradient2_tableau(palette = "Orange-Blue Diverging",limits = c(-3,3)) #limits sets the scale bar for fill colour))

#### Cell-type marker ####

cell_markers <- read.table(file = "cellmarkers.txt", sep = "\t", header = T)

#### MOCK 

prot <- read.table(file = "lfq_matrix.txt", sep = "\t", header = T)
prot <- prot[,c(1:11)] #Mock

input_prot <- data.frame((subset(prot, rownames(prot) %in% cell_markers$Protein))) #Subset only genes in list
input_prot$gene <- row.names(input_prot)
input_prot <- reshape2::melt(input_prot)
input_prot <- input_prot[order(input_prot$gene),] #order by gene

input_prot$Age <- ifelse(grepl("CH", input_prot$variable), "Child", 
                         ifelse(grepl("YA", input_prot$variable), "Young Adult", 
                                ifelse(grepl("OA", input_prot$variable), "Older Adult", "Bronchial")))
input_prot$Age <- factor(x = input_prot$Age, levels = c("Older Adult","Young Adult", "Child", "Bronchial")) #Keeps plot in order of child > YA > OA
input_prot <- merge(input_prot, cell_markers, by.x = "gene", by.y = "Protein")

ggplot(input_prot, aes (x = Celltype, y = log2(value), fill = Age)) + 
  geom_violin(alpha = 1) +
  theme_minimal() +
  ylim(22.5,37.5) +
  scale_fill_manual(values = c("#BFBADA", "#FCF6B5", "#8ED2C6", "#F89521"))+ 
  theme(legend.position = "none") +
  ylab("log2 LFQ intensity") +
  stat_summary(fun=mean, geom="point", shape=16, size=2, position = position_dodge(0.9))#5.5x7

stat_test <- input_prot %>% group_by(Celltype) %>%
  wilcox_test(value ~ Age)
print(stat_test)

### VIRUS

prot <- read.table(file = "lfq_matrix.txt", sep = "\t", header = T)
prot <- prot[,c(12:24)] #Virus

input_prot <- data.frame((subset(prot, rownames(prot) %in% cell_markers$Protein))) #Subset only genes in list
input_prot$gene <- row.names(input_prot)
input_prot <- reshape2::melt(input_prot)
input_prot <- input_prot[order(input_prot$gene),] #order by gene

input_prot$Age <- ifelse(grepl("CH", input_prot$variable), "Child", 
                         ifelse(grepl("YA", input_prot$variable), "Young Adult", 
                                ifelse(grepl("OA", input_prot$variable), "Older Adult", "Bronchial")))
input_prot$Age <- factor(x = input_prot$Age, levels = c("Older Adult","Young Adult", "Child", "Bronchial")) #Keeps plot in order of child > YA > OA
input_prot <- merge(input_prot, cell_markers, by.x = "gene", by.y = "Protein")

ggplot(input_prot, aes (x = Celltype, y = log2(value), fill = Age)) + 
  geom_violin(alpha = 1) +
  theme_minimal() +
  ylim(22.5,37.5) +
  scale_fill_manual(values = c("#BFBADA", "#FCF6B5", "#8ED2C6", "#F89521"))+ 
  theme(legend.position = "none") +
  ylab("log2 LFQ intensity") +
  stat_summary(fun.y=mean, geom="point", shape=16, size=2, position = position_dodge(0.9))#5.5x7

stat_test <- input_prot %>% group_by(Celltype) %>%
  wilcox_test(value ~ Age)
print(stat_test)



##### Boxplots & Heatmaps ####

prot <- read.table(file = "lfq_matrix.txt", sep = "\t", header = T)
prot <- prot[,c(3:11,15:24)] #Mock nasal

heatmap_average <- function(genes){
  heatmap_plot <- data.frame((subset(prot, rownames(prot) %in% genes))) #Subset only genes in list
  average_plot <- data.frame(row.names = row.names(heatmap_plot))
  average_plot$Mock_OA <- rowMeans(heatmap_plot[,5:7])
  average_plot$Mock_YA <- rowMeans(heatmap_plot[,8:11])
  average_plot$Mock_CH <- rowMeans(heatmap_plot[,3:4])
  average_plot$Mock_B <- rowMeans(heatmap_plot[,1:2])
  average_plot$Virus_OA <- rowMeans(heatmap_plot[,18:20])
  average_plot$Virus_YA <- rowMeans(heatmap_plot[,21:24])
  average_plot$Virus_CH <- rowMeans(heatmap_plot[,15:17])
  average_plot$Virus_B <- rowMeans(heatmap_plot[,12:14])
  Mean <- summarise_all(average_plot, mean)
  average_plot <- rbind(average_plot, Mean)
  rownames(average_plot)[rownames(average_plot) == "1"] <- "Mean"
  heatmap.2(as.matrix(average_plot),scale = "row", Rowv = F, Colv = F, breaks=seq(-1.5,1.5,0.2), trace = 'none', density.info = 'none', col = viridis)
 }#Heatmap of average for each condition + age group 
boxplot_zscore <- function(genes, scales, cols){
  agegroups <- read.table(file = "agegroups.txt", sep = "\t", header = T)
  input_prot <- data.frame((subset(prot, rownames(prot) %in% genes))) #Subset only genes in list
  input_prot <- input_prot[complete.cases(input_prot),]
  zscore <- data.frame(scale(t(input_prot)))
  zscore <- data.frame(t(zscore))
  zscore$geneID <- row.names(zscore)
  zscore_long <- zscore %>% gather(ID, zscore, M_B_376AC:V_YA_SB)
  zscore_long <- merge(zscore_long, agegroups, by = "ID")
  zscore_long$AgeGroup <- factor(zscore_long$AgeGroup, levels = c("OA_Mock","YA_Mock","CH_Mock","B_Mock","OA_Virus","YA_Virus","CH_Virus","B_Virus"))
  
  .GlobalEnv$stat_test <- zscore_long %>% 
    wilcox_test(zscore ~ AgeGroup, p.adjust.method = "bonferroni", paired = F)
  print(stat_test)
  
  plot<-  ggplot(zscore_long, aes(x = AgeGroup, y = zscore, fill = AgeGroup)) + 
    geom_boxplot() +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    ylim(-3,4)
  #export 4.5x6
  return(list(plot=plot, stat_test = stat_test))
}

boxplot_prot(genes = c("ACE","ACE2","TMPRSS2","TMPRSS4","CTSL"),scales = "free", cols = 4)

boxplot_prot(genes = c("EPPK1","TRAF3IP1","DHRS9","DDOST","CES2","LAD1","GAS8","CFAP53","GAR1","GBA","EHD1","LSM5","TBC1D5","NUP133","EPS8L1","CFAP20","SLC11A2","OAS2","RPS17","TMED9"),scales = "free", cols = 4) 
boxplot_prot(genes = c("SCGB1A1"),scales = "free", cols = 4)
