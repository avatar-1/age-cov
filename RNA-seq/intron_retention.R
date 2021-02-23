#Code is used for analysis of rMATS output for retained introns for COVID-19 Age paper
#Code developed by Alexander Capraro, February 2021

library(tidyverse)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(viridis)
library(rmulti)
library(rstatix)
library(dplyr)
library(ggpubr)
library(matrixStats)

#### Intron retention ratios ####

id_geneID <- read.table(file = "id_geneID.txt", sep = "\t", header = T)
agegroups <- read.table(file = "boxplot_zscore_groups_All.txt", sep = "\t", header = T)
norm_vector <- c(0.93493296,1.015311128,0.921509557,0.908189626,0.976351798,0.897252521,1.081556438,
                 1.124373232,1.178587361,1.035609372,1.007446121,0.852932689,0.874325679,0.904038206,
                 1.019615276,0.908572972,1.167623835,0.945156585,0.927774692,0.950327209,1.117936928,
                 1.091069299,0.924244042,0.900502024,1.036471925,1.060013502,1.151830817,0.947322641,
                 0.929802136,1.117854478,1.071082279,1.053706555,1.138821439,0.994728251,0.979801373,
                 0.93082647,1.049918095,1.009018197,0.965349277,0.994720103,1.11106527,0.903699422) #Normalise based on readcounts per sample (total readcount/mean readcount)

RI_positive <- read.table(file = "RI_positive.txt", sep = "\t", header = T, row.names = 1)
RI_negative <-read.table(file = "RI_negative.txt", sep = "\t", header = T, row.names = 1)

positive_matrix <- as.matrix(RI_positive[,2:43]) #Create matrix and remove geneIDs
negative_matrix <- as.matrix(RI_negative[,2:43]) #Create matrix and remove geneIDs

positive_matrix_norm <- positive_matrix*norm_vector[col(positive_matrix)]
negative_matrix_norm <- negative_matrix*norm_vector[col(negative_matrix)]

#filt <- as.matrix(positive_matrix[rowSums(positive_matrix,na.rm=TRUE)>9,]) #filter RI sites with 10 or more counts across all samples
#positive_matrix <- positive_matrix[row.names(filt),]
#negative_matrix <- negative_matrix[row.names(filt),]

RI_ratio <- positive_matrix_norm / (negative_matrix_norm + positive_matrix_norm)
RI_ratio <- as.data.frame(RI_ratio)

ratio_genes <- merge(RI_ratio, id_geneID, by.x = "row.names", by.y = "ID") #re-add geneID
ratio_genes <- ratio_genes[,c(44,2:43)]

ratio_plot <- reshape2::melt(ratio_genes)
ratio_plot <- ratio_plot[order(ratio_plot$GeneID),] #order by gene

ratio_plot <- merge(ratio_plot, agegroups, by.x = "variable", by.y = "ID")
ratio_plot$Condition <- factor(ratio_plot$Condition, levels = c("OA_Mock","YA_Mock","CH_Mock","B_Mock",
                                                                  "OA_Virus","YA_Virus","CH_Virus","B_Virus"))

stat_test <- ratio_plot %>%
  wilcox_test(value ~ Condition, p.adjust.method = "bonferroni")
stat_test

ggplot(ratio_plot, aes(x = Condition, y = value, fill = Condition)) + 
  geom_violin() +
  theme_minimal() +
  scale_fill_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7","#F7941D","#BEBADA","#FFFFB3","#8DD3C7","#F7941D"))  +
  stat_summary(fun.y=median, geom="point", shape=16, size=2, position = position_dodge(0.9)) +
  ylab("Retained Intron Ratio per gene (%)")

#### Correlatin with viral load ####

mean_ratio <- as.data.frame(colMeans(RI_ratio, na.rm = TRUE))

viral_reads <- read.table(file = "viral_tpm.txt", sep = "\t", header = T)

mean_ratio <- merge(mean_ratio, viral_reads, by.x = "row.names", by.y = "ID")

corr_test <- cor.test(x=mean_ratio[,3], y=mean_ratio[,2], method = 'spearman')
corr_test

ggplot(mean_ratio, aes(x = mean_ratio[,3], y = mean_ratio[,2])) +
  geom_point(color = "#00AFBB") +
  theme_minimal() +
  geom_smooth(method=lm, color="#00AFBB", fill = "#00AFBB") +
  xlim(0,15) + ylim(0.35, 0.5) #4.5x4
