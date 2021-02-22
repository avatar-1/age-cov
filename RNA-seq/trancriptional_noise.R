# This script was used to determine transcriptional noise within each age group for the COVID-19 paper #
# Code developed by Alexander Capraro, September 2020
library(dplyr)
library(Hmisc)
library(EnvStats)
library(ggpubr)

#### Mock ####

# 1. Get CPM for all genes and filter out no expression and unwanted samples
cpm <- read.table(file = "cpm_for_sc.txt", sep = "\t", header = T, row.names = 1)
cpm <- cpm[,c(1:16,18:40,42:44)] #Remove BC_B1
cpm_plots <- cpm[,grepl('Mock', colnames(cpm))] #only keep Mock

# 2. Divide genes into 10 equally sized bins based on mean expression and exclude bin 1 and 10

cpm_mean <- data.frame(rowMeans(cpm_plots))
row.names(cpm_mean) <- row.names(cpm_plots)
cpm_mean <- data.frame(cpm_mean[order(-cpm_mean$rowMeans.cpm_plots.), , drop = FALSE])
cpm_mean$bins <- ntile(cpm_mean,10) #creates 10 equally sized bins
cpm_mean <- cpm_mean[cpm_mean$bins != "10", ] #remove top bin
cpm_mean <- cpm_mean[cpm_mean$bins != "1", ] #remove bottom bin

# 3. Calculate co-efficent of variation for all genes

cpm_retained <- cpm_plots[rownames(cpm_mean), ]
cpm_mean$sd <- apply(cpm_retained[,1:14],1,sd)
cpm_mean$cv <- (cpm_mean$sd / cpm_mean$rowMeans.cpm_plots. * 100)
cpm_mean <- data.frame(cpm_mean[order(-cpm_mean$cv), , drop = FALSE])
cpm_mean$cv_bins <- ntile(cpm_mean$cv,10) #split into 10 groups based on co-efficient of variation
low_cv <- cpm_mean[cpm_mean$cv_bins == "1", ] #only keep 10% of genes with lowest CV

# 4. Subset cpm with new list

cpm_retained <- cpm_retained[rownames(low_cv), ] 

# 4. Calculate euclidean distance between each sample & the sample mean within each age group

cpm_B <- (cpm_retained[, 1:5])
cpm_B$mean <- rowMeans(cpm_B)
cpm_CH <- (cpm_retained[, 6:10])
cpm_CH$mean <- rowMeans(cpm_CH)
cpm_YA <- (cpm_retained[, 16:19])
cpm_YA$mean <- rowMeans(cpm_YA)
cpm_OA <- (cpm_retained[, 11:15])
cpm_OA$mean <- rowMeans(cpm_OA)

euclideandistance_B <- dist(t(cpm_B))
euclideandistance_CH <- dist(t(cpm_CH))
euclideandistance_YA <- dist(t(cpm_YA))
euclideandistance_OA <- dist(t(cpm_OA))

euclideandistance_B
euclideandistance_CH
euclideandistance_YA
euclideandistance_OA


#### Virus Infected ####

# 1. Get CPM for all genes and filter out no expression and unwanted samples
cpm <- read.table(file = "cpm_for_sc.txt", sep = "\t", header = T, row.names = 1)
cpm <- cpm[,c(1:16,18:40,42:44)] #Remove BC_B1
cpm_plots <- cpm[,grepl('Virus', colnames(cpm))] #only keep Mock

# 2. Divide genes into 10 equally sized bins based on mean expression and exlude bin 1 and 10

cpm_mean <- data.frame(rowMeans(cpm_plots))
row.names(cpm_mean) <- row.names(cpm_plots)
cpm_mean <- data.frame(cpm_mean[order(-cpm_mean$rowMeans.cpm_plots.), , drop = FALSE])
cpm_mean$bins <- ntile(cpm_mean,10) #creates 10 equally sized bins
cpm_mean <- cpm_mean[cpm_mean$bins != "10", ] #remove top bin
cpm_mean <- cpm_mean[cpm_mean$bins != "1", ] #remove bottom bin

# 3. Calculate co-efficent of variation for all genes

cpm_retained <- cpm_plots[rownames(cpm_mean), ]
cpm_mean$sd <- apply(cpm_retained[,1:14],1,sd)
cpm_mean$cv <- (cpm_mean$sd / cpm_mean$rowMeans.cpm_plots. * 100)
cpm_mean <- data.frame(cpm_mean[order(-cpm_mean$cv), , drop = FALSE])
cpm_mean$cv_bins <- ntile(cpm_mean$cv,10) #creates 10 equally sized bins
low_cv <- cpm_mean[cpm_mean$cv_bins == "1", ] #only keep 10% of genes with lowest CV

# 4. Subset cpm with new list

cpm_retained <- cpm_retained[rownames(low_cv), ] 

# 4. Calculate euclidean distance between each sample & hte sample mean within each age group

cpm_B <- (cpm_retained[, 1:5])
cpm_B$mean <- rowMeans(cpm_B)
cpm_CH <- (cpm_retained[, 6:13])
cpm_CH$mean <- rowMeans(cpm_CH)
cpm_YA <- (cpm_retained[, 20:23])
cpm_YA$mean <- rowMeans(cpm_YA)
cpm_OA <- (cpm_retained[, 14:19])
cpm_OA$mean <- rowMeans(cpm_OA)

euclideandistance_B <- dist(t(cpm_B))
euclideandistance_CH <- dist(t(cpm_CH))
euclideandistance_YA <- dist(t(cpm_YA))
euclideandistance_OA <- dist(t(cpm_OA))

euclideandistance_B
euclideandistance_CH
euclideandistance_YA
euclideandistance_OA

#### Combined ####

# 1. Get CPM for all genes and filter out no expression and unwanted samples
cpm <- read.table(file = "cpm_for_sc.txt", sep = "\t", header = T, row.names = 1)
cpm_plots <- cpm[,c(1:16,18:40,42:44)] #Remove BC_B1

# 2. Divide genes into 10 equally sized bins based on mean expression and exlude bin 1 and 10

cpm_mean <- data.frame(rowMeans(cpm_plots))
row.names(cpm_mean) <- row.names(cpm_plots)
cpm_mean <- data.frame(cpm_mean[order(-cpm_mean$rowMeans.cpm_plots.), , drop = FALSE])
cpm_mean$bins <- ntile(cpm_mean,10) #creates 10 equally sized bins
cpm_mean <- cpm_mean[cpm_mean$bins != "10", ] #remove top bin
cpm_mean <- cpm_mean[cpm_mean$bins != "1", ] #remove bottom bin

# 3. Calculate co-efficent of variation for all genes

cpm_retained <- cpm_plots[rownames(cpm_mean), ]
cpm_mean$sd <- apply(cpm_retained[,1:14],1,sd)
cpm_mean$cv <- (cpm_mean$sd / cpm_mean$rowMeans.cpm_plots. * 100)
cpm_mean <- data.frame(cpm_mean[order(-cpm_mean$cv), , drop = FALSE])
cpm_mean$cv_bins <- ntile(cpm_mean$cv,10) #creates 10 equally sized bins
low_cv <- cpm_mean[cpm_mean$cv_bins == "1", ] #only keep 10% of genes with lowest CV

# 4. Subset cpm with new list

cpm_retained <- cpm_retained[rownames(low_cv), ] 

# 4. Calculate euclidean distance between each sample & the sample mean within each age group

cpm_B_mock <- log2(cpm_retained[, 1:5])
cpm_B_mock$mean <- rowMeans(cpm_B_mock)
cpm_CH_mock <- log2(cpm_retained[, 6:10])
cpm_CH_mock$mean <- rowMeans(cpm_CH_mock)
cpm_YA_mock <- log2(cpm_retained[, 16:19])
cpm_YA_mock$mean <- rowMeans(cpm_YA_mock)
cpm_OA_mock <- log2(cpm_retained[, 11:15])
cpm_OA_mock$mean <- rowMeans(cpm_OA_mock)

cpm_B_virus <- log2(cpm_retained[, 20:24])
cpm_B_virus$mean <- rowMeans(cpm_B_virus)
cpm_CH_virus <- log2(cpm_retained[, 25:32])
cpm_CH_virus$mean <- rowMeans(cpm_CH_virus)
cpm_YA_virus <- log2(cpm_retained[, 39:42])
cpm_YA_virus$mean <- rowMeans(cpm_YA_virus)
cpm_OA_virus <- log2(cpm_retained[, 33:38])
cpm_OA_virus$mean <- rowMeans(cpm_OA_virus)

euclideandistance_B_mock <- dist(t(cpm_B_mock))
euclideandistance_CH_mock <- dist(t(cpm_CH_mock))
euclideandistance_YA_mock <- dist(t(cpm_YA_mock))
euclideandistance_OA_mock <- dist(t(cpm_OA_mock))

euclideandistance_B_virus <- dist(t(cpm_B_virus))
euclideandistance_CH_virus <- dist(t(cpm_CH_virus))
euclideandistance_YA_virus <- dist(t(cpm_YA_virus))
euclideandistance_OA_virus <- dist(t(cpm_OA_virus))

euclideandistance_B_mock
euclideandistance_CH_mock
euclideandistance_YA_mock
euclideandistance_OA_mock
euclideandistance_B_virus
euclideandistance_CH_virus
euclideandistance_YA_virus
euclideandistance_OA_virus


#### Boxplots ####

#The distance between the group mean was used for further analyses

mock_distance <- read.table(file = "mock_all_distances.txt", sep = "\t", header = T)
mock_distance$AgeGroup <- factor(x = mock_distance$AgeGroup, levels = c("OA", "YA", "CH", "B")) 

virus_distance <- read.table(file = "virus_all_distances.txt", sep = "\t", header = T)
virus_distance$AgeGroup <- factor(x = virus_distance$AgeGroup, levels = c("OA", "YA", "CH", "B")) #Keeps plot in order of child > YA > OA

#Mock

stat_test <- mock_distance %>%
  wilcox_test(Log ~ AgeGroup, p.adjust.method = "none")
stat_test

ggplot(mock_distance, aes(x = AgeGroup, y = Log, fill = AgeGroup)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.75), dotsize = 0.5) +
  scale_fill_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7","#F7941D")) +
  theme_minimal() +
  ylim(6,10) #export pdf 6x9

#Virus

stat_test <- virus_distance %>%
  wilcox_test(Log ~ AgeGroup, p.adjust.method = "bonferroni")
stat_test

ggplot(virus_distance, aes(x = AgeGroup, y = Log, fill = AgeGroup)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.75), dotsize = 0.5) +
  scale_fill_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7","#F7941D")) +
  theme_minimal() +
  ylim(6,10) #export pdf 6x9
