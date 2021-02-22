# This code is used to analyse lysate proteomics for the COVID-19 study
# Code developed by Alexander Capraro August 2020
library(ggplot2)
library(gplots)
library(tidyr)
library(gridExtra)
library(RColorBrewer)
library(ggfortify)
library(DEP)
library(factoextra)
library(UpSetR)
library(EnhancedVolcano)
options(scipen=999) #removes scientific notation

#### Child Mock vs. Virus ####

data <- read.table(file = "proteinGroups.txt", sep = "\t", header = TRUE)

data <- data[data$Reverse != "+" & data$Potential.contaminant != "+" & data$Unique.peptides >= 2,]
data <- data[data$Gene.names != "",] #removes rows with no Gene names
data <- data[complete.cases(data),] #removes NAs at the end 
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";") 

#Make summarised experiment
columns <- grep("LFQ.intensity.", colnames(data_unique))
exp_design <- read.table(file = "experimental_design_CH.txt", sep = "\t", header = TRUE) #Import experimental design - the same as LFQ-Analyst
se <- make_se(data_unique, columns, exp_design)

# Filter and normalise
filt <- filter_proteins(se, "fraction", min = 0.5)
#plot_missval(se)
#plot_detect(se)
plot_coverage(filt)
norm <- normalize_vsn(filt)
meanSdPlot(norm)
plot_numbers(norm)

#Imputation and DEP analysis

bpca <- impute(norm, fun = "bpca")
diff <- test_diff(bpca, type = "all")

imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
diff <- test_diff(imputed_manual, type = "all")

dep <- add_rejections(diff, alpha = 0.05, lfc = 0.264)

# Plot PCA
plot_pca(dep, indicate = c("condition","replicate"))
plot_volcano(dep, contrast = "Virus_CH_vs_Mock_CH")

data_results <- get_results(dep)
data_results <- data_results[, c(1,7,3)]
sig <- data_results[data_results$Virus_CH_vs_Mock_CH_p.val < 0.05,]
EnhancedVolcano(data_results,
                lab = data_results$name,
                selectLab = c('SPIKE'),
                x='Mock_CH_vs_Virus_CH_ratio',
                y='Mock_CH_vs_Virus_CH_p.val',
                pCutoff = 0.05,
                FCcutoff = 0.264,
                xlim = c(-4, 4),
                ylim = c(0,4),
                col=c('darkgrey', 'darkgrey', 'darkgrey', '#00B0F6'),
                colAlpha = 1,
                legendPosition = 'right'
) #export pdf as 8x6

write.table(data_results, file = "ChildNasal_MockvVirus_DEP.txt", sep = "\t", quote = F, row.names = F)
write.table(data_results, file = "Perseus-like/CH_MockvVirus_PL.txt", sep = "\t", quote = F, row.names = F)


#### Young Adult Mock vs. Virus ####

data <- read.table(file = "proteinGroups.txt", sep = "\t", header = TRUE)

data <- data[data$Reverse != "+" & data$Potential.contaminant != "+" & data$Unique.peptides >= 2,]
data <- data[data$Gene.names != "",] #removes rows with no Gene names
data <- data[complete.cases(data),] #removes NAs at the end 
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";") 

#Make summarised experiment
columns <- grep("LFQ.intensity.", colnames(data_unique))
exp_design <- read.table(file = "experimental_design_YA.txt", sep = "\t", header = TRUE) #Import experimental design - the same as LFQ-Analyst
se <- make_se(data_unique, columns, exp_design)

# Filter and normalise
filt <- filter_proteins(se, "fraction", min = 0.5)
#plot_missval(se)
#plot_detect(se)
plot_coverage(filt)
norm <- normalize_vsn(filt)
meanSdPlot(norm)
plot_numbers(norm)

#Imputation and DEP analysis

bpca <- impute(norm, fun = "bpca")
diff <- test_diff(bpca, type = "all")

imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
diff <- test_diff(imputed_manual, type = "all")

dep <- add_rejections(diff, alpha = 0.05, lfc = 0.264)

# Plot PCA
plot_pca(dep, indicate = c("condition","replicate"))
plot_volcano(dep, contrast = "Mock_vs_Virus")

data_results <- get_results(dep)
data_results <- data_results[, c(1,7,3)]
sig <- data_results[data_results$Mock_vs_Virus_p.val < 0.05, ] 
EnhancedVolcano(data_results,
                lab = data_results$name,
                selectLab = c('SPIKE'),
                x='Mock_vs_Virus_ratio',
                y='Mock_vs_Virus_p.val',
                pCutoff = 0.05,
                FCcutoff = 0.264,
                xlim = c(-4, 4),
                ylim = c(0,4),
                col=c('darkgrey', 'darkgrey', 'darkgrey', '#00B0F6'),
                colAlpha = 1,
                legendPosition = 'right'
) #export pdf as 8x6

write.table(data_results, file = "YA_MockvVirus_DEP.txt", sep = "\t", quote = F, row.names = F)

#### Older Adult Mock vs. Virus ####

data <- read.table(file = "proteinGroups.txt", sep = "\t", header = TRUE)

data <- data[data$Reverse != "+" & data$Potential.contaminant != "+" & data$Unique.peptides >= 2,]
data <- data[data$Gene.names != "",] #removes rows with no Gene names
data <- data[complete.cases(data),] #removes NAs at the end 
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";") 

#Make summarised experiment
columns <- grep("LFQ.intensity.", colnames(data_unique))
exp_design <- read.table(file = "experimental_design_OA.txt", sep = "\t", header = TRUE) #Import experimental design - the same as LFQ-Analyst
se <- make_se(data_unique, columns, exp_design)

# Filter and normalise
filt <- filter_proteins(se, "fraction", min = 0.5)
#plot_missval(se)
#plot_detect(se)
plot_coverage(filt)
norm <- normalize_vsn(filt)
meanSdPlot(norm)
plot_numbers(norm)

#Imputation and DEP analysis

bpca <- impute(norm, fun = "bpca")
diff <- test_diff(bpca, type = "all")

imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
diff <- test_diff(imputed_manual, type = "all")

dep <- add_rejections(diff, alpha = 0.05, lfc = 0.264)

# Plot PCA
plot_pca(dep, indicate = c("condition","replicate"))
plot_volcano(dep, contrast = "Mock_OA_vs_Virus_OA")

data_results <- get_results(dep)
data_results <- data_results[, c(1,7,3)]
sig <- data_results[data_results$Mock_OA_vs_Virus_OA_p.val < 0.05, ]

EnhancedVolcano(data_results,
                lab = data_results$name,
                selectLab = c('SPIKE'),
                x='Mock_OA_vs_Virus_OA_ratio',
                y='Mock_OA_vs_Virus_OA_p.val',
                pCutoff = 0.05,
                FCcutoff = 0.264,
                xlim = c(-4, 4),
                ylim = c(0,4),
                col=c('darkgrey', 'darkgrey', 'darkgrey', '#00B0F6'),
                colAlpha = 1,
                legendPosition = 'right'
) #export pdf as 8x6

write.table(data_results, file = "OA_MockvVirus_DEP.txt", sep = "\t", quote = F, row.names = F)



#### Child vs. YA vs. OA MOCK####

data <- read.table(file = "proteinGroups.txt", sep = "\t", header = TRUE)

data <- data[data$Reverse != "+" & data$Potential.contaminant != "+" & data$Unique.peptides >= 2,]
data <- data[data$Gene.names != "",] #removes rows with no Gene names
data <- data[complete.cases(data),] #removes NAs at the end 
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

#Make summarised experiment
columns <- grep("LFQ.intensity.", colnames(data_unique))
#Only HC and nasal samples. JW_Mock removed
exp_design <- read.table(file = "experimental_design_CHvAll.txt", sep = "\t", header = TRUE) #Import experimental design - the same as LFQ-Analyst
se <- make_se(data_unique, columns, exp_design)

# Filter and normalise
filt <- filter_proteins(se, "fraction", min = 0.5)
#plot_missval(se)
#plot_detect(se)
plot_coverage(filt)
norm <- normalize_vsn(filt)
meanSdPlot(norm)
plot_numbers(norm)

#Imputation and DEP analysis

bpca <- impute(norm, fun = "bpca")
diff <- test_diff(bpca, type = "all")

imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
diff <- test_diff(imputed_manual, type = "all")

dep <- add_rejections(diff, lfc = 0.264)

# Plot PCA
plot_pca(dep, indicate = c("condition","replicate"))

#Plot Heatmap
plot_heatmap(dep)
plot_volcano(dep, contrast = "CH_vs_All")

data_results <- get_results(dep)
data_results <- data_results[, c(1,14,4)]
sig <- data_results[data_results$CH_vs_All_p.val < 0.05, ]

EnhancedVolcano(data_results,
                lab = data_results$name,
                selectLab = c('SPIKE'),
                x='CH_vs_All_ratio',
                y='CH_vs_All_p.val',
                pCutoff = 0.05,
                FCcutoff = 0.264,
                xlim = c(-4, 4),
                ylim = c(0,4),
                col=c('darkgrey', 'darkgrey', 'blue', '#00B0F6'),
                colAlpha = 1,
                legendPosition = 'right'
)

write.table(data_results, file = "ChildvAll_Mock_DEP.txt", sep = "\t", quote = F, row.names = F)
write.table(data_results, file = "Perseus-like/ChildvAll_Mock_PL.txt", sep = "\t", quote = F, row.names = F)

#### YA vs. CH + OA MOCK ####

data <- read.table(file = "proteinGroups.txt", sep = "\t", header = TRUE)

data <- data[data$Reverse != "+" & data$Potential.contaminant != "+" & data$Unique.peptides >= 2,]
data <- data[data$Gene.names != "",] #removes rows with no Gene names
data <- data[complete.cases(data),] #removes NAs at the end 
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";") 

#Make summarised experiment
columns <- grep("LFQ.intensity.", colnames(data_unique))
#Only HC and nasal samples. JW_Mock removed
exp_design <- read.table(file = "experimental_design_YAvAll.txt", sep = "\t", header = TRUE) #Import experimental design - the same as LFQ-Analyst
se <- make_se(data_unique, columns, exp_design)

# Filter and normalise
filt <- filter_proteins(se, "fraction", min = 0.5)
#plot_missval(se)
#plot_detect(se)
plot_coverage(filt)
norm <- normalize_vsn(filt)
meanSdPlot(norm)
plot_numbers(norm)

#Imputation and DEP analysis

bpca <- impute(norm, fun = "bpca")
diff <- test_diff(bpca, type = "all")

imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
diff <- test_diff(imputed_manual, type = "all")

dep <- add_rejections(diff, lfc = 0.264)

# Plot PCA
plot_pca(dep, indicate = c("condition","replicate"))
plot_volcano(dep, contrast = "YA_vs_All")

data_results <- get_results(dep)
data_results <- data_results[, c(1,14,4)]
sig <- data_results[data_results$YA_vs_All_p.val < 0.05, ]

EnhancedVolcano(data_results,
                lab = data_results$name,
                selectLab = c('SPIKE'),
                x='YA_vs_All_ratio',
                y='YA_vs_All_p.val',
                pCutoff = 0.05,
                FCcutoff = 0.264,
                xlim = c(-4, 4),
                ylim = c(0,4),
                col=c('darkgrey', 'darkgrey', 'blue', '#00B0F6'),
                colAlpha = 1,
                legendPosition = 'right'
                )

write.table(data_results, file = "YAvAll_Mock_DEP.txt", sep = "\t", quote = F, row.names = F)
write.table(data_results, file = "Perseus-like/YAvAll_Mock_PL.txt", sep = "\t", quote = F, row.names = F)


#### OA vs. CH + YA MOCK####

data <- read.table(file = "proteinGroups.txt", sep = "\t", header = TRUE)

data <- data[data$Reverse != "+" & data$Potential.contaminant != "+" & data$Unique.peptides >= 2,]
data <- data[data$Gene.names != "",] #removes rows with no Gene names
data <- data[complete.cases(data),] #removes NAs at the end 
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";") 

#Make summarised experiment
columns <- grep("LFQ.intensity.", colnames(data_unique))
#Only HC and nasal samples. JW_Mock removed
exp_design <- read.table(file = "experimental_design_OAvAll.txt", sep = "\t", header = TRUE) #Import experimental design - the same as LFQ-Analyst
se <- make_se(data_unique, columns, exp_design)

# Filter and normalise
filt <- filter_proteins(se, "fraction", min = 0.5)
#plot_missval(se)
#plot_detect(se)
plot_coverage(filt)
norm <- normalize_vsn(filt)
meanSdPlot(norm)
plot_numbers(norm)

#Imputation and DEP analysis

bpca <- impute(norm, fun = "bpca")
diff <- test_diff(bpca, type = "all")

imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
diff <- test_diff(imputed_manual, type = "all")

dep <- add_rejections(diff, lfc = 0.264)

# Plot PCA
plot_pca(dep, indicate = c("condition","replicate"))
plot_volcano(dep, contrast = "OA_vs_All")

data_results <- get_results(dep)
data_results <- data_results[, c(1,14,4)]
sig <- data_results[data_results$OA_vs_All_p.val < 0.05,]

EnhancedVolcano(data_results,
                lab = data_results$name,
                selectLab = c('SPIKE'),
                x='OA_vs_All_ratio',
                y='OA_vs_All_p.val',
                pCutoff = 0.05,
                FCcutoff = 0.264,
                xlim = c(-4, 4),
                ylim = c(0,4),
                col=c('darkgrey', 'darkgrey', 'blue', '#00B0F6'),
                colAlpha = 1,
                legendPosition = 'right'
)

write.table(data_results, file = "OAvAll_Mock_DEP.txt", sep = "\t", quote = F, row.names = F)
write.table(data_results, file = "OAvAll_Mock_PL.txt", sep = "\t", quote = F, row.names = F)

#### Mock vs. Virus All ####

data <- read.table(file = "proteinGroups.txt", sep = "\t", header = TRUE)

data <- data[data$Reverse != "+" & data$Potential.contaminant != "+" & data$Unique.peptides >= 2,]
data <- data[data$Gene.names != "",] #removes rows with no Gene names
data <- data[complete.cases(data),] #removes NAs at the end 
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";") 

#Make summarised experiment
columns <- grep("LFQ.intensity.", colnames(data_unique))
exp_design <- read.table(file = "experimental_design_MockvVirus_All.txt", sep = "\t", header = TRUE) #Import experimental design - the same as LFQ-Analyst
se <- make_se(data_unique, columns, exp_design)

# Filter and normalise
filt <- filter_proteins(se, "fraction", min = 0.5)
#plot_missval(se)
#plot_detect(se)
plot_coverage(filt)
norm <- normalize_vsn(filt)
meanSdPlot(norm)
plot_numbers(norm)

#Imputation and DEP analysis

bpca <- impute(norm, fun = "bpca")
diff <- test_diff(bpca, type = "all")

imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
diff <- test_diff(imputed_manual, type = "all")

dep <- add_rejections(diff, alpha = 0.05, lfc = 0.264)

# Plot PCA
plot_pca(dep, indicate = c("condition","replicate"))

data_results <- get_results(dep)
data_results <- data_results[, c(1,7,3)]
sig <- data_results[data_results$Virus_vs_Mock_p.val < 0.05,]

EnhancedVolcano(data_results,
                lab = data_results$name,
                selectLab = c('SPIKE'),
                x='Virus_vs_Mock_ratio',
                y='Virus_vs_Mock_p.val',
                pCutoff = 0.05,
                FCcutoff = 0.264,
                xlim = c(-4, 4),
                ylim = c(0,4),
                col=c('darkgrey', 'darkgrey', 'darkgrey', '#00B0F6'),
                colAlpha = 1,
                legendPosition = 'right'
) #7x5.5 export

sig_results <- data_results[data_results$Mock_vs_Virus_p.val < 0.05, ]

write.table(data_results, file = "MockvVirus_DEP.txt", sep = "\t", quote = F, row.names = F)
write.table(data_results, file = "Perseus-like/MockvVirus_PL.txt", sep = "\t", quote = F, row.names = F)


#### Child vs. YA vs. OA VIRUS ####

data <- read.table(file = "proteinGroups.txt", sep = "\t", header = TRUE)

data <- data[data$Reverse != "+" & data$Potential.contaminant != "+" & data$Unique.peptides >= 2,]
data <- data[data$Gene.names != "",] #removes rows with no Gene names
data <- data[complete.cases(data),] #removes NAs at the end 
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

#Make summarised experiment
columns <- grep("LFQ.intensity.", colnames(data_unique))
#Only HC and nasal samples. JW_Mock removed
exp_design <- read.table(file = "experimental_design_CHvAll_Virus.txt", sep = "\t", header = TRUE) #Import experimental design - the same as LFQ-Analyst
se <- make_se(data_unique, columns, exp_design)

# Filter and normalise
filt <- filter_proteins(se, "fraction", min = 0.5)
#plot_missval(se)
#plot_detect(se)
plot_coverage(filt)
norm <- normalize_vsn(filt)
meanSdPlot(norm)
plot_numbers(norm)

#Imputation and DEP analysis

bpca <- impute(norm, fun = "bpca")
diff <- test_diff(bpca, type = "all")

imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
diff <- test_diff(imputed_manual, type = "all")

dep <- add_rejections(diff, lfc = 0.264)

# Plot PCA
plot_pca(dep, indicate = c("condition","replicate"))

#Plot Heatmap
plot_heatmap(dep)
plot_volcano(dep, contrast = "CH_vs_All")

data_results <- get_results(dep)
data_results <- data_results[, c(1,14,4)]
sig <- data_results[data_results$CH_vs_All_p.val < 0.05, ]

EnhancedVolcano(data_results,
                lab = data_results$name,
                selectLab = c('SPIKE'),
                x='CH_vs_All_ratio',
                y='CH_vs_All_p.val',
                pCutoff = 0.05,
                FCcutoff = 0.264,
                xlim = c(-4, 4),
                ylim = c(0,4),
                col=c('darkgrey', 'darkgrey', 'blue', '#00B0F6'),
                colAlpha = 1,
                legendPosition = 'right'
)

write.table(data_results, file = "ChildvAll_Virus_DEP.txt", sep = "\t", quote = F, row.names = F)
write.table(data_results, file = "Perseus-like/ChildvAll_Mock_PL.txt", sep = "\t", quote = F, row.names = F)



#### YA vs. CH + OA VIRUS ####

data <- read.table(file = "proteinGroups.txt", sep = "\t", header = TRUE)

data <- data[data$Reverse != "+" & data$Potential.contaminant != "+" & data$Unique.peptides >= 2,]
data <- data[data$Gene.names != "",] #removes rows with no Gene names
data <- data[complete.cases(data),] #removes NAs at the end 
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";") 

#Make summarised experiment
columns <- grep("LFQ.intensity.", colnames(data_unique))
#Only HC and nasal samples. JW_Mock removed
exp_design <- read.table(file = "experimental_design_YAvAll_Virus.txt", sep = "\t", header = TRUE) #Import experimental design - the same as LFQ-Analyst
se <- make_se(data_unique, columns, exp_design)

# Filter and normalise
filt <- filter_proteins(se, "fraction", min = 0.5)
#plot_missval(se)
#plot_detect(se)
plot_coverage(filt)
norm <- normalize_vsn(filt)
meanSdPlot(norm)
plot_numbers(norm)

#Imputation and DEP analysis

bpca <- impute(norm, fun = "bpca")
diff <- test_diff(bpca, type = "all")

imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
diff <- test_diff(imputed_manual, type = "all")

dep <- add_rejections(diff, lfc = 0.264)

# Plot PCA
plot_pca(dep, indicate = c("condition","replicate"))
plot_volcano(dep, contrast = "YA_vs_All")

data_results <- get_results(dep)
data_results <- data_results[, c(1,14,4)]
sig <- data_results[data_results$YA_vs_All_p.val < 0.05, ]

EnhancedVolcano(data_results,
                lab = data_results$name,
                selectLab = c('SPIKE'),
                x='YA_vs_All_ratio',
                y='YA_vs_All_p.val',
                pCutoff = 0.05,
                FCcutoff = 0.264,
                xlim = c(-4, 4),
                ylim = c(0,4),
                col=c('darkgrey', 'darkgrey', 'blue', '#00B0F6'),
                colAlpha = 1,
                legendPosition = 'right'
)

write.table(data_results, file = "YAvAll_Virus_DEP.txt", sep = "\t", quote = F, row.names = F)
write.table(data_results, file = "Perseus-like/YAvAll_Mock_PL.txt", sep = "\t", quote = F, row.names = F)



#### OA vs. CH + OA VIRUS ####

data <- read.table(file = "proteinGroups.txt", sep = "\t", header = TRUE)

data <- data[data$Reverse != "+" & data$Potential.contaminant != "+" & data$Unique.peptides >= 2,]
data <- data[data$Gene.names != "",] #removes rows with no Gene names
data <- data[complete.cases(data),] #removes NAs at the end 
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";") 

#Make summarised experiment
columns <- grep("LFQ.intensity.", colnames(data_unique))
#Only HC and nasal samples. JW_Mock removed
exp_design <- read.table(file = "experimental_design_OAvAll.txt", sep = "\t", header = TRUE) #Import experimental design - the same as LFQ-Analyst
se <- make_se(data_unique, columns, exp_design)

# Filter and normalise
filt <- filter_proteins(se, "fraction", min = 0.5)
#plot_missval(se)
#plot_detect(se)
plot_coverage(filt)
norm <- normalize_vsn(filt)
meanSdPlot(norm)
plot_numbers(norm)

#Imputation and DEP analysis

bpca <- impute(norm, fun = "bpca")
diff <- test_diff(bpca, type = "all")

imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
diff <- test_diff(imputed_manual, type = "all")

dep <- add_rejections(diff, lfc = 0.264)

# Plot PCA
plot_pca(dep, indicate = c("condition","replicate"))
plot_volcano(dep, contrast = "OA_vs_All")

data_results <- get_results(dep)
data_results <- data_results[, c(1,14,4)]
sig <- data_results[data_results$OA_vs_All_p.val < 0.05,]

EnhancedVolcano(data_results,
                lab = data_results$name,
                selectLab = c('SPIKE'),
                x='OA_vs_All_ratio',
                y='OA_vs_All_p.val',
                pCutoff = 0.05,
                FCcutoff = 0.264,
                xlim = c(-4, 4),
                ylim = c(0,4),
                col=c('darkgrey', 'darkgrey', 'blue', '#00B0F6'),
                colAlpha = 1,
                legendPosition = 'right'
)

write.table(data_results, file = "OAvAll_Virus_DEP.txt", sep = "\t", quote = F, row.names = F)
write.table(data_results, file = "Perseus-like/OAvAll_Mock_PL.txt", sep = "\t", quote = F, row.names = F)


