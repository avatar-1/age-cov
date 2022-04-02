# This code is used to analyse lysate proteomics for the COVID-19 study
# Code developed by Alexander Capraro August 2021
library(ggplot2)
library(gplots)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(ggfortify)
library(DEP)
library(factoextra)
library(UpSetR)
library(EnhancedVolcano)
library(magrittr)
library(data.table)
library(SummarizedExperiment)
options(scipen=999) #removes scientific notation

# Read in files ####

data_unique <- read.table(file = "proteinGroups_unique.txt", sep = "\t", header = TRUE)
sinfo <- read.table(file= "sampleinfo.txt", sep="\t", header=T)
log2lfq <- read.table(file="log2_lfq.txt", sep = "\t", header=T)

# Create LFG intensity file ####
data <- read.table(file="proteinGroups.txt",sep="\t",header=TRUE)

sinfo <- read.table(file= "sampleinfo.txt", sep="\t", header=T)

data <- data[data$Reverse != "+" & data$Potential.contaminant != "+" & data$Unique.peptides >= 2,]
data <- data[data$Gene.names != "",] #removes rows with no Gene names
data <- data[complete.cases(data),] #removes NAs at the end 
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

lfq_intensity <- data_unique %>%
  set_rownames(.$Gene.names) %>%
  select(contains("LFQ.intensity.")) %>%
  log2()
lfq_intensity[lfq_intensity == "-Inf"] <- 0

colnames(lfq_intensity) <- sub("LFQ.intensity.","", colnames(lfq_intensity)) #subsitutes colnames with LFQ.intensity removed
mm <- match(names(lfq_intensity), sinfo$label) #match names with label
names(lfq_intensity)[!is.na(mm)] <- as.character(sinfo$sname[na.omit(mm)])

lfq_intensity <- lfq_intensity %>%
  select(-contains("glufib-")) %>% #remove glufib
  select(-contains(".NA.")) #remove non-SARS samples

write.table(lfq_intensity,file = "log2_lfq.txt",sep="\t",quote=F)
write.table(data_unique,file="proteinGroups_unique.txt", sep="\t", quote=F)

# Differential protein expression analysis ####

make_se <- function(proteins_unique, columns, expdesign) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(proteins_unique),
                          is.integer(columns),
                          is.data.frame(expdesign))
  
  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)),
         "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns",
         "are not present in the experimental design",
         call. = FALSE)
  }
  if(any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric",
         "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }
  
  # If input is a tibble, convert to data.frame
  if(tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)
  if(tibble::is_tibble(expdesign))
    expdesign <- as.data.frame(expdesign)
  
  # Select the assay data
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)
  
  # Generate the colData from the experimental design
  # and match these with the assay data
  expdesign <- mutate(expdesign, condition = make.names(condition)) %>%
    unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID
  
  matched <- match(make.names(expdesign$label),
                   make.names(colnames(raw)))
  if(any(is.na(matched))) {
    stop("None of the labels in the experimental design match ",
         "with column names in 'proteins_unique'",
         "\nRun make_se() with the correct labels in the experimental design",
         "and/or correct columns specification")
  }
  
  colnames(raw)[matched] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
  
  # Select the rowData
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name
  
  # Generate the SummarizedExperiment object
  se <- SummarizedExperiment(assays = as.matrix(raw),
                             colData = expdesign,
                             rowData = row_data)
  return(se)
} #Slightly edited from DEP as it wasn't reading names correctly
diffex <- function(cond1,cond2) {
  exp_design <- sinfo %>%
    select(label,condition,replicate) %>%
    filter(condition == cond1 | condition == cond2) %>%
    mutate(label=as.character(label))
  
  exp_design$label <- paste("LFQ.intensity.",exp_design$label,sep="")
  matched <- match(make.names(exp_design$label),
                   make.names(colnames(data_unique)))
  #exp_design$label <- gsub('\\.', '-', exp_design$label)

  columns <- grep("LFQ.intensity.", colnames(data_unique))
  se <- make_se(data_unique, columns, exp_design)
  
  # Filter and normalise
  filt <- filter_proteins(se, "fraction", min = 0.5)
  #plot_missval(se)
  #plot_detect(se)
  #plot_coverage(filt)
  norm <- normalize_vsn(filt)
  #meanSdPlot(norm)
  #plot_numbers(norm)
  
  #Imputation and DEP analysis
  
  imputed <- impute(norm, fun = "bpca")
  diff <- test_diff(imputed, type = "all", test = cond1, control = cond2)
  dep <- add_rejections(diff, lfc = 0.264)
  data_results <- get_results(dep)
  data_results <- data_results[, c(1,7,3,4)]
  names(data_results) <- c("Protein","logFC","Pvalue","P.adj")
  sig <- data_results %>% filter(Pvalue <= 0.05, logFC > 0.264 | logFC < -0.264)
  print(sig)
  print(paste("Upregulated:", nrow(sig[sig$logFC > 0,]), sep=" "))
  print(paste("Downregulated:", nrow(sig[sig$logFC < 0,]), sep=" "))
  
  comp <- paste(cond1,"_vs_",cond2,sep = "")
  
  keyvals <- ifelse(
    data_results$logFC <= -0.264 & data_results$Pvalue < 0.05, '#E8495C',
    ifelse(data_results$logFC >= 0.264 & data_results$Pvalue < 0.05, '#38ACE2',
           'grey'))
  names(keyvals)[keyvals == '#38ACE2'] <- 'Upregulated'
  names(keyvals)[keyvals == 'grey'] <- 'Not significant'
  names(keyvals)[keyvals == '#E8495C'] <- 'Downregulated'
  
  volcanoplot <- EnhancedVolcano(data_results,
                                 lab = row.names(data_results),
                                 selectLab = c('SPIKE'),
                                 x='logFC',
                                 y='Pvalue',
                                 pCutoff = 0.05,
                                 FCcutoff = 0.264,
                                 colCustom = keyvals,
                                 colAlpha = 1,
                                 title = paste(cond1,"vs",cond2, sep=" "),
                                 subtitle = "",
                                 legendPosition = 'none')
  plot(volcanoplot)

  write.table(sig,file = paste0("DEP.",cond1,".vs.",cond2,".txt"), sep = "\t", quote = F, row.names = F)
} #Cond1 = numerator, cond2 = denominator

# Virus 0.2 vs. Mock
diffex(cond1 = "V2.CH.N", cond2 = "M.CH.N")
diffex(cond1 = "V2.YA", cond2 = "M.YA")
diffex(cond1 = "V2.OA", cond2 = "M.OA")
diffex(cond1 = "V2.CH.B", cond2 = "M.CH.B")

#Virus 0.6 vs. Virus 0.2
diffex(cond1 = "V6.CH", cond2 = "V2.CH.N")
diffex(cond1 = "V6.YA", cond2 = "V2.YA")
diffex(cond1 = "V6.OA", cond2 = "V2.OA")

#Virus 0.6 vs. Virus 0.2
diffex(cond1 = "Nafa2.CH", cond2 = "V2.CH.N")
diffex(cond1 = "Nafa2.YA", cond2 = "V2.YA")
diffex(cond1 = "Nafa2.OA", cond2 = "V2.OA")

#Copper 0.2 vs. Virus 0.2
diffex(cond1 = "Cu2.CH", cond2 = "V2.CH.N")
diffex(cond1 = "Cu2.YA", cond2 = "V2.YA")
diffex(cond1 = "Cu2.OA", cond2 = "V2.OA")

#Copper 0.6 vs. Virus 0.6
diffex(cond1 = "Cu6.CH", cond2 = "V6.CH")
diffex(cond1 = "Cu6.YA", cond2 = "V6.YA")
diffex(cond1 = "Cu6.OA", cond2 = "V6.OA")

#Copper Mock vs. Mock
diffex(cond1 = "CuM.CH", cond2 = "M.CH.N")
diffex(cond1 = "CuM.YA", cond2 = "M.YA")
diffex(cond1 = "CuM.OA", cond2 = "M.OA")

#Mock Age vs. Age
diffex(cond1 = "M.CH.N", cond2 = "M.YA")
diffex(cond1 = "M.CH.N", cond2 = "M.OA")
diffex(cond1 = "M.OA", cond2 = "M.YA")

#Nafa 0.2 Age vs. Age
diffex(cond1 = "Nafa2.CH", cond2 = "Nafa2.YA")
diffex(cond1 = "Nafa2.CH", cond2 = "Nafa2.OA")
diffex(cond1 = "Nafa2.OA", cond2 = "Nafa2.YA")

# Differential protein expression analysis - ALL AGES ####

make_se <- function(proteins_unique, columns, expdesign) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(proteins_unique),
                          is.integer(columns),
                          is.data.frame(expdesign))
  
  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)),
         "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns",
         "are not present in the experimental design",
         call. = FALSE)
  }
  if(any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric",
         "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }
  
  # If input is a tibble, convert to data.frame
  if(tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)
  if(tibble::is_tibble(expdesign))
    expdesign <- as.data.frame(expdesign)
  
  # Select the assay data
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)
  
  # Generate the colData from the experimental design
  # and match these with the assay data
  expdesign <- mutate(expdesign, condition = make.names(condition)) %>%
    unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID
  
  matched <- match(make.names(expdesign$label),
                   make.names(colnames(raw)))
  if(any(is.na(matched))) {
    stop("None of the labels in the experimental design match ",
         "with column names in 'proteins_unique'",
         "\nRun make_se() with the correct labels in the experimental design",
         "and/or correct columns specification")
  }
  
  colnames(raw)[matched] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
  
  # Select the rowData
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name
  
  # Generate the SummarizedExperiment object
  se <- SummarizedExperiment(assays = as.matrix(raw),
                             colData = expdesign,
                             rowData = row_data)
  return(se)
} #Slightly edited from DEP as it wasn't reading names correctly
diffex <- function(cond1,cond2) {
  exp_design <- sinfo %>%
    select(label,all_condition,all_replicate) %>%
    filter(all_condition == cond1 | all_condition == cond2) %>%
    mutate(label=as.character(label))
  
  exp_design$label <- paste("LFQ.intensity.",exp_design$label,sep="")
  names(exp_design) <- c("label","condition","replicate")
  matched <- match(make.names(exp_design$label),
                   make.names(colnames(data_unique)))
  #exp_design$label <- gsub('\\.', '-', exp_design$label)
  
  columns <- grep("LFQ.intensity.", colnames(data_unique))
  se <- make_se(data_unique, columns, exp_design)
  
  # Filter and normalise
  filt <- filter_proteins(se, "fraction", min = 0.5)
  #plot_missval(se)
  #plot_detect(se)
  #plot_coverage(filt)
  norm <- normalize_vsn(filt)
  #meanSdPlot(norm)
  #plot_numbers(norm)
  
  #Imputation and DEP analysis
  
  imputed <- impute(norm, fun = "bpca")
  diff <- test_diff(imputed, type = "all", test = cond1, control = cond2)
  dep <- add_rejections(diff, lfc = 0.264)
  data_results <- get_results(dep)
  data_results <- data_results[, c(1,7,3,4)]
  names(data_results) <- c("Protein","logFC","Pvalue","P.adj")
  sig <- data_results %>% filter(Pvalue <= 0.05, logFC > 0.264 | logFC < -0.264)
  print(sig)
  print(paste("Upregulated:", nrow(sig[sig$logFC > 0,]), sep=" "))
  print(paste("Downregulated:", nrow(sig[sig$logFC < 0,]), sep=" "))
  
  comp <- paste(cond1,"_vs_",cond2,sep = "")
  
  keyvals <- ifelse(
    data_results$logFC <= -0.264 & data_results$Pvalue < 0.05, '#E8495C',
    ifelse(data_results$logFC >= 0.264 & data_results$Pvalue < 0.05, '#38ACE2',
           'grey'))
  names(keyvals)[keyvals == '#38ACE2'] <- 'Upregulated'
  names(keyvals)[keyvals == 'grey'] <- 'Not significant'
  names(keyvals)[keyvals == '#E8495C'] <- 'Downregulated'
  
  volcanoplot <- EnhancedVolcano(data_results,
                                 lab = row.names(data_results),
                                 selectLab = c('SPIKE'),
                                 x='logFC',
                                 y='Pvalue',
                                 pCutoff = 0.05,
                                 FCcutoff = 0.264,
                                 colCustom = keyvals,
                                 colAlpha = 1,
                                 title = paste(cond1,"vs",cond2, sep=" "),
                                 subtitle = "",
                                 legendPosition = 'none')
  plot(volcanoplot)
  
  write.table(sig,file = paste0("DEP.",cond1,".vs.",cond2,".txt"), sep = "\t", quote = F, row.names = F)
} #Cond1 = numerator, cond2 = denominator

diffex(cond1="CuM",cond2="Mock")
diffex(cond1="Cu2",cond2="V2")
diffex(cond1="Cu6",cond2="V6")
diffex(cond1="Nafa2",cond2="V2")
diffex(cond1="V6",cond2="V2")
