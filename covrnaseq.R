library(tidyverse)
library(edgeR)
library(plotly)
library(gprofiler2)
library(gplots)

sampleinfo <- read.table("metadata.txt", header = T, sep = "\t")
##dashes replaced with dots and start digits removed
sampleinfo <- sampleinfo %>% mutate(sname = gsub("^\\d+-", "", sampleid, perl = T)) %>% mutate(sname = gsub("-",".",sname))
sampleinfo <- sampleinfo %>% mutate(fname = gsub("^\\d+-", "", filebase, perl = T)) %>% mutate(fname = gsub("-",".",fname)) %>% mutate(fname = gsub("_R1$", "", fname, perl = T))

countsummary <- read.table("covseq.counts.summary", header = T, sep = "\t", row.names = 1)
colnames(countsummary) <- gsub(".bam", "", colnames(countsummary))
countsummary <- data.frame(t(countsummary))
countsummary <- countsummary[,colnames(countsummary)[colSums(countsummary) > 0]]
countsummary$fname <- row.names(countsummary)
countsummary <- mutate(countsummary, fname = gsub("_R1$", "", fname, perl = T)) %>% mutate(fname=gsub("^X\\d+\\.", "", fname, perl = T))

sampleinfo <- left_join(sampleinfo, countsummary, by = "fname")

##figure for proportion of reads mapped and (un)assigned
sampleinfo %>% select(fname,Assigned,Unassigned_Unmapped,Unassigned_NoFeatures) %>% gather("feature", "counts", -fname) %>% mutate(fname = as.factor(fname), feature = as.factor(feature)) %>% left_join(sampleinfo, by = "fname") %>% mutate(proportion = counts/readcounts) %>% select(fname,feature,proportion) %>% plot_ly(x=~fname, y=~proportion, color = ~feature, type="bar") %>% layout(barmode = "stack")

##read counts file
counts <- read.table("covseq.counts", header = T, sep = "\t")
colnames(counts) <- gsub("_R1.bam$", "", colnames(counts), perl = T)
colnames(counts) <- gsub("^X\\d+\\.", "", colnames(counts), perl = T)

lcounts <- counts %>% select(!(Chr:Length)) %>% filter(rowMeans(select(.,-Geneid))>=5)
row.names(lcounts) <- lcounts$Geneid
lcounts <- select(lcounts, -Geneid)
##this is important to do
##arrange sampleinfo table to match the order of colnames in lcounts
##this will ensure that you can use the table as is
sampleinfo <- sampleinfo[match(colnames(lcounts), sampleinfo$fname),]

dge <- DGEList(counts = as.matrix(lcounts), lib.size = sampleinfo$readcounts, group = sampleinfo$laneid)
dge <- calcNormFactors(dge)
design <- model.matrix(~ 0 + sampleinfo$laneid + sampleinfo$treatment + sampleinfo$technicalrepid)
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = 5)
diffex <- topTags(lrt, n = nrow(lcounts))
dx <- diffex$table
dx$gid <- row.names(dx)
positivefc <- gost(dx %>% filter(FDR <= 0.05 & logFC > 0) %>% .$gid, organism = "hsapiens", custom_bg = dx$gid)
negativefc <- gost(dx %>% filter(FDR <= 0.05 & logFC < 0) %>% .$gid, organism = "hsapiens", custom_bg = dx$gid)


if (! is.null(positivefc)){
  write.table(positivefc$result %>% select(-parents), file = "positiveFC.functionalanalysis.txt", col.names = T, row.names = F, sep = "\t", quote = F)
}
if (! is.null(negativefc)){
  write.table(negativefc$result %>% select(-parents), file = "negativeFC.functionalanalysis.txt", col.names = T, row.names = F, sep = "\t", quote = F)
}
  
write.table(dx, file = "diffex.txt", col.names = T, row.names = T, quote = F, sep = "\t")

plotMDS(dge, labels = sampleinfo[match(row.names(dge$samples), sampleinfo$fname),"treatment"])

x <- cpm.DGEList(dge, normalized.lib.sizes = T, log = T)
dxgenes <- x[dx %>% filter(FDR <= 0.05) %>% .$gid,]
hm <- heatmap.2(as.matrix(x[dx %>% filter(FDR <= 0.05) %>% .$gid,]))
heatmap <- plot_ly(z = hm$carpet, x = colnames(hm$carpet), y = row.names(hm$carpet), type = "heatmap")
htmlwidgets::saveWidget(as_widget(heatmap), "diffexheatmap.html")








               
               