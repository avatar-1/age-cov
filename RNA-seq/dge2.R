.libPaths( c( "/g/data/te53/software/Rpackages/4.0.0", .libPaths()))
library(tidyverse)
library(edgeR)
library(plotly)
library(gprofiler2)
library(RUVSeq)
library(Hmisc)
library(pvclust)
library(ggpubr)
library(ggrepel)


sampleinfo <- read_delim("metadata//sampleinfo.txt", delim="\t")
countsummary <- read_delim("./covseqcounts.2.txt.summary", delim="\t")
counts <- read_delim("covseqcounts.2.txt", delim = "\t", comment = "#")
viruscounts <- read_delim("NC_045512.2.counts.txt", delim=",", col_names = F)

# Sanitise headers to remove file paths.
# Organise countsummary information to the long format for merging with sampleinfo.
# Remove CF treatment samples from all files
# Create gene info table from counts table
# Remove gene info (chr, start, end, strand, length) from the counts table

colnames(countsummary)[-1] <- str_match(colnames(countsummary)[-1], "\\S+\\/(\\S+)\\.bam")[,2]
colnames(counts)[-c(1:6)] <- str_match(colnames(counts)[-c(1:6)], "\\S+\\/(\\S+)\\.bam")[,2]
colnames(viruscounts) <- c("runid", "virusmapped")
countsummary <- countsummary %>% pivot_longer(!Status, names_to = "runid", values_to = "counts") %>% 
                pivot_wider(names_from = "Status", values_from = "counts")

sampleinfo <- left_join(sampleinfo, countsummary, by = "runid")
sampleinfo <- left_join(sampleinfo, viruscounts, by = "runid")
sampleinfo <- mutate(sampleinfo, mapped = readcount - Unassigned_Unmapped, humanmapped = mapped - virusmapped)
geneinfo <- select(counts,1:6)
counts <- select(counts,-c(2:6))

sampleinfo <- filter(sampleinfo, treatment != "CF")
counts <- select(counts, c(Geneid,sampleinfo$runid))

sinfo <- left_join(select(sampleinfo, 8:18) %>% distinct(), 
                   select(sampleinfo, c(8,19,22,35:37)) %>% group_by(libraryid) %>% summarise(across(readcount:humanmapped, sum)), by = "libraryid") %>%
            mutate(t = case_when(str_detect(treatment, "Virus6") ~ "V6", str_detect(treatment,"Virus") ~ "V2", str_detect(treatment, "Mock") ~ "M"),
            sname = paste(t,agegroup,sex,tissue,individualid,technicalrepid, sep = "-"))

sampleinfo <- left_join(sampleinfo, select(sinfo, libraryid, sname), by = "libraryid")

cinfo <- left_join(pivot_longer(counts, !Geneid, names_to = "runid", values_to = "counts"), select(sampleinfo, runid,sname), by = "runid") %>% 
            group_by(sname, Geneid) %>% 
            summarise(counts = sum(counts)) %>% 
            pivot_wider(Geneid, names_from = "sname", values_from = "counts")

##this step is extremely essential. Order of columns in counts table if same as they appear in the sinfo, then life becomes easy in pulling things out of the sinfo when needed
cinfo <- select(cinfo, Geneid, sinfo$sname)

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

## function to hard filter to retain genes expressed in topthreshold % in any one sample

seltop <- function (x,threshold=0.99) 
    {
        slist <- cumsum(x[order(x,decreasing=TRUE)])/sum(x)<=threshold;
        bottom<-labels(slist)[slist==FALSE];
        x[bottom]<-NA;
        return(x)
    }
      
drawmds <- function(rawcounts, sdetails, topthreshold = 0.995, rmvirus = T, plotlabel = "all") {
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
 
    dge <- DGEList(counts = as.matrix(fcounts), lib.size = sdetails$readcount, group = sdetails$treatment)
    dge <- calcNormFactors(dge)
    z <- plotMDS.DGEList(dge,top = 1000,plot = F)$cmdscale.out
    mdsval <- as.tibble(z) %>% mutate(sname = row.names(z))
    colnames(mdsval) <- c("leading logFC dim1", "leading logFC dim2", "sname")
    mdsval <- left_join(mdsval, sdetails, by = "sname")
    mdsval %>% ggplot() + 
    geom_point(aes(x=`leading logFC dim1`, y = `leading logFC dim2`, colour = treatment, shape = agegroup), size = 3) + 
    geom_label_repel(aes(x=`leading logFC dim1`, y = `leading logFC dim2`,label = sex, fill=experiment), alpha = 0.5, size = 2, segment.color = 'black') +
      theme_classic()

} 

diffx <- function(rawcounts, sdetails, topthreshold = 0.995, rmvirus = T, plotlabel = "all") {
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
 
    dge <- DGEList(counts = as.matrix(fcounts), lib.size = sdetails$readcount, group = sdetails$treatment)
    dge <- calcNormFactors(dge)
    tpm <- cpm.DGEList(dge, log = T)
  
    design <- model.matrix(~sdetails$individualid + sdetails$treatment)
    y <- estimateDisp(dge,design)
    fit <- glmFit(y, design)
    qlf <- glmLRT(fit)
    diffex <- topTags(qlf, n = nrow(fcounts))
    dx <- diffex$table
    dx$gid <- row.names(dx)
    dxgenes <- tpm[dx %>% filter(FDR <= 0.05) %>% .$gid,]
    write_delim(dx, file = paste("diffex_glm",plotlabel, topthreshold, "txt", sep ="."), delim="\t")
    heatmap.2(as.matrix(tpm[dx %>% filter(FDR <= 0.05) %>% .$gid,]), trace = 'none', density.info = 'none', col = viridis)
} 

sdetails <- sinfo
rawcounts <- select(cinfo, Geneid, sdetails$sname)
drawmds(rawcounts = rawcounts, sdetails = sdetails)

sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1"))
rawcounts <- select(cinfo, Geneid, sdetails$sname)
drawmds(rawcounts = rawcounts, sdetails = sdetails)

sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1"), tissue == "N")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
drawmds(rawcounts = rawcounts, sdetails = sdetails)

sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1"), tissue == "B")
rawcounts <- select(cinfo, Geneid, sdetails$sname)
drawmds(rawcounts = rawcounts, sdetails = sdetails)

#Batch effect was visible on the second dimension (B1,B2, B3 colored square boxes). There is a tendency of older adults to be different as well. However, this can be the biological signal that we are after. RUVseq can potentially eliminate batch effect. RUVseq offers few ways of removing unwanted variation. First is to use gene expression based information where spike-in control genes or empirically estimated non-differentially expressed genes can be used to remove noise. Second is to use replicate samples and third is to use residuals from design matrix to remove the noise. We will attempt the third method.

rmvirus <- TRUE
topthreshold <- 0.995
sdetails <- filter(sinfo, !(individualid == "BC" & technicalrepid == "t1"))
rawcounts <- select(cinfo, Geneid, sdetails$sname)
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
normset <- RUVr(fcounts, row.names(fcounts), k=1, res)
z <- plotMDS(log(normset$normalizedCounts),top = 1000,plot = F)$cmdscale.out
mdsval <- as.tibble(z) %>% mutate(sname = row.names(z))
colnames(mdsval) <- c("leading logFC dim1", "leading logFC dim2", "sname")
mdsval <- left_join(mdsval, sdetails, by = "sname")
    mdsval %>% ggplot() + 
    geom_point(aes(x=`leading logFC dim1`, y = `leading logFC dim2`, shape = agegroup, colour = treatment), size = 3) + 
    geom_label_repel(aes(x=`leading logFC dim1`, y = `leading logFC dim2`,label = sex, fill=experiment), alpha = 0.5, size = 2, segment.color = 'black') +
      theme_classic()

mdsval <- as.tibble(z) %>% mutate(sname = row.names(z))
colnames(mdsval) <- c("leading logFC dim1", "leading logFC dim2", "sname")
mdsval <- left_join(mdsval, sdetails, by = "sname")
    mdsval %>% ggplot() + 
    geom_point(aes(x=`leading logFC dim1`, y = `leading logFC dim2`, shape = agegroup, colour = treatment), size = 3) + 
    geom_label_repel(aes(x=`leading logFC dim1`, y = `leading logFC dim2`,label = sname, fill=experiment), alpha = 0.5, size = 2, segment.color = 'black') +
      theme_classic()
