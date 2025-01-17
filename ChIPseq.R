library(RColorBrewer)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
files<-dir(pattern = "bed")
names(files) <- c("H3K27Ac", "H3K27me3", "H3K4me1")
#read the first bed file
peaks <- readPeakFile(files[1])
peaks

#library(GenoGAM)
#q=GRanges(seqnames = c("chr1", "chr2"), ranges = NULL, strand = "*")
#peaks2 <-peaks[seqnames(peaks) %in% c("chr1", "chr2")]

## covplot(peaks, weightCol = "V5") # all chromosomes
covplot(peaks, weightCol = "V5", chrs = paste0("chr", c(1:19,"X","Y")), 
        title = "H3K27Ac Peaks over Chromosomes")

promoter <- getPromoters(TxDb = TxDb, upstream = 3000, downstream = 3000)
TagMatrix <- getTagMatrix(peaks, windows = promoter)


# plot heatmap
tagHeatmap(TagMatrix)
peakHeatmap(files[[1]], TxDb=TxDb, upstream=3000, downstream=3000)

# plot profile
plotAvgProf(TagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
# one step from bed file to average profile plot
plotAvgProf2(files[[2]], TxDb=TxDb, upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
# ChIP binding profiles with confidence interval
#plotAvgProf(TagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)

# Peak Annotation
peakAnno <- annotatePeak(files[[1]], tssRegion=c(-3000, 3000),
                         TxDb=TxDb, annoDb="org.Mm.eg.db")

#  Genomic Annotation by pieplot & barplot
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno) # vennpie
upsetplot(peakAnno) # upsetplot
upsetplot(peakAnno, vennpie=TRUE) # upsetplot with vennpie

#Visualize distribution of TF-binding loci relative to TSS
plotDistToTSS(peakAnno,
              title="Distribution of H3K27Ac relative to TSS")

### Functional enrichment analysis
library(ReactomePA)
pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId, organism = "mouse")
head(pathway1, 2)
gene <- seq2gene(peaks, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=TxDb)
pathway2 <- enrichPathway(gene, organism = "mouse")
head(pathway2, 2)
dotplot(pathway2)
emapplot(pathway2)

### ChIP peak data set comparison
# Profile of several ChIP peak data binding to TSS region
promoter <- getPromoters(TxDb=TxDb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), facet = "column")

tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=brewer.pal(4, "Set1"))

peakAnnoList <- lapply(files, annotatePeak, TxDb=TxDb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)

### Functional profiles comparison
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")

### Overlap of peaks and annotated genes
genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)

### Statistical testing of ChIP seq overlap
## Peak overlap enrichment analysis
enrichPeakOverlap(queryPeak     = files[[4]],
                  targetPeak    = unlist(files[1:3]),
                  TxDb          = TxDb,
                  pAdjustMethod = "BH",
                  nShuffle      = 1000,
                  chainFile     = NULL,
                  verbose       = FALSE)

devtools::session_info()
