
#PC3: https://www.encodeproject.org/experiments/ENCSR359LOD/

#PC3 (PC-3) is a human prostate cancer cell line used in prostate cancer research and drug development. PC3 cells are useful in investigating biochemical changes in advanced prostate cancer cells and in assessing their response to chemotherapeutic agents.

#PC9: https://www.encodeproject.org/experiments/ENCSR243INX/

#PC-9 is lung adenocarcinoma cell line with a deletion in exon 19 of the EGFR gene that exhibits high sensitivity to TKIs.

Goals:
  
  
  1. Make a venn diagram comparing the overlap of binding sites for your two ChIP-Seq datasets
2. Make a metaplot of your two datasets around the TSS.
3. Annotate the peaks for genomic features such as intron, exon, 3’UTR, etc and compare the annotations between your two datasets.
4. Assign peaks to genes – then perform pathway enrichment.
6. What genes and pathways/genesets shared between these datasets? What pathways differ? What is your interpretation of these results? What future directions could you propose to follow up on these findings? (No right answers to these questions, just important to think through this part).

overlaps.anno <- annotatePeakInBatch(overlaps, 
                                     AnnotationData=annoData, 
                                     output="nearestBiDirectionalPromoters",
                                     bindingRegion=c(-2000, 500))
library(org.Hs.eg.db)
library(DBI)

overlaps.anno <- addGeneIDs(overlaps.anno,
                            "org.Hs.eg.db",
                            IDs2Add = "entrez_id")

over <- getEnrichedGO(overlaps.anno, orgAnn="org.Hs.eg.db", condense=TRUE)
enrichmentPlot(over,n= 10)

peaks <- GRangesList(rep1=FOXA1,
                     rep2=ER)
genomicElementDistribution(peaks, 
                           TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000))


library(ChIPseeker)
library(clusterProfiler)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#Second
pc9<- toGRanges("/Users/cordell/Downloads/ENCFF267BQD.bed.gz" , format="narrowPeak", header=FALSE)
pc3<- toGRanges("/Users/cordell/Downloads/ENCFF450JNO.bed.gz" , format="narrowPeak", header=FALSE)

PC9 <- unique(pc9)
PC3 <- unique(pc3)

ol <- findOverlapsOfPeaks(PC3,PC9)
ol <- addMetadata(ol, colNames="score", FUN=mean) 

#Thrid
makeVennDiagram(ol, fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2")) # label color, keep same as circle border color

#fourth
library(EnsDb.Hsapiens.v75) ##(hg19)
## create annotation file from EnsDb or TxDb
annoData <- toGRanges(EnsDb.Hsapiens.v75, feature="gene")

#Fifth
pie1(table(ol$overlappingPeaks[["PC3///PC9"]]$overlapFeature))

#Sixth
overlaps <- ol$peaklist[["PC3///PC9"]]

binOverFeature(overlaps, annotationData=annoData,
               radius=5000, nbins=200, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS")

#Seventh
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

peaks <- GRangesList(rep1=PC9,
                     rep2=PC3)

genomicElementDistribution(peaks, 
                           TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000))

#Eighth
out <- genomicElementDistribution(overlaps, 
                                  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                  promoterRegion=c(upstream=2000, downstream=500),
                                  geneDownstream=c(upstream=0, downstream=2000),
                                  promoterLevel=list(
                                    # from 5' -> 3', fixed precedence 3' -> 5'
                                    breaks = c(-2000, -1000, -500, 0, 500),
                                    labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
                                               "upstream <500b", "TSS - 500b"),
                                    colors = c("#FFE5CC", "#FFCA99", 
                                               "#FFAD65", "#FF8E32")))










