
#######################################################
################# test GDCRNATools ####################
sessionInfo()
ls()
rm(list=ls())

######## set working directory #######
setwd('~/Documents/Research/Cancer_Genomics/Data/Bladder/ceTCGA/newPackage3/')


#remove.packages('ceTCGAtools')

####### installation #######

install.packages('GDCRNATools_0.99.12.tar.gz', repos = NULL, type='source')
library(GDCRNATools)


### if the last step doesn't work, try this
#install.packages(devtools)
#library(devtools)
#install('GDCRNATools')

sessionInfo()



####### load required libraries (optional) #######
#library(survival)
#library(survminer)
#library(pathview)
#library(rjson)
#library(edgeR)
#library(limma)
#library(ggplot2)
#library(gplots)
#library(shiny)
#library(clusterProfiler)
#library(org.Hs.eg.db)
#library(biomaRt)


# Quick start

library(DT)
### load RNA counts data
data(rnaCounts)

### load miRNAs counts data
data(mirCounts)


####### RNAseq data #######
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)

####### miRNAs data #######
mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)


metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
datatable(as.data.frame(metaMatrix.RNA[1:5,]), extensions = 'Scroller',
          options = list(scrollX = TRUE, deferRender = TRUE, scroller = TRUE))



####### Parse miRNAs metadata #######
metaMatrix.MIR <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'miRNAs', 
                                   write.meta = FALSE)

metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)
metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)




## Differential gene expression analysis

DEGAll <- gdcDEAnalysis(counts     = rnaCounts, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma')
datatable(as.data.frame(DEGAll), 
          options = list(scrollX = TRUE, pageLength = 5))

#### DE long-noncoding
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding')

#### DE protein coding genes
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')



## ceRNAs network analysis

ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = 'starBase', 
                          pc.targets  = 'starBase', 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)

datatable(as.data.frame(ceOutput), 
          options = list(scrollX = TRUE, pageLength = 5))


## Report ceRNA network to Cytoscape

ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 
                      & ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]

edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
datatable(as.data.frame(edges), 
          options = list(scrollX = TRUE, pageLength = 5))



nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
datatable(as.data.frame(nodes), 
          options = list(scrollX = TRUE, pageLength = 5))


################################################################################
################################################################################

# Case study: TCGA-CHOL {#case}

## Data download

### Automatic download  
  
project <- 'TCGA-CHOL'
rnadir <- paste(project, 'RNAseq', sep='/')
mirdir <- paste(project, 'miRNAs', sep='/')
mirdir
####### Download RNAseq data #######
gdcRNADownload(project.id     = 'TCGA-CHOL', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               directory      = rnadir)

####### Download miRNAs data #######
gdcRNADownload(project.id     = 'TCGA-CHOL', 
               data.type      = 'miRNAs', 
               write.manifest = FALSE,
               directory      = mirdir)


## Data organization and DE analysis

### Parse metadata

####### Parse RNAseq metadata #######
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)


####### Filter duplicated samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)



####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)





####### Parse miRNAs metadata #######
metaMatrix.MIR <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'miRNAs', 
                                   write.meta = FALSE)



####### Filter duplicated samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)




####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)



### Merge raw counts data

####### Merge RNAseq data #######
rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                         path      = rnadir, 
                         data.type = 'RNAseq')

####### Merge miRNAs data #######
mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                         path      = mirdir,
                         data.type = 'miRNAs')


### TMM normalization and voom transformation


####### RNAseq data #######
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)

####### miRNAs data #######
mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)



### Differential gene expression analysis


DEGAll <- gdcDEAnalysis(counts     = rnaCounts, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma')


#data(DEGAll)


### All DEGs
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all')


#### DE long-noncoding
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding')



#### DE protein coding genes
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')



## Competing endogenous RNAs network analysis


### Hypergeometric test

ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = 'starBase', 
                          pc.targets  = 'starBase', 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)


#### ceRNAs network analysis using user-provided datasets

#### load miRNA-lncRNA interactions
data(lncTarget)

#### load miRNA-mRNA interactions
data(pcTarget)
pcTarget[1:3]



ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = lncTarget, 
                          pc.targets  = pcTarget, 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)





### Network visulization in Cytoscape


ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 & 
                          ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]

edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')

write.table(edges, file='edges.txt', sep='\t', quote=F)
write.table(nodes, file='nodes.txt', sep='\t', quote=F)




### Correlation plot on a local webpage

shinyCorPlot(gene1    = rownames(deLNC), 
             gene2    = rownames(dePC), 
             rna.expr = rnaExpr, 
             metadata = metaMatrix.RNA)

## Other downstream analyses


### Univariate survival analysis



#### CoxPH analysis

####### CoxPH analysis #######
survOutput <- gdcSurvivalAnalysis(gene     = rownames(deALL), 
                                  method   = 'coxph', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA)



#### KM analysis

####### KM analysis #######
survOutput <- gdcSurvivalAnalysis(gene     = rownames(deALL), 
                                  method   = 'KM', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA, 
                                  sep      = 'median')



#### KM plot on a local webpage by shinyKMPlot

shinyKMPlot(gene = rownames(deALL), rna.expr = rnaExpr, 
            metadata = metaMatrix.RNA)



### Functional enrichment analysis

enrichOutput <- gdcEnrichAnalysis(gene = rownames(deALL), simplify = TRUE)


#### Barplot

#data(enrichOutput)



gdcEnrichPlot(enrichOutput, type = 'bar', category = 'GO', num.terms = 10)



#### Bubble plot

gdcEnrichPlot(enrichOutput, type='bubble', category='GO', num.terms = 10)



#### View pathway maps on a local webpage

library(pathview)
deg <- deALL$logFC
names(deg) <- rownames(deALL)
pathways <- as.character(enrichOutput$Terms[enrichOutput$Category=='KEGG'])



shinyPathview(deg, pathways = pathways, directory = 'pathview')

