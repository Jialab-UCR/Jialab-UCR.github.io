
###########################################################
# *****************     GDCRNATools     ***************** #
###########################################################

# GDCRNATools is an R package which provides a standard, 
# easy-to-use and comprehensive pipeline for downloading, 
# organizing, and integrative analyzing RNA expression data 
# in the GDC portal with an emphasis on deciphering the 
# lncRNA-mRNA related ceRNAs regulatory network in cancer.


# Here we provide code of the basic steps for data analysis
# by GDCRNATools. Detailed instructions can be found here:
# http://htmlpreview.github.io/?https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/GDCRNATools_manual.html




#=========================================================#
#          1. GDCRNATools package installation            #
#=========================================================#

# Get the current working directory, make sure that it is 
# writable, otherwise, change to a new directory
getwd()
#setwd(workingDirectory)

# installation of GDCRNATools from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("GDCRNATools")

library(GDCRNATools)






#=========================================================#
#                      2. Quick start                     #
#=========================================================#

# A small internal dataset is used here to show the most basic 
# steps for ceRNAs network analysis in GDCRNATools




###########################################################
##         2.1 Normalization of HTSeq-Counts data        ##


### load RNA counts data
data(rnaCounts)
rnaCounts[1:5,1:5]
    
### load miRNAs counts data
data(mirCounts)
mirCounts[1:5,1:5]
    
### Normalization of RNAseq data
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)
rnaExpr[1:5,1:5]
    
### Normalization of miRNAs data
mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)
mirExpr[1:5,1:5]
    

###########################################################





###########################################################
##          2.2 Parse and filter RNAseq metadata         ##

metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
metaMatrix.RNA[1:5,]


###########################################################





###########################################################
##             2.3 ceRNAs network analysis               ##


### Identification of differentially expressed genes ###

DEGAll <- gdcDEAnalysis(counts     = rnaCounts, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma')
DEGAll[1:5,]
    
# All DEGs
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all')
deALL[1:5,]
    
# DE long-noncoding genes
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding')
deLNC[1:5,]

# DE protein coding genes
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')
dePC[1:5,]

    
    
########### ceRNAs network analysis of DEGs ############
ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = 'starBase', 
                          pc.targets  = 'starBase', 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)

ceOutput[1:5,]

    
### Export ceRNAs network to Cytoscape
ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 
    & ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]

# Export edges
edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
edges[1:5,]

# Export nodes
nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
nodes[1:5,]


###########################################################






#=========================================================#
#               3. Case study: TCGA-CHOL                  #
#=========================================================#


###########################################################
##                   3.1 Download data                   ##


# set up directories for downloaded data
project <- 'TCGA-CHOL'
rnadir <- paste(project, 'RNAseq', sep='/')
mirdir <- paste(project, 'miRNAs', sep='/')
    
### Download RNAseq data
gdcRNADownload(project.id     = 'TCGA-CHOL', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               directory      = rnadir)

### Download miRNAs data
gdcRNADownload(project.id     = 'TCGA-CHOL', 
               data.type      = 'miRNAs', 
               write.manifest = FALSE,
               directory      = mirdir)

###########################################################






###########################################################
##                 3.2 Data organization                 ##


### Parse RNAseq metadata
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

# Filter duplicated samples in RNAseq metadata
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
# Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)

    

### Parse miRNAs metadata
metaMatrix.MIR <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'miRNAs', 
                                   write.meta = FALSE)

# Filter duplicated samples in miRNAs metadata
metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)
# Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata
metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)

    


### Merge raw counts data
# Merge RNAseq data
rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                         path      = rnadir, 
                         data.type = 'RNAseq')

# Merge miRNAs data
mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                         path      = mirdir,
                         data.type = 'miRNAs')


### TMM normalization and voom transformation
# Normalization of RNAseq data
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)

# Normalization of miRNAs data
mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)


### Differential gene expression analysis
DEGAll <- gdcDEAnalysis(counts     = rnaCounts, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma')
#data(DEGAll)

# All DEGs
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all')

# DE long-noncoding
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding')

# DE protein coding genes
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')


###########################################################






###########################################################
##    3.3 Competing endogenous RNAs network analysis     ##


### The 3 steps of ceRNAs network analysis:
# Hypergeometric test
# Pearson correlation analysis
# Regulation pattern analysis

### All of the 3 steps can be performed in a single function


### ceRNAs network analysis using internal databases
ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = 'starBase', 
                          pc.targets  = 'starBase', 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)



### ceRNAs network analysis using user-provided datasets
# load miRNA-lncRNA interactions
data(lncTarget)
lncTarget[1:3]

# load miRNA-mRNA interactions
data(pcTarget)
pcTarget[1:3]

ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = lncTarget, 
                          pc.targets  = pcTarget, 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)
    

### Network visulization in Cytoscape

# Filter potential ceRNA interactions
ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 & 
    ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]


# Edges and nodes can be simply imported into Cytoscape 
# for network visualization
edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
edges[1:5,]

nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
nodes[1:5,]

write.table(edges, file='edges.txt', sep='\t', quote=F) ### Network of Cytoscape
write.table(nodes, file='nodes.txt', sep='\t', quote=F) ### Table of Cytoscape


### Correlation plot on a local webpage

shinyCorPlot(gene1    = rownames(deLNC), 
             gene2    = rownames(dePC), 
             rna.expr = rnaExpr, 
             metadata = metaMatrix.RNA)


###########################################################
    




###########################################################
##             3.4 Other downstream analyses             ##


############### Univariate survival analysis ##############

# CoxPH analysis
survOutput <- gdcSurvivalAnalysis(gene     = rownames(deALL), 
                                  method   = 'coxph', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA)

    
# KM analysis
survOutput <- gdcSurvivalAnalysis(gene     = rownames(deALL), 
                                  method   = 'KM', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA, 
                                  sep      = 'median')

# KM plot on a local webpage by shinyKMPlot
shinyKMPlot(gene = rownames(deALL), rna.expr = rnaExpr, 
            metadata = metaMatrix.RNA)

 


############## Functional enrichment analysis #############

### All the functional enrichment analyses can be 
### performed in a single function, including:
# Gene Ontology (BP, CC, MF) analysis
# KEGG pathway analysis
# Disease Ontology analysis

enrichOutput <- gdcEnrichAnalysis(gene = rownames(deALL), simplify = TRUE)

#data(enrichOutput)

# Barplot
gdcEnrichPlot(enrichOutput, type = 'bar', category = 'GO', num.terms = 10)

# Bubble plot
gdcEnrichPlot(enrichOutput, type='bubble', category='GO', num.terms = 10)
    

# View pathway maps on a local webpage
library(pathview)

deg <- deALL$logFC
names(deg) <- rownames(deALL)
pathways <- as.character(enrichOutput$Terms[enrichOutput$Category=='KEGG'])

shinyPathview(deg, pathways = pathways, directory = 'pathview')


###########################################################


