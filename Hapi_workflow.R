

###########################################################
# ********************     Hapi     ********************* #
###########################################################


######### Introduction #########

# Hapi a novel easy-to-use and high-efficient algorithm that 
# only requires 3 to 5 gametes to reconstruct accurate and 
# high-resolution haplotypes of an individual. The gamete 
# genotype data may be generated from various platforms 
# including genotyping arrays and next generation sequencing 
# even with low-coverage. Hapi simply takes genotype data 
# of known hetSNPs in single gamete cells as input and report 
# the high-resolution haplotypes as well as the confidence level 
# of each phased hetSNPs. The package also includes a module 
# allowing downstream analyses and visualization of identified 
# crossovers in the gametes. 



#=========================================================#
#              1. Hapi package installation               #
#=========================================================#

setwd('~/Documents/Research/Git')
options(stringsAsFactors=FALSE)

### installation of Hapi locally
#Install 'HMM' package ahead
install.packages('HMM')
install.packages('Hapi_0.0.1.tar.gz', repos = NULL, type='source')
library(Hapi)

### installation of Hapi from Github
#install.packages('devtools')
#install.packages('HMM')

#devtools::install_github('Jialab-UCR/Hapi')






#=========================================================#
#                2. Haplotype phasing Module              #
#=========================================================#





###########################################################
##        2.1 Quick start (automatic phasing mode)       ##

### load example genotype data
data(gmt)
rownames(gmt) <- gmt$pos
head(gmt)


### automatic haplotype phasing
hapOutput <- hapiAutoPhase(gmt = gmt, code = 'atcg')
head(hapOutput)

###########################################################





###########################################################
##           2.2 Haplotype phasing step by step          ##


###### 2.2.1 Data preprocessing  ######

### covert A/T/C/G to 0/1
hetDa <- gmt[,1:4]
ref <- hetDa$ref
alt <- hetDa$alt

gmtDa <- gmt[,-(1:4)]
gmtDa <- base2num(gmt = gmtDa, ref = ref, alt = alt)
head(gmtDa)


### define HMM probabilities
library(HMM)
hmm = initHMM(States=c("S","D"), Symbols=c("s","d"), 
              transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
              emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
              startProbs = c(0.5,0.5))
hmm

### filter out genotyping errors
gmtDa <- hapiFilterError(gmt = gmtDa, hmm = hmm)


### select a subset of high-quality markers
gmtFrame <- hapiFrameSelection(gmt = gmtDa, n = 3) ###


### imputation
imputedFrame <- hapiImupte(gmt = gmtFrame, nSPT = 2, allowNA = 0)
head(imputedFrame)



###### 2.2.2 Draft haplotype inference  ######


### majority voting
draftHap <- hapiPhase(gmt = imputedFrame) ###
head(draftHap)

### check positions with cv-links
draftHap[draftHap$cvlink>=1,]

### identification of clusters of cv-links
cvCluster <- hapiCVCluster(draftHap = draftHap, cvlink = 2)
cvCluster

### determine hetSNPs in small regions involving multiple cv-links
filter <- c()
for (i in 1:nrow(cvCluster)) {
    filter <- c(filter, which (rownames(draftHap) >= cvCluster$left[i] & 
                                   rownames(draftHap) <= cvCluster$right[i]))
}

length(filter)

### filter out hetSNPs in complex regions and infer new draft haplotypes
if (length(filter) > 0) {
    imputedFrame <- imputedFrame[-filter, ]
    draftHap <- hapiPhase(imputedFrame)
} 


### proofreading by MPR
finalDraft <- hapiBlockMPR(draftHap = draftHap, gmtFrame = gmtFrame, cvlink = 1)
head(finalDraft)




###### 2.2.3 High-resolution haplotype assembly  ######

### high-resolution haplotype assembly
consensusHap <- hapiAssemble(draftHap = finalDraft, gmt = gmtDa)
head(consensusHap)


### assembly of regions at the ends of a chromosome
consensusHap <- hapiAssembleEnd(gmt = gmtDa, draftHap = finalDraft, 
                                consensusHap = consensusHap, k = 300)

### Haplotype 1
hap1 <- sum(consensusHap$hap1==0)
hap1
### Haplotype 2
hap2 <- sum(consensusHap$hap1==1)
hap2
### Number of unphased hetSNPs
hap7 <- sum(consensusHap$hap1==7)
hap7

### Accuracy
max(hap1, hap2)/sum(hap1, hap2)



### find hetSNP overlaps
snp <- which(rownames(hetDa) %in% rownames(consensusHap))

ref <- hetDa$ref[snp]
alt <- hetDa$alt[snp]

### convert back to A/T/C/G
consensusHap <- num2base(hap = consensusHap, ref = ref, alt = alt)
head(consensusHap)

### output all the information
hapOutput <- data.frame(gmt[snp,], consensusHap)
head(hapOutput)



###### 2.2.4 Visualization of haplotypes in single gamete cells  ######


### load haplotypes in each gamete cell
data(gamete11)
head(gamete11)

### load chromosome information of the genome
data(hg19)
head(hg19)


### view gamete cells 
hapiGameteView(chr = hg19, hap = gamete11)


###########################################################



#=========================================================#
#                3. Crossover analysis Module             #
#=========================================================#


###########################################################
##            3.1 Identification of crossovers           ##


### haplotypes
hap <- hapOutput[,10:11]
head(hap)

### gametes
gmt <- hapOutput[,5:9]
head(gmt)

### identify crossover
cvOutput <- hapiIdentifyCV(hap = hap, gmt = gmt)
cvOutput


###########################################################
##              3.2 Crossover visualization              ##



###### 3.2.1 Visualization of crossover resolution  ######


### load crossover table
data(crossover)
head(crossover)

hapiCVResolution(cv = crossover)


###### 3.2.2 Visualization of crossover distances  ######

hapiCVDistance(cv = crossover)


###### 3.2.3 Visualization of crossover map  ######

data(hg19)
head(hg19)

hapiCVMap(chr = hg19, cv = crossover, step = 5)
