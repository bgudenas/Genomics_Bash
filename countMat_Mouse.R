#!/usr/bin/env Rscript
library(stringr)
library(DESeq2)

readspergene = list.files()[grepl("ReadsPerGene", list.files())]
paste(length(readspergene), " SAMPLES")

geneMat = c()
IDs =c()
for ( i in readspergene){
  vals = read.table(i, header = FALSE, skip = 4)
  print(i)
  print(dim(vals))
  print(sum(vals$V4) > sum(vals$V3))
  genes = vals$V1
  ### Select 4th column for Illumina stranded
  vals = vals[ ,c(4)]
  geneMat = cbind(geneMat, vals)
  
  IDs = c(IDs,  unlist(strsplit(i, "_Reads" ))[1] )
}
rownames(geneMat) = genes 
colnames(geneMat) = IDs


meta = as.data.frame(colnames(geneMat))

DE =DESeqDataSetFromMatrix(countData = geneMat, colData = meta, design = ~ 1) 
denames = rownames(mcols(DE, use.names = TRUE))

vecraw = read.csv("/home/bgudenas/Annots/Mouse/Gene_lengths_mm38v93.csv")
vec = as.vector(vecraw$x)
names(vec) = vecraw$X

vec0 = vec[order(match(names(vec), denames))]
all(denames == names(vec0))
mcols(DE)$basepairs = vec0
FPKMs = fpkm(DE)


dir =basename(getwd())

setwd("../../")
write.csv(geneMat, paste0("./Results/", dir, "_CountMat.csv"))
write.csv(FPKMs, paste0("./Results/", dir, "_FPKMMat.csv"))
