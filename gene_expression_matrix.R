#source("https://bioconductor.org/biocLite.R")
#biocLite("genefilter")
#biocLite("ballgown")
library(ballgown)
library(genefilter)


ball_path = list.files()
bg <- ballgown::ballgown(samples = ball_path)
genes = ballgown::gexpr(bg)
genes = genes[(grepl("ENS", rownames(genes))), ]
colnames(genes) = ball_path

bg_filt = subset(bg,"genefilter::rowVars(texpr(bg)) >= 1",genomesubset=TRUE)
trans = ballgown::texpr(bg_filt, "all")

write.csv(trans, "../Transcripts_Matrix.csv")
write.csv(genes, "../Genes_FPKM.csv")
