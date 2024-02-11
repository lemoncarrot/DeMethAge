library(FlowSorted.BloodExtended.EPIC)
library(EpiDISH)
library(data.table)
#install.packages("/Users/kchen/OneDrive/Documents/deconvolution_old/FlowSorted.BloodExtended.EPIC_1.1.2.tar.gz", repos = NULL, type = "source")
data("FlowSorted.BloodExtended.450klegacy.compTable")
reference <- FlowSorted.BloodExtended.450klegacy.compTable
filter <- data("IDOLOptimizedCpGsBloodExtended450k")

x <- read.table(gzfile("GSE62992_100ids.485578probes.raw.Signal_B.NA.txt.gz"), sep="\t")
b <- read.table(gzfile("GSE62992_100ids.485578probes.raw.Signal_A.NA.txt.gz"), sep="\t")
#process into beta values

x <- readRDS("HANNUM_DATA.rds")
x <- readRDS("GSE174555_proc.rds")
colnam <- colnames(x)
x <- as.data.table(t(x))
rownames(x) <- colnam

filter <- intersect(rownames(x), IDOLOptimizedCpGsBloodExtended450k)

to_deconvolve <- x[filter, ]

#methylation data <- row as cpg site, column as sample
#reference <- row as cpg site, column as sample, colname as cell type

out.l <- epidish(beta.m = hannum, ref.m = reference, method = "RPC") 
