args = commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
	source_path = args[1]
	outfile_name = args[2]
} else {
	stop("Incorrect number of arguments was supplied (input file).n", call.=FALSE)
}
library(DESeq2)
library(GenomicFeatures)
counts = read.table(file.path(source_path,"counts","counts.txt"))
cond = read.table(file.path(source_path,"counts","coldata"),header = TRUE)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = cond, design = ~condition)
dds <- DESeq(dds)
dds$condition <-relevel(dds$condition, ref = "control")
res <- results(dds)
write.csv(res, file = file.path(source_path,outfile_name))
