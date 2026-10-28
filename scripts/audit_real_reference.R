# TEST-ONLY diagnostic extension of scripts/real_bench.R.
#!/usr/bin/env Rscript
# Real-data benchmark: FPKM -> log10(x+1) -> MAD top-N -> WGCNA -> module-trait.
# Matches the Rust `real_bench` example exactly (minModuleSize=30,
# mergeCutHeight=0.25, deepSplit=2, detectCutHeight=0.995).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: real_bench.R <fpkm.tsv> <trait.tsv> <outdir> [n_keep]")
fpkmpath <- args[1]
traitpath <- args[2]
outdir <- args[3]
n_keep <- if (length(args) >= 4) as.integer(args[4]) else 8000

suppressMessages(library(WGCNA))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
timings <- list()

t <- proc.time()
expr <- read.table(fpkmpath, header = TRUE, row.names = 1, check.names = FALSE)
stopifnot(!anyDuplicated(rownames(expr)), !anyDuplicated(colnames(expr)),all(is.finite(as.matrix(expr))))
datExpr <- t(expr)  # samples x genes
datExpr <- log10(datExpr + 1)
timings[["log10_normalize"]] <- (proc.time() - t)[["elapsed"]]

t <- proc.time()
mads <- apply(datExpr, 2, stats::mad)
keep <- order(mads, decreasing = TRUE)[seq_len(n_keep)]
datExpr <- datExpr[, keep, drop = FALSE]
genes <- colnames(datExpr)
timings[["mad_topN"]] <- (proc.time() - t)[["elapsed"]]

t <- proc.time()
cor_mat <- cor(datExpr, use = "p")
timings[["cor"]] <- (proc.time() - t)[["elapsed"]]

t <- proc.time()
tom <- TOMsimilarityFromExpr(datExpr, networkType = "unsigned", TOMType = "signed",
    TOMDenom = "min", power = 6)
timings[["tom"]] <- (proc.time() - t)[["elapsed"]]

write_mat <- function(m,file) write.table(data.frame(ID=rownames(m),m,check.names=FALSE),file,sep="\t",row.names=FALSE,quote=FALSE)
probe<-seq(1,ncol(datExpr),by=max(1,ncol(datExpr)%/%64));dimnames(tom)<-list(genes,genes)
write_mat(cor_mat[probe,,drop=FALSE],file.path(outdir,"cor_probe.tsv"));write_mat(tom[probe,,drop=FALSE],file.path(outdir,"tom_probe.tsv"))
tree<-hclust(as.dist(1-tom),method="average")
write.table(data.frame(ID=seq_along(tree$height),left=tree$merge[,1],right=tree$merge[,2],height=tree$height),file.path(outdir,"tree.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
dynamic<-cutreeDynamic(tree,distM=1-tom,deepSplit=2,minClusterSize=30,cutHeight=.995,pamRespectsDendro=TRUE,verbose=0)
write_mat(matrix(dynamic,ncol=1,dimnames=list(genes,"label")),file.path(outdir,"dynamic.tsv"))
sft<-pickSoftThreshold(datExpr,networkType="unsigned",powerVector=c(1:10,seq(12,20,2)),verbose=0)$fitIndices
write.table(sft,file.path(outdir,"sft.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
t <- proc.time()
source("validation/audit/instrument_reference.R"); reference_diagnostic <- instrument_blockwise()
net <- reference_diagnostic(datExpr, power = 6, networkType = "unsigned", TOMType = "signed",
    TOMDenom = "min", minModuleSize = 30, mergeCutHeight = 0.25,
    maxBlockSize = ncol(datExpr), deepSplit = 2, detectCutHeight = 0.995,
    numericLabels = FALSE, verbose = 0)
timings[["module_detection_total"]] <- (proc.time() - t)[["elapsed"]]

t <- proc.time()
trait_df <- read.table(traitpath, header = TRUE, row.names = 1, check.names = FALSE)
stopifnot(!anyDuplicated(rownames(trait_df)),setequal(rownames(datExpr),rownames(trait_df)))
datTraits <- trait_df[rownames(datExpr), , drop = FALSE]
me <- moduleEigengenes(datExpr, net$colors)
MEs <- as.matrix(me$eigengenes)
mt_cor <- cor(MEs, datTraits, use = "p")
mt_p <- corPvalueStudent(mt_cor, nrow(datExpr))
timings[["module_trait"]] <- (proc.time() - t)[["elapsed"]]

# Export actual blockwiseModules premerge and final states. unmergedColors is
# after cleanup/reassignment, NOT the raw dynamic tree cut.
write.table(data.frame(ID=genes,color=net$unmergedColors),file.path(outdir,"premerge.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
colnames(MEs)<-substring(colnames(MEs),3)
write_mat(MEs,file.path(outdir,"eigengenes.tsv"))
km<-as.matrix(signedKME(datExpr,MEs));colnames(km)<-colnames(MEs)
write_mat(km,file.path(outdir,"kme.tsv"))
adj<-abs(cor_mat)^6;diag(adj)<-0
within<-vapply(seq_along(genes),function(i)sum(adj[i,net$colors==net$colors[i]]),numeric(1))
write_mat(matrix(within,ncol=1,dimnames=list(genes,"within")),file.path(outdir,"intramodular.tsv"));rm(adj)
# outputs (match Rust format)
writeLines(genes, file.path(outdir, "filtered_genes.txt"))
write.table(data.frame(gene = genes, color = net$colors), file.path(outdir, "module_colors.txt"),
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(data.frame(ID = rownames(mt_cor), mt_cor, check.names = FALSE),
    file.path(outdir, "module_trait_cor.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(ID = rownames(mt_p), mt_p, check.names = FALSE),
    file.path(outdir, "module_trait_p.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(step = names(timings), seconds = unlist(timings)),
    file.path(outdir, "timings.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
sizes <- table(net$colors)
write.table(data.frame(module = names(sizes), size = as.integer(sizes)),
    file.path(outdir, "module_sizes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("R: modules =", paste(names(sizes), collapse = ","), "\n")
