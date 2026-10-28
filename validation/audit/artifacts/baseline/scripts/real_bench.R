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

t <- proc.time()
net <- blockwiseModules(datExpr, power = 6, networkType = "unsigned", TOMType = "signed",
    TOMDenom = "min", minModuleSize = 30, mergeCutHeight = 0.25,
    maxBlockSize = ncol(datExpr), deepSplit = 2, detectCutHeight = 0.995,
    numericLabels = FALSE, verbose = 0)
timings[["module_detection_total"]] <- (proc.time() - t)[["elapsed"]]

t <- proc.time()
trait_df <- read.table(traitpath, header = TRUE, row.names = 1, check.names = FALSE)
datTraits <- trait_df[rownames(datExpr), , drop = FALSE]
me <- moduleEigengenes(datExpr, net$colors)
MEs <- as.matrix(me$eigengenes)
mt_cor <- cor(MEs, datTraits, use = "p")
mt_p <- corPvalueStudent(mt_cor, nrow(datExpr))
timings[["module_trait"]] <- (proc.time() - t)[["elapsed"]]

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
