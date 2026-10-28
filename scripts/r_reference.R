#!/usr/bin/env Rscript
# R reference computation matching the Rust harness defaults:
#   networkType = unsigned, TOMType = signed, TOMDenom = min, power = 6.
# Outputs cor / adjacency / TOM / SFT / eigengenes / kME / module-trait.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("usage: r_reference.R <expr.tsv> <outdir> <power> <colors.txt> <trait.tsv>")
}
expr_file <- args[1]
outdir <- args[2]
power <- as.numeric(args[3])
colors_file <- if (length(args) >= 4) args[4] else NA
trait_file <- if (length(args) >= 5) args[5] else NA

suppressMessages(library(WGCNA))
nt <- if(length(args)>=6) args[6] else "unsigned"
tt <- if(length(args)>=7) args[7] else "signed"
td <- if(length(args)>=8) args[8] else "min"
deep <- if(length(args)>=9) as.integer(args[9]) else 2
merge_height <- if(length(args)>=10) as.numeric(args[10]) else .15

expr <- read.table(expr_file, header = TRUE, row.names = 1, check.names = FALSE)
datExpr <- t(expr)  # samples x genes
genes <- colnames(datExpr)
n_samples <- nrow(datExpr)

stopifnot(!anyDuplicated(genes), !anyDuplicated(rownames(datExpr)), all(is.finite(as.matrix(datExpr))))
cor_mat <- cor(datExpr, use = "p")
adj <- adjacency(datExpr, type = nt, power = power)
tom <- TOMsimilarityFromExpr(
    datExpr,
    networkType = nt, TOMType = tt, TOMDenom = td, power = power
)
dimnames(cor_mat) <- list(genes, genes)
dimnames(adj) <- list(genes, genes)
dimnames(tom) <- list(genes, genes)

write_square <- function(m, file) {
    write.table(
        data.frame(ID = rownames(m), m, check.names = FALSE),
        file = file, sep = "\t", row.names = FALSE, quote = FALSE
    )
}
write_mat <- function(m, file) {
    write.table(
        data.frame(ID = rownames(m), m, check.names = FALSE),
        file = file, sep = "\t", row.names = FALSE, quote = FALSE
    )
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write_square(cor_mat, file.path(outdir, "cor.tsv"))
write_square(adj, file.path(outdir, "adjacency.tsv"))
write_square(tom, file.path(outdir, "tom.tsv"))

# Extra labelled test diagnostics, without changing reference algorithms.
tree <- hclust(as.dist(1-tom),method="average")
write.table(data.frame(ID=seq_along(tree$height),left=tree$merge[,1],right=tree$merge[,2],height=tree$height),file.path(outdir,"tree.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
dyn <- cutreeDynamic(tree,distM=1-tom,deepSplit=deep,minClusterSize=20,cutHeight=.995,pamRespectsDendro=TRUE,verbose=0)
write_mat(matrix(dyn,ncol=1,dimnames=list(genes,"label")),file.path(outdir,"dynamic.tsv"))
write_mat(matrix(rowSums(adj)-diag(adj),ncol=1,dimnames=list(genes,"total")),file.path(outdir,"connectivity.tsv"))

# --- soft threshold ---
power_vector <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
sft <- pickSoftThreshold(
    datExpr, networkType = nt, powerVector = power_vector,
    RsquaredCut = 0.85, verbose = 0
)
write.table(
    sft$fitIndices, file = file.path(outdir, "sft.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

# --- eigengenes + kME (require colors) ---
if (!is.na(colors_file)) {
    colors <- scan(colors_file, what = "", quiet = TRUE)
    stopifnot(length(colors) == length(genes))
    me <- moduleEigengenes(datExpr, colors)
    MEs <- as.matrix(me$eigengenes)
    colnames(MEs) <- substring(colnames(MEs), 3)  # strip "ME" prefix
    write_mat(MEs, file.path(outdir, "eigengenes.tsv"))

    kME <- as.matrix(signedKME(datExpr, me$eigengenes))
    colnames(kME) <- substring(colnames(kME), 4)  # strip "kME" prefix
    write_mat(kME, file.path(outdir, "kme.tsv"))

    # --- module-trait (require trait) ---
    if (!is.na(trait_file)) {
        trait_df <- read.table(trait_file, header = TRUE, row.names = 1, check.names = FALSE)
        stopifnot(!anyDuplicated(rownames(trait_df)), setequal(rownames(trait_df),rownames(datExpr)))
        datTraits <- trait_df[rownames(datExpr), , drop = FALSE]
        mt_cor <- cor(MEs, datTraits, use = "p")
        mt_p <- corPvalueStudent(mt_cor, n_samples)
        write_mat(mt_cor, file.path(outdir, "module_trait_cor.tsv"))
        write_mat(mt_p, file.path(outdir, "module_trait_p.tsv"))
    }
}

# --- full module detection (blockwiseModules, single block) ---
net <- blockwiseModules(
    datExpr,
    power = power, networkType = nt, TOMType = tt, TOMDenom = td,
    minModuleSize = min(20, ncol(datExpr) / 2), mergeCutHeight = merge_height,
    maxBlockSize = ncol(datExpr), deepSplit = deep, detectCutHeight = 0.995,
    numericLabels = FALSE, verbose = 0
)
writeLines(net$colors, file.path(outdir, "module_colors.tsv"))

cat("wrote reference outputs to", outdir, "\n")
