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

expr <- read.table(expr_file, header = TRUE, row.names = 1, check.names = FALSE)
datExpr <- t(expr)  # samples x genes
genes <- colnames(datExpr)
n_samples <- nrow(datExpr)

cor_mat <- cor(datExpr, use = "p")
adj <- adjacency(datExpr, type = "unsigned", power = power)
tom <- TOMsimilarityFromExpr(
    datExpr,
    networkType = "unsigned", TOMType = "signed", TOMDenom = "min", power = power
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

# --- soft threshold ---
power_vector <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
sft <- pickSoftThreshold(
    datExpr, networkType = "unsigned", powerVector = power_vector,
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
    power = power, networkType = "unsigned", TOMType = "signed", TOMDenom = "min",
    minModuleSize = min(20, ncol(datExpr) / 2), mergeCutHeight = 0.15,
    maxBlockSize = ncol(datExpr), deepSplit = 2, detectCutHeight = 0.995,
    numericLabels = FALSE, verbose = 0
)
writeLines(net$colors, file.path(outdir, "module_colors.tsv"))

cat("wrote reference outputs to", outdir, "\n")
