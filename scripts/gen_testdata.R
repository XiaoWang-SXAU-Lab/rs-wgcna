#!/usr/bin/env Rscript
# Generate synthetic test data: expression matrix (genes x samples), a module
# color vector, and a trait table (sample x traits). All written to testdata/.

args <- commandArgs(trailingOnly = TRUE)
outfile <- if (length(args) >= 1) args[1] else "testdata/expr.tsv"
outdir <- dirname(outfile)

n_genes <- 150
n_samples <- 40
n_modules <- 5
genes_per_module <- n_genes / n_modules

set.seed(42)
latent <- matrix(rnorm(n_modules * n_samples), nrow = n_modules, ncol = n_samples)
expr <- matrix(0, nrow = n_genes, ncol = n_samples)
for (i in seq_len(n_genes)) {
    m <- ((i - 1) %/% genes_per_module) + 1
    expr[i, ] <- latent[m, ] + rnorm(n_samples, sd = 0.3)
}
expr <- round(expr, 6)
rownames(expr) <- paste0("gene", seq_len(n_genes))
colnames(expr) <- paste0("sample", seq_len(n_samples))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# expression matrix
write.table(
    data.frame(ID = rownames(expr), expr, check.names = FALSE),
    file = outfile, sep = "\t", row.names = FALSE, quote = FALSE
)

# module colors (one per gene, in gene order)
colors <- rep(c("turquoise", "blue", "brown", "yellow", "green"), each = genes_per_module)
writeLines(colors, file.path(outdir, "colors.txt"))

# trait table: trait1 binary correlated with module 1; trait2 continuous noise
set.seed(7)
trait1 <- ifelse(latent[1, ] + rnorm(n_samples, sd = 0.5) > 0, 1, 0)
trait2 <- round(rnorm(n_samples), 6)
trait_df <- data.frame(sample = colnames(expr), trait1 = trait1, trait2 = trait2)
write.table(
    trait_df, file = file.path(outdir, "trait.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

cat(
    "wrote", outfile, "(", n_genes, "genes x", n_samples, "samples),",
    "colors.txt, trait.tsv\n"
)
