#!/usr/bin/env Rscript
# Regenerate the synthetic benchmark matrices used by `examples/bench.rs`.
# (These are intentionally not committed; the real demo data lives in Demo/.)

args <- commandArgs(trailingOnly = TRUE)
dir.create("benchdata", showWarnings = FALSE)

gen <- function(G, N, M, file, seed) {
    set.seed(seed)
    latent <- matrix(rnorm(M * N), M, N)
    e <- matrix(0, G, N)
    for (i in seq_len(G)) {
        m <- ((i - 1) %/% (G / M)) + 1
        e[i, ] <- latent[m, ] + rnorm(N, sd = 0.4)
    }
    e <- round(e, 6)
    rownames(e) <- paste0("g", seq_len(G))
    colnames(e) <- paste0("s", seq_len(N))
    write.table(data.frame(ID = rownames(e), e, check.names = FALSE),
        file = file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("wrote", file, "(", G, "x", N, ")\n")
}

gen(3000, 150, 30, "benchdata/expr3000.tsv", 1)
gen(8000, 200, 40, "benchdata/expr8000.tsv", 2)
