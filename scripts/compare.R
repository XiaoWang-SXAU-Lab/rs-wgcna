#!/usr/bin/env Rscript
# Compare Rust outputs against R reference outputs with a tolerance.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("usage: compare.R <rust_dir> <r_dir> [tol]")
}
rust_dir <- args[1]
r_dir <- args[2]
tol <- if (length(args) >= 3) as.numeric(args[3]) else 1e-9

read_mat <- function(file) {
    as.matrix(read.table(file, header = TRUE, row.names = 1, check.names = FALSE))
}

files <- c(
    "cor.tsv", "adjacency.tsv", "tom.tsv", "sft.tsv",
    "eigengenes.tsv", "kme.tsv", "module_trait_cor.tsv", "module_trait_p.tsv"
)
ok <- TRUE
for (f in files) {
    rf <- file.path(rust_dir, f)
    rr <- file.path(r_dir, f)
    if (!file.exists(rf) || !file.exists(rr)) {
        cat(sprintf("%-24s SKIP (missing)\n", f))
        next
    }
    r <- read_mat(rf)
    ref <- read_mat(rr)
    if (!identical(dim(r), dim(ref))) {
        cat(sprintf("%-24s FAIL (dim %dx%d vs %dx%d)\n", f, nrow(r), ncol(r), nrow(ref), ncol(ref)))
        ok <- FALSE
        next
    }
    d <- abs(r - ref)
    maxd <- max(d, na.rm = TRUE)
    nbad <- sum(d > tol, na.rm = TRUE)
    pass <- maxd <= tol
    if (!pass) ok <- FALSE
    cat(sprintf(
        "%-24s %s  max_abs_diff=%.3e  (>tol: %d/%d entries)\n",
        f, if (pass) "PASS" else "FAIL", maxd, nbad, length(d)
    ))
}

# Module colors (discrete labels): compare exact match rate.
cf_r <- file.path(rust_dir, "module_colors.tsv")
cf_ref <- file.path(r_dir, "module_colors.tsv")
if (file.exists(cf_r) && file.exists(cf_ref)) {
    r <- readLines(cf_r)
    ref <- readLines(cf_ref)
    nmatch <- sum(r == ref)
    pass <- nmatch == length(r) && length(r) == length(ref)
    if (!pass) ok <- FALSE
    cat(sprintf(
        "%-24s %s  %d/%d genes identical\n",
        "module_colors.tsv", if (pass) "PASS" else "FAIL", nmatch, length(r)
    ))
}

if (!ok) {
    cat("RESULT: FAIL\n")
    quit(status = 1)
} else {
    cat("RESULT: PASS (all within tol =", format(tol, scientific = TRUE), ")\n")
}
