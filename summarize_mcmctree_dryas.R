#!/usr/bin/env Rscript

# Report each MCMCTree chain separately and both chains pooled by chromosome.
# Usage: Rscript summarize_mcmctree_dryas.R OUTPUT_DIR

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript summarize_mcmctree_dryas.R OUTPUT_DIR")
if (!requireNamespace("coda", quietly = TRUE)) stop("Install coda: install.packages('coda')")
root <- normalizePath(args[[1]], mustWork = TRUE)

read_chain <- function(path) {
  x <- utils::read.table(path, header = TRUE, check.names = FALSE,
                         comment.char = "", stringsAsFactors = FALSE)
  x <- x[, vapply(x, is.numeric, logical(1)), drop = FALSE]
  x <- x[, !tolower(names(x)) %in%
           c("gen", "iter", "iteration", "generation", "sample", "index"), drop = FALSE]
  if (!ncol(x)) stop("No numeric posterior columns in ", path)
  x
}

summarize_matrix <- function(x, chromosome, chain) {
  q <- apply(x, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  data.frame(
    chromosome = chromosome, chain = chain, parameter = names(x),
    mean = colMeans(x, na.rm = TRUE),
    sd = vapply(x, stats::sd, numeric(1), na.rm = TRUE),
    median = q["50%", ], lower_95 = q["2.5%", ], upper_95 = q["97.5%", ],
    stringsAsFactors = FALSE, row.names = NULL
  )
}

chain_rows <- pooled_rows <- parameter_diag_rows <- chromosome_diag_rows <- list()
for (chr in 1:9) {
  paths <- file.path(root, sprintf("chr%02d", chr), c("chain1", "chain2"), "mcmc.txt")
  if (any(!file.exists(paths))) stop("Missing chain output for chromosome ", chr)
  chains <- lapply(paths, read_chain)
  shared <- Reduce(intersect, lapply(chains, names))
  chains <- lapply(chains, function(x) x[, shared, drop = FALSE])
  ess <- lapply(chains, function(x) coda::effectiveSize(coda::mcmc(x)))

  for (chain in 1:2) {
    z <- summarize_matrix(chains[[chain]], chr, chain)
    z$ESS <- as.numeric(ess[[chain]][z$parameter])
    chain_rows[[length(chain_rows) + 1L]] <- z
  }

  ml <- coda::mcmc.list(coda::mcmc(chains[[1]]), coda::mcmc(chains[[2]]))
  rhat <- tryCatch(
    coda::gelman.diag(ml, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1],
    error = function(e) setNames(rep(NA_real_, length(shared)), shared)
  )
  z <- summarize_matrix(do.call(rbind, chains), chr, "pooled")
  z$Rhat <- as.numeric(rhat[z$parameter])
  pooled_rows[[length(pooled_rows) + 1L]] <- z

  pd <- data.frame(
    chromosome = chr, parameter = shared,
    chain1_ESS = as.numeric(ess[[1]][shared]),
    chain2_ESS = as.numeric(ess[[2]][shared]),
    Rhat = as.numeric(rhat[shared]), stringsAsFactors = FALSE
  )
  parameter_diag_rows[[length(parameter_diag_rows) + 1L]] <- pd
  chromosome_diag_rows[[length(chromosome_diag_rows) + 1L]] <- data.frame(
    chromosome = chr, parameters = length(shared),
    chain1_minimum_ESS = min(pd$chain1_ESS, na.rm = TRUE),
    chain2_minimum_ESS = min(pd$chain2_ESS, na.rm = TRUE),
    maximum_Rhat = max(pd$Rhat, na.rm = TRUE),
    both_chains_ESS_ge_200 = all(pd$chain1_ESS >= 200 & pd$chain2_ESS >= 200),
    all_Rhat_le_1.05 = all(is.finite(pd$Rhat)) && all(pd$Rhat <= 1.05), stringsAsFactors = FALSE
  )
}

chain_summary <- do.call(rbind, chain_rows)
pooled_summary <- do.call(rbind, pooled_rows)
parameter_diagnostics <- do.call(rbind, parameter_diag_rows)
chromosome_diagnostics <- do.call(rbind, chromosome_diag_rows)
write_tsv <- function(x, filename) utils::write.table(
  x, file.path(root, filename), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write_tsv(chain_summary, "mcmctree_chain_posterior_summary.tsv")
write_tsv(pooled_summary, "mcmctree_pooled_posterior_summary.tsv")
write_tsv(parameter_diagnostics, "mcmctree_parameter_convergence.tsv")
write_tsv(chromosome_diagnostics, "mcmctree_convergence_summary.tsv")
write_tsv(pooled_summary, "mcmctree_posterior_summary.tsv")
print(chromosome_diagnostics, row.names = FALSE)
cat("\nWrote chain-specific, pooled, and convergence summaries to", root, "\n")
