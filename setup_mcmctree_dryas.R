#!/usr/bin/env Rscript

# Prepare nine simple MCMCTree analyses for Dryas.
# Design: optional chromosome-specific BPP crown-Dryas calibration, independent
# rates (clock=2), and two replicate chains per chromosome. Set the environment
# fourth argument selects the calibration scenario; see README.md.
#
# Usage:
#   Rscript setup_mcmctree_dryas.R ALIGNMENT_DIR TREE_DIR OUTPUT_DIR SCENARIO
#
# Input files must contain chr01 ... chr09 in their basenames. TREE_DIR should
# contain one constrained IQ-TREE .treefile per chromosome. ALIGNMENT_DIR should
# contain one PHYLIP alignment (.phy, .phylip, .txt, or .aln) per chromosome.

suppressPackageStartupMessages({
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required: install.packages('ape')")
  }
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript setup_mcmctree_dryas.R ALIGNMENT_DIR TREE_DIR OUTPUT_DIR SCENARIO")
}

dataset <- match.arg(Sys.getenv("DATASET", "15tip"), c("15tip", "p28"))
alignment_dir <- normalizePath(args[[1]], mustWork = TRUE)
tree_dir <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- args[[3]]
if (dir.exists(output_dir) && length(list.files(output_dir, all.files=TRUE, no..=TRUE))) stop("Use a new empty output directory")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

scenario <- match.arg(args[[4]], c("fossil", "combined", "baseline", "japan_stem", "outgroup_only"))
use_bpp_crown_dryas <- scenario %in% c("combined", "baseline")
use_fossils <- scenario %in% c("fossil", "combined", "japan_stem")
prior_only <- Sys.getenv("PRIOR_ONLY", "0") == "1"

dryas_tips <- c(
  "drum_BABY63", "sp_POB10", "cauc_ODP7438", "oxyd_ODP43154",
  "octo_END01", "cren_ODP53249", "ajan_ODP48916", "alas_ODP43138",
  "hook_GLA03", "intg_ODP7653", "octo_JAT08"
)
outgroup_tips <- c("Cermon_DW847", "Chafol_1611", "Pursta_DW848", "Purtri_DW846b")
core_tips <- setdiff(dryas_tips, "drum_BABY63")
# AIH includes the ajanensis/alaskensis and integrifolia/hookeriana branches.
aih_tips <- c("ajan_ODP48916", "alas_ODP43138", "hook_GLA03", "intg_ODP7653")
if (dataset == "p28") {
  reps <- read.delim(Sys.getenv("REPRESENTATIVES"), stringsAsFactors=FALSE)
  if (!all(c("sampleID", "population") %in% names(reps)) || nrow(reps)!=28 || anyDuplicated(reps$sampleID)) stop("Invalid representatives table")
  # Explicit membership is mandatory: do not infer clades from population names.
  membership <- read.delim(Sys.getenv("CLADE_MEMBERSHIP"), stringsAsFactors=FALSE)
  if (!all(c("sampleID", "calibration_clade") %in% names(membership)) ||
      anyDuplicated(membership$sampleID) || !setequal(membership$sampleID,reps$sampleID) ||
      any(!membership$calibration_clade %in% c("drummondii", "AIH", "other_core")))
    stop("clade_membership.tsv must use calibration_clade to assign every p28 sample to drummondii, AIH or other_core")
  dryas_tips <- reps$sampleID
  core_tips <- membership$sampleID[membership$calibration_clade != "drummondii"]
  aih_tips <- membership$sampleID[membership$calibration_clade == "AIH"]
  if (length(core_tips)>=28 || length(aih_tips)<2) stop("Missing drummondii or AIH membership")
}
expected_tips <- c(dryas_tips, outgroup_tips)

# BPP H4d posterior means and reported intervals, converted with
# u=3e-9/site/generation and g=10 yr.
bpp_tau_mean <- c(0.00362, 0.00358, 0.00350, 0.00357, 0.00342,
                  0.00363, 0.00387, 0.00383, 0.00370)
bpp_tau_lower <- c(0.00342, 0.00336, 0.00327, 0.00332, 0.00319,
                   0.00336, 0.00360, 0.00354, 0.00345)
bpp_tau_upper <- c(0.00382, 0.00379, 0.00372, 0.00381, 0.00366,
                   0.00389, 0.00414, 0.00411, 0.00395)
tau_to_ma <- 10 / 3e-9 / 1e6
bpp_dryas_mean_ma <- bpp_tau_mean * tau_to_ma
bpp_dryas_lower_ma <- bpp_tau_lower * tau_to_ma
bpp_dryas_upper_ma <- bpp_tau_upper * tau_to_ma

# MCMCTree calibration numbers below use 100 Ma as one time unit.
# Crown Dryadoideae spans the Cercocarpus fossil minimum to the upper 95% age
# estimate reported by Zhang et al. (2017).
root_min_ma <- 34.07
root_max_ma <- 48.70
cerc_min_ma <- 34.07

extract_chr <- function(path) {
  hit <- regmatches(basename(path), regexec("chr0*([1-9])", basename(path), ignore.case = TRUE))[[1]]
  if (length(hit) == 2L) as.integer(hit[[2]]) else NA_integer_
}

find_one_per_chr <- function(directory, pattern, kind) {
  files <- list.files(directory, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  chr <- vapply(files, extract_chr, integer(1))
  files <- files[!is.na(chr)]
  chr <- chr[!is.na(chr)]
  answer <- character(9)
  for (i in 1:9) {
    candidates <- files[chr == i]
    if (length(candidates) != 1L) {
      stop("Expected exactly one ", kind, " containing chr", sprintf("%02d", i),
           " in ", directory, "; found ", length(candidates))
    }
    answer[[i]] <- candidates
  }
  answer
}

alignment_files <- find_one_per_chr(
  alignment_dir, "\\.(phy|phylip|txt|aln)$", "alignment"
)
tree_files <- if (dataset == "p28") rep(normalizePath(Sys.getenv("P28_TREE"), mustWork=TRUE),9) else find_one_per_chr(tree_dir, "\\.treefile$", "tree")

read_phylip_header <- function(path) {
  first <- readLines(path, n = 1L, warn = FALSE)
  fields <- strsplit(trimws(first), "[[:space:]]+")[[1]]
  values <- suppressWarnings(as.integer(fields[1:2]))
  if (length(values) != 2L || anyNA(values) || any(values <= 0)) {
    stop("Invalid PHYLIP header in ", path)
  }
  values
}

parent_of_tip <- function(tree, label) {
  tip <- match(label, tree$tip.label)
  rows <- which(tree$edge[, 2] == tip)
  if (length(rows) != 1L) stop("Cannot identify parent of ", label)
  tree$edge[rows, 1]
}

calibration_label <- function(lower_ma, upper_ma) {
  # Single quotes protect the comma-containing B() label in Newick.
  sprintf("'B(%.6f,%.6f)'", lower_ma / 100, upper_ma / 100)
}

write_tree_file <- function(tree, path, dryas_lower_ma, dryas_upper_ma,
                            calibrate_dryas = TRUE) {
  ntip <- length(tree$tip.label)
  tree$node.label <- rep("", tree$Nnode)
  root_node <- ntip + 1L
  dryas_node <- ape::getMRCA(tree, dryas_tips)
  cerc_node <- parent_of_tip(tree, "Cermon_DW847")
  if (is.null(dryas_node) || cerc_node == root_node) {
    stop("Could not identify calibration nodes")
  }

  # ape sanitizes punctuation inside node labels. Write safe placeholders first,
  # then replace them in the finished Newick string with MCMCTree B() labels.
  tree$node.label[root_node - ntip] <- "CAL_ROOT"
  if (calibrate_dryas) {
    tree$node.label[dryas_node - ntip] <- "CAL_DRYAS"
  }
  tree$node.label[cerc_node - ntip] <- "CAL_CERCOCARPUS"
  core_node <- ape::getMRCA(tree, core_tips)
  aih_node <- ape::getMRCA(tree, aih_tips)
  japan_node <- if (scenario == "japan_stem") dryas_node else core_node
  if (use_fossils) {
    if (!ape::is.monophyletic(tree, core_tips) || !ape::is.monophyletic(tree, aih_tips))
      stop("Core Dryas or AIH tip set is not monophyletic; inspect topology")
    tree$node.label[japan_node - ntip] <- "CAL_JAPAN"
    tree$node.label[aih_node - ntip] <- "CAL_AIH"
  }
  tree$edge.length <- NULL
  newick <- ape::write.tree(tree)
  newick <- gsub(
    "CAL_ROOT",
    calibration_label(root_min_ma, root_max_ma),
    newick, fixed = TRUE
  )
  if (calibrate_dryas) {
    newick <- gsub(
      "CAL_DRYAS",
      calibration_label(dryas_lower_ma, dryas_upper_ma),
      newick, fixed = TRUE
    )
  }
  newick <- gsub(
    "CAL_CERCOCARPUS",
    sprintf("'L(%.6f)'", cerc_min_ma / 100),
    newick, fixed = TRUE
  )
  newick <- gsub("CAL_JAPAN", "'L(0.060000)'", newick, fixed=TRUE)
  newick <- gsub("CAL_AIH", "'L(0.034000)'", newick, fixed=TRUE)
  if (grepl("CAL_", newick, fixed = TRUE) || grepl("'B-", newick, fixed = TRUE)) {
    stop("Internal error while writing MCMCTree calibration labels")
  }
  writeLines(c(sprintf("%d 1", ntip), newick), path)
  invisible(c(root = root_node, dryas = dryas_node, cercocarpus = cerc_node, core = core_node, aih = aih_node, japan = japan_node))
}

common_ctl <- function(seqfile, treefile, usedata, seed, outfile, mcmcfile) {
  c(
    paste("seed =", seed),
    paste("seqfile =", seqfile),
    paste("treefile =", treefile),
    paste("mcmcfile =", mcmcfile),
    paste("outfile =", outfile),
    "ndata = 1",
    "seqtype = 0",
    if (usedata == 2) "usedata = 2 in.BV" else paste("usedata =", usedata),
    "clock = 2",
    "RootAge = <0.487",
    "model = 4",
    "alpha = 0.5",
    "ncatG = 5",
    "cleandata = 0",
    # PAML >=4.10.10 requires C (conditional) or M (multiplicative).
    "BDparas = 1 1 0.1 C",
    "kappa_gamma = 6 2",
    "alpha_gamma = 1 1",
    # Mean overall rate = 2/40 = 0.05 substitutions/site/100 Ma.
    "rgene_gamma = 2 40 1",
    "sigma2_gamma = 1 10 1",
    "print = 1",
    "burnin = 100000",
    "sampfreq = 20",
    "nsample = 100000"
  )
}

manifest <- list()
for (chr in 1:9) {
  chr_name <- sprintf("chr%02d", chr)
  chr_dir <- file.path(output_dir, chr_name)
  prepare_dir <- file.path(chr_dir, "prepare")
  chain1_dir <- file.path(chr_dir, "chain1")
  chain2_dir <- file.path(chr_dir, "chain2")
  dir.create(prepare_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(chain1_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(chain2_dir, recursive = TRUE, showWarnings = FALSE)

  header <- read_phylip_header(alignment_files[[chr]])
  if (header[[1]] != length(expected_tips)) {
    stop(alignment_files[[chr]], " declares ", header[[1]],
         " sequences; expected ", length(expected_tips))
  }
  alignment_out <- file.path(chr_dir, "alignment.phy")
  file.copy(alignment_files[[chr]], alignment_out, overwrite = TRUE)

  tree <- ape::read.tree(tree_files[[chr]])
  if (inherits(tree, "multiPhylo")) tree <- tree[[1]]
  if (!setequal(tree$tip.label, expected_tips)) {
    stop(tree_files[[chr]], " does not contain the expected sample labels")
  }
  tree <- ape::root(tree, outgroup = outgroup_tips, resolve.root = TRUE)
  if (!ape::is.binary(tree) || !ape::is.monophyletic(tree, dryas_tips)) {
    stop("Required clades are not monophyletic in ", tree_files[[chr]])
  }
  tree_out <- file.path(chr_dir, "calibrated_tree.tre")
  nodes <- write_tree_file(
    tree, tree_out, bpp_dryas_lower_ma[[chr]], bpp_dryas_upper_ma[[chr]],
    calibrate_dryas = use_bpp_crown_dryas
  )

  writeLines(
    common_ctl("../alignment.phy", "../calibrated_tree.tre", 3, -1,
               "prepare.out", "prepare.mcmc.txt"),
    file.path(prepare_dir, "mcmctree.ctl")
  )
  for (chain in 1:2) {
    chain_dir <- if (chain == 1) chain1_dir else chain2_dir
    seed <- 20261002L + chr * 100L + chain
    writeLines(
      common_ctl("../alignment.phy", "../calibrated_tree.tre", if (prior_only) 0 else 2, seed,
                 "mcmctree.out", "mcmc.txt"),
      file.path(chain_dir, "mcmctree.ctl")
    )
  }

  manifest[[chr]] <- data.frame(
    scenario = scenario, prior_only = prior_only,
    core_node = nodes[["core"]], aih_node = nodes[["aih"]],
    japan_node = if (use_fossils) nodes[["japan"]] else NA_integer_,
    japan_min_ma = if (use_fossils) 6 else NA_real_,
    aih_min_ma = if (use_fossils) 3.4 else NA_real_,
    chromosome = chr, alignment_source = alignment_files[[chr]],
    tree_source = tree_files[[chr]], alignment_sites = header[[2]],
    crown_dryas_calibration = if (use_bpp_crown_dryas) "BPP_interval" else "none",
    crown_dryas_mean_ma = if (use_bpp_crown_dryas) bpp_dryas_mean_ma[[chr]] else NA_real_,
    crown_dryas_lower_ma = if (use_bpp_crown_dryas) bpp_dryas_lower_ma[[chr]] else NA_real_,
    crown_dryas_upper_ma = if (use_bpp_crown_dryas) bpp_dryas_upper_ma[[chr]] else NA_real_,
    root_min_ma = root_min_ma, root_max_ma = root_max_ma,
    cercocarpus_min_ma = cerc_min_ma,
    root_node = nodes[["root"]], dryas_node = nodes[["dryas"]],
    cercocarpus_node = nodes[["cercocarpus"]], stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, manifest)
utils::write.table(manifest, file.path(output_dir, "mcmctree_manifest.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE)

message("Prepared nine chromosome analyses in: ", output_dir)
message("Crown-Dryas calibration: ",
        if (use_bpp_crown_dryas) "chromosome-specific BPP interval" else "none")
message("Next: bash run_mcmctree_dryas.sh ", output_dir)
