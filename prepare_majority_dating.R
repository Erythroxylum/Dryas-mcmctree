#!/usr/bin/env Rscript
# Build an explicitly resolved majority-rule p28 topology for fossil-only dating.
# Usage: Rscript prepare_majority_dating.R ROOTED_TREE_DIR BASE32 REPS MEMBERSHIP OUTPUT_DIR
suppressPackageStartupMessages(library(ape))
args <- commandArgs(TRUE)
if (length(args)!=5L) stop("Need ROOTED_TREE_DIR BASE32 REPS MEMBERSHIP OUTPUT_DIR")
out <- args[5]
paths_out <- file.path(out,c("p28_majority_unresolved.populations.nwk",
  "p28_majority_resolved.populations.nwk","p28_majority_resolved.samples32.tre",
  "p28_majority_clade_audit.tsv","p28_majority_resolution_choices.tsv"))
if (any(file.exists(paths_out))) stop("Topology outputs already exist; use a new output directory.")
files <- list.files(args[1],pattern="\\.(tre|tree|nwk|txt)$",full.names=TRUE)
chr <- vapply(basename(files),function(s) {
  m <- regmatches(s,regexec("(chr|ch)0*([1-9])([._]|$)",s))[[1]]
  if (length(m)) as.integer(m[3]) else NA_integer_
},0L)
files <- files[!is.na(chr)];chr <- chr[!is.na(chr)]
if (length(files)!=9L || anyDuplicated(chr) || !setequal(chr,1:9))
  stop("Need exactly one rooted tree for chromosomes 1-9.")
files <- files[order(chr)]
trees <- lapply(files,function(f) {
  tr <- read.tree(f)
  if (inherits(tr,"multiPhylo")) stop("One tree required per file: ",f)
  tr$node.label <- NULL
  tr <- collapse.singles(tr)
  if (Ntip(tr)!=28L || anyDuplicated(tr$tip.label) || !"drum" %in% tr$tip.label)
    stop("Need 28 unique population tips, including drum.")
  r <- setdiff(tr$edge[,1],tr$edge[,2]);ch <- tr$edge[tr$edge[,1]==r,2]
  if (length(ch)!=2 || !match("drum",tr$tip.label) %in% ch)
    stop("Input trees must be rooted on drum.")
  tr
})
tips <- trees[[1]]$tip.label
if (!all(vapply(trees,function(tr) setequal(tr$tip.label,tips),TRUE)))
  stop("Chromosome taxon sets differ.")
clade_key <- function(tips) paste(sort(tips,method="radix"),collapse=";")
desc_function <- function(tr) {
  n <- Ntip(tr);ch <- split(tr$edge[,2],tr$edge[,1])
  desc <- function(z) if (z<=n) tr$tip.label[z] else
    unlist(lapply(ch[[as.character(z)]],desc),use.names=FALSE)
  desc
}
clades <- lapply(trees,function(tr) {
  desc <- desc_function(tr);n <- Ntip(tr)
  sets <- lapply(seq.int(n+1,n+tr$Nnode),desc)
  unique(vapply(sets[lengths(sets)>=2L & lengths(sets)<28L],clade_key,""))
})
votes <- table(unlist(clades,use.names=FALSE))
sets <- setNames(lapply(names(votes),function(s) strsplit(s,";",fixed=TRUE)[[1]]),names(votes))
conflicts <- function(a,b) length(intersect(a,b))>0 &&
  length(setdiff(a,b))>0 && length(setdiff(b,a))>0
compatible <- function(k,selected)
  !any(vapply(selected,function(z) conflicts(sets[[k]],sets[[z]]),TRUE))
majority <- names(votes)[as.integer(votes)>=5L]
selected <- majority
choices <- list()
# Deterministic greedy extension: frequency descending, then lexical clade key.
candidates <- names(votes)[order(-as.integer(votes),names(votes),method="radix")]
for (k in candidates) {
  if (length(selected)==26L) break
  if (k %in% selected || !compatible(k,selected)) next
  available <- candidates[!candidates %in% selected]
  tied <- available[as.integer(votes[available])==as.integer(votes[k])]
  tied <- tied[vapply(tied,function(z) compatible(z,selected),TRUE)]
  choices[[length(choices)+1L]] <- data.frame(
    added_order=length(choices)+1L,votes=as.integer(votes[k]),
    proportion=as.integer(votes[k])/9,selected_clade=k,
    equally_frequent_compatible_candidates=paste(tied,collapse=" | "))
  selected <- c(selected,k)
}
if (length(selected)!=26L) stop("Observed clades cannot fully resolve majority tree; no unsupported resolutions added.")
if (!all(vapply(selected,function(k) compatible(k,setdiff(selected,k)),TRUE)))
  stop("Selected clades are incompatible.")
render <- function(keys) {
  clusters <- c(unname(sets[keys]),as.list(tips))
  build <- function(a) {
    if (length(a)==1L) return(a)
    sub <- clusters[vapply(clusters,function(b)
      length(b)<length(a) && all(b %in% a),TRUE)]
    children <- sub[!vapply(sub,function(b) any(vapply(sub,function(z)
      length(z)>length(b) && all(b %in% z),TRUE)),TRUE)]
    if (!setequal(unlist(children,use.names=FALSE),a)) stop("Invalid clade hierarchy.")
    children <- children[order(vapply(children,function(b) min(match(b,tips)),0L))]
    paste0("(",paste(vapply(children,build,""),collapse=","),")")
  }
  read.tree(text=paste0(build(tips),";"))
}
unresolved <- render(majority)
resolved <- render(selected)
if (!all(table(resolved$edge[,1])==2L)) stop("Resolved tree is not bifurcating.")

reps <- read.delim(args[3],stringsAsFactors=FALSE)
if (!all(c("sampleID","population") %in% names(reps)) || nrow(reps)!=28L ||
    anyDuplicated(reps$sampleID) || anyDuplicated(reps$population) ||
    !setequal(reps$population,tips)) stop("Representatives must map exactly the 28 original population IDs to samples.")
resolved_samples <- resolved
resolved_samples$tip.label <- reps$sampleID[match(resolved$tip.label,reps$population)]
outgroups <- c("Cermon_DW847","Chafol_1611","Pursta_DW848","Purtri_DW846b")
scaffold <- read.tree(args[2])
if (inherits(scaffold,"multiPhylo") || !setequal(scaffold$tip.label,c(reps$sampleID,outgroups)))
  stop("BASE32 must contain the 28 sample IDs and four expected outgroups.")
scaffold$node.label <- NULL
scaffold <- collapse.singles(root(scaffold,outgroup=outgroups,resolve.root=TRUE))
desc <- desc_function(scaffold)
dryas_node <- getMRCA(scaffold,reps$sampleID)
if (!setequal(desc(dryas_node),reps$sampleID)) stop("Dryas is not monophyletic in the outgroup scaffold.")
dryas_newick <- sub(";$","",write.tree(resolved_samples))
n <- Ntip(scaffold);ch <- split(scaffold$edge[,2],scaffold$edge[,1])
replace_dryas <- function(z) {
  if (z==dryas_node) return(dryas_newick)
  if (z<=n) return(scaffold$tip.label[z])
  paste0("(",paste(vapply(ch[[as.character(z)]],replace_dryas,""),collapse=","),")")
}
r <- setdiff(scaffold$edge[,1],scaffold$edge[,2])
full <- read.tree(text=paste0(replace_dryas(r),";"))
if (!all(table(full$edge[,1])==2L)) stop("Full scaffold tree is not bifurcating.")
# Check calibration groups by exact descendant sample sets.
membership <- read.delim(args[4],stringsAsFactors=FALSE)
if (!all(c("sampleID","calibration_clade") %in% names(membership)) ||
    anyDuplicated(membership$sampleID) ||
    !setequal(membership$sampleID,reps$sampleID) ||
    any(!membership$calibration_clade %in% c("drummondii","AIH","other_core")))
  stop("Invalid explicit calibration membership table.")
drum_sample <- reps$sampleID[match("drum",reps$population)]
if (!setequal(membership$sampleID[membership$calibration_clade=="drummondii"],drum_sample))
  stop("drummondii membership does not match drum.")
calibration_sets <- list(core_Dryas=membership$sampleID[membership$calibration_clade!="drummondii"],
                        AIH=membership$sampleID[membership$calibration_clade=="AIH"])
desc_full <- desc_function(full)
for (name in names(calibration_sets)) {
  a <- calibration_sets[[name]];z <- getMRCA(full,a)
  if (length(a)<2L || is.null(z) || !setequal(desc_full(z),a))
    stop("Fossil calibration clade is not retained: ",name)
}
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.tree(unresolved,paths_out[1])
write.tree(resolved,paths_out[2])
write.tree(full,paths_out[3])
audit <- data.frame(descendants=selected,chromosomes=as.integer(votes[selected]),
  fraction=as.integer(votes[selected])/9,
  origin=ifelse(selected %in% majority,"majority","greedy_resolution"))
write.table(audit,paths_out[4],sep="\t",row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,choices),paths_out[5],sep="\t",row.names=FALSE,quote=TRUE)
writeLines(c("Resolved majority-rule p28 dating reference.",
  "Retain all clades present in at least 5/9 chromosomes.",
  "Add remaining observed compatible clades by frequency, breaking ties by lexical clade key.",
  "This is not an MCC tree and need not be observed in a chromosome.",
  "The original 32-tip tree supplies only the outgroup scaffold; its Dryas subtree is replaced.",
  "Core Dryas and AIH fossil-calibration membership verified exactly.",
  paste("Original scaffold:",normalizePath(args[2])),
  paste("Chromosome inputs:",paste(normalizePath(files),collapse=";"))),
  file.path(out,"README_majority_resolution.txt"))
message("Prepared unresolved and resolved population trees, sample-labelled 32-tip input, and clade audit.")
print(do.call(rbind,choices)[,c("added_order","votes","selected_clade")],row.names=FALSE)
message("No MCMCTree run has been started.")
