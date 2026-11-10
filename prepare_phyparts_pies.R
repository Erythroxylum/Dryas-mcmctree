#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ape))
args <- commandArgs(trailingOnly=TRUE)
if (!length(args)) stop("Usage: Rscript prepare_phyparts_pies.R PREFIX [OUTPUT_CSV]")
prefix <- args[1]
output <- if (length(args)>=2) args[2] else "data/phyparts_p28_pies.csv"
clade_key <- function(labels) paste(sort(labels), collapse=";")
clade_counts <- function(tree) {
  n <- Ntip(tree)
  children <- split(tree$edge[,2],tree$edge[,1])
  desc <- function(z) {
    if (z<=n) return(tree$tip.label[z])
    unlist(lapply(children[[as.character(z)]],desc),use.names=FALSE)
  }
  nodes <- seq.int(n+1,n+tree$Nnode)
  setNames(suppressWarnings(as.numeric(tree$node.label)),
    vapply(nodes,function(z) clade_key(desc(z)),""))
}
lines <- readLines(paste0(prefix,".concon.tre"),warn=FALSE)
lines <- lines[nzchar(trimws(lines))]
if (length(lines)<2) stop("Missing concordance/conflict trees.")
concord <- clade_counts(read.tree(text=lines[1]))
conflict <- clade_counts(read.tree(text=lines[2]))
key_lines <- readLines(paste0(prefix,".node.key"),warn=FALSE)
key_lines <- key_lines[nzchar(trimws(key_lines))]
ids <- sub("^([0-9]+).*","\\1",key_lines)
subtrees <- sub("^[0-9]+[[:space:]]+","",key_lines)
keys <- setNames(vapply(subtrees,function(s) {
  if (!grepl(";[[:space:]]*$",s)) s <- paste0(s,";")
  clade_key(read.tree(text=s)$tip.label)
},""),ids)
hist <- readLines(paste0(prefix,".hist"),warn=FALSE)
hist <- hist[nzchar(trimws(hist))]
rows <- lapply(strsplit(hist,",",fixed=TRUE),function(a) {
  id <- sub("^Node","",a[1])
  v <- as.numeric(a[-1]); total <- tail(v,1); c <- v[1]
  alternatives <- if (length(v)>2) v[2:(length(v)-1)] else numeric()
  m <- max(c(0,alternatives))
  clade <- unname(keys[id]); d <- unname(conflict[clade])
  if (is.na(clade) || is.na(d) || is.na(concord[clade]) ||
      concord[clade]!=c) stop("Missing/mismatched clade: ",id)
  data.frame(phyparts_node=as.integer(id),descendants=clade,total=total,
    concordant=c,top_conflict=m,other_conflict=d-m,
    uninformative=total-c-d,
    plot=length(strsplit(clade,";",fixed=TRUE)[[1]])<27)
})
result <- do.call(rbind,rows)
if (nrow(result)!=26 || anyDuplicated(result$descendants) ||
    any(result$total!=9) ||
    any(as.matrix(result[c("concordant","top_conflict","other_conflict","uninformative")])<0))
  stop("Unexpected clade counts or invalid proportions; inspect PhyParts outputs.")
if (any(result$uninformative!=0))
  stop("Uninformative signal requires a fourth pie category; no three-category table written.")
dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
write.csv(result,output,row.names=FALSE)
message("Wrote ",nrow(result)," clades to ",output)
