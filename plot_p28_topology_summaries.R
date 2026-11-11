#!/usr/bin/env Rscript
# Equal chromosome weights; same full-taxon rooted-bipartition counts as PhyParts.
suppressPackageStartupMessages(library(ape))
args <- commandArgs(trailingOnly=TRUE)
if (!length(args)) stop("Usage: Rscript plot_p28_topology_summaries.R ROOTED_TREE_DIR [OUTPUT_DIR]")
input <- args[1]
output <- if (length(args)>1L) args[2] else "figures/topology_summaries"
dir.create(output,recursive=TRUE,showWarnings=FALSE)
files <- list.files(input,pattern="\\.(tre|tree|nwk|txt)$",full.names=TRUE)
chr <- vapply(basename(files),function(s) {
  m <- regmatches(s,regexec("(chr|ch)0*([1-9])([._]|$)",s))[[1]]
  if (length(m)) as.integer(m[3]) else NA_integer_
},0L)
files <- files[!is.na(chr)]; chr <- chr[!is.na(chr)]
if (length(files)!=9L || anyDuplicated(chr) || !setequal(chr,1:9))
  stop("Require one tree per chromosome 1-9; only individual trees belong in this directory.")
files <- files[order(chr)]
trees <- lapply(files,function(f) {
  tr <- read.tree(f)
  if (inherits(tr,"multiPhylo")) stop("Expected exactly one tree in ",f)
  if (Ntip(tr)!=28L || anyDuplicated(tr$tip.label) || !"drum" %in% tr$tip.label)
    stop("Expected 28 unique populations including drum: ",f)
  root <- setdiff(tr$edge[,1],tr$edge[,2])
  children <- tr$edge[tr$edge[,1]==root,2]
  if (length(children)!=2L || !match("drum",tr$tip.label) %in% children)
    stop("Tree must already be rooted with drum as one root daughter: ",f)
  tr
})
if (!all(vapply(trees,function(tr) setequal(tr$tip.label,trees[[1]]$tip.label),TRUE)))
  stop("Chromosome taxon sets differ.")
class(trees) <- "multiPhylo"
clade_key <- function(tips) paste(sort(tips),collapse=";")
tree_clades <- function(tr) {
  n <- Ntip(tr); children <- split(tr$edge[,2],tr$edge[,1])
  desc <- function(z) if (z<=n) tr$tip.label[z] else
    unlist(lapply(children[[as.character(z)]],desc),use.names=FALSE)
  nodes <- seq.int(n+1,n+tr$Nnode)
  data.frame(node=nodes,descendants=vapply(nodes,function(z) clade_key(desc(z)),""),
    size=vapply(nodes,function(z) length(desc(z)),0L))
}
chrom_clades <- lapply(trees,tree_clades)
all_keys <- unlist(lapply(chrom_clades,function(d) d$descendants[d$size<28]),use.names=FALSE)
frequency <- table(all_keys)
sets <- lapply(names(frequency),function(s) strsplit(s,";",fixed=TRUE)[[1]])
incompatible <- function(a,b) length(intersect(a,b))>0L &&
  length(setdiff(a,b))>0L && length(setdiff(b,a))>0L
counts_for <- function(tr) {
  d <- tree_clades(tr); root <- setdiff(tr$edge[,1],tr$edge[,2])
  d <- d[d$node!=root,,drop=FALSE]
  values <- lapply(d$descendants,function(s) {
    a <- strsplit(s,";",fixed=TRUE)[[1]]
    c <- sum(vapply(chrom_clades,function(x) s %in% x$descendants,TRUE))
    conflict_keys <- names(frequency)[vapply(sets,function(b) incompatible(a,b),TRUE)]
    # Count distinct trees in conflict, not the sum of alternative clade frequencies.
    conflict <- sum(vapply(chrom_clades,function(x)
      any(x$descendants %in% conflict_keys),TRUE))
    top <- if (length(conflict_keys)) max(frequency[conflict_keys]) else 0L
    c(concordant=c,top_conflict=top,other_conflict=conflict-top,
      uninformative=9-c-conflict)
  })
  result <- cbind(d,as.data.frame(do.call(rbind,values)))
  result$total <- 9L
  result$plot <- result$size<27L
  if (any(as.matrix(result[c("concordant","top_conflict","other_conflict","uninformative")])<0) ||
      any(result$uninformative!=0)) stop("Unexpected counts: unresolved input signal requires a fourth slice.")
  result
}
# These trees are topology summaries; supplied CASTER branch lengths are not dates.
references <- list(chromosome1=trees[[1]],
  majority_rule=consensus(trees,p=0.5,check.labels=TRUE,rooted=TRUE))
key <- if (file.exists("population_species_key.csv"))
  read.csv("population_species_key.csv",stringsAsFactors=FALSE) else NULL
pie_colors <- c("black","white","#999999")
draw_tree <- function(tr,counts,title) {
  n <- Ntip(tr); children <- split(tr$edge[,2],tr$edge[,1])
  root <- setdiff(tr$edge[,1],tr$edge[,2])
  x <- y <- numeric(n+tr$Nnode)
  leaf_order <- function(z) if (z<=n) z else unlist(lapply(children[[as.character(z)]],leaf_order))
  order <- leaf_order(root); y[order] <- rev(seq_len(n))
  height <- numeric(length(x))
  position <- function(z) {
    if (z<=n) return(invisible(NULL))
    ch <- children[[as.character(z)]]
    for (a in ch) position(a)
    height[z] <<- 1+max(height[ch])
    y[z] <<- mean(range(y[leaf_order(z)]))
  }
  position(root); x <- height[root]-height
  labels <- tr$tip.label
  if (!is.null(key)) {
    ix <- match(labels,key$tree_tip)
    if (anyNA(ix)) stop("Population key does not cover all tree tips.")
    labels <- key$population[ix]
  }
  par(mar=c(1,1,3,1),xpd=FALSE)
  plot.new(); plot.window(xlim=c(-0.5,height[root]+4),ylim=c(0,n+1),xaxs="i",yaxs="i")
  title(main=title,cex.main=1)
  for (z in unique(tr$edge[,1])) {
    ch <- children[[as.character(z)]]
    segments(x[z],min(y[ch]),x[z],max(y[ch]),lwd=0.8)
    segments(x[z],y[ch],x[ch],y[ch],lwd=0.8)
  }
  text(x[seq_len(n)]+0.15,y[seq_len(n)],gsub("_"," ",labels),adj=0,cex=0.75)
  rx <- xinch(0.075); ry <- yinch(0.075)
  for (j in which(counts$plot)) {
    z <- counts$node[j]
    values <- unlist(counts[j,c("concordant","top_conflict","other_conflict")],use.names=FALSE)
    bounds <- c(0,cumsum(values/9))*2*pi+pi/2
    for (k in seq_along(values)) {
      if (!values[k]) next
      a <- seq(bounds[k],bounds[k+1],length.out=80)
      polygon(c(x[z],x[z]+rx*cos(a)),c(y[z],y[z]+ry*sin(a)),
              col=pie_colors[k],border=NA)
    }
    a <- seq(0,2*pi,length.out=100)
    lines(x[z]+rx*cos(a),y[z]+ry*sin(a),lwd=0.6)
  }
}
for (name in names(references)) {
  tr <- references[[name]]
  counts <- counts_for(tr)
  tr$node.label <- NULL; tr$edge.length <- NULL; tr$root.edge <- NULL
  prefix <- file.path(output,paste0("p28_",name))
  write.tree(tr,paste0(prefix,".nwk"))
  write.csv(counts,paste0(prefix,"_pies.csv"),row.names=FALSE)
  title <- if (name=="chromosome1") "Chromosome 1 topology" else
    "Majority-rule consensus of nine chromosome trees"
  pdf(paste0(prefix,".pdf"),width=9,height=11,useDingbats=FALSE)
  draw_tree(tr,counts,title); dev.off()
  png(paste0(prefix,".png"),width=2700,height=3300,res=300)
  draw_tree(tr,counts,title); dev.off()
}
message("Wrote chromosome 1 and majority-rule topology figures to ",output)
