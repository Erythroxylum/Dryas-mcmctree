#!/usr/bin/env Rscript
library(ape)
a <- commandArgs(TRUE); stopifnot(length(a)==3)
root<-a[1]; dataset<-a[2]
outgroups<-c("Cermon_DW847","Chafol_1611","Pursta_DW848","Purtri_DW846b")
# Use the numbered topology printed by PAML, not ape's node-number order.
numbered_tree <- function(path, tips, columns) {
 lines<-readLines(path,warn=FALSE)
 for (line in lines[grepl(";",lines,fixed=TRUE)]) {
  pos<-regexpr("(",line,fixed=TRUE)[1]; if(pos<1) next
  txt<-substring(line,pos)
  tr<-tryCatch(read.tree(text=txt),error=function(e) NULL)
  if(is.null(tr) || inherits(tr,"multiPhylo")) next
  tr$tip.label<-sub("^[0-9]+_", "",tr$tip.label)
  if(!setequal(tr$tip.label,tips) || is.null(tr$node.label)) next
  ids<-sub("^t_n", "",tr$node.label)
  if(any(!grepl("^[0-9]+$",ids)) || !all(paste0("t_n",ids) %in% columns)) next
  tr$node.label<-ids; return(tr)
 }
 stop("Could not recover PAML numbered topology from ",path,". No trees exported; inspect numbered tree format.")
}
alltrees<-list(); set.seed(20261104)
for(chr in 1:9) for(chain in 1:2) {
 d<-file.path(root,sprintf("chr%02d",chr)); cd<-file.path(d,paste0("chain",chain))
 x<-read.table(file.path(cd,"mcmc.txt"),header=TRUE,check.names=FALSE)
 template<-read.tree(text=paste(readLines(file.path(d,"calibrated_tree.tre"))[-1],collapse=""))
 tr<-numbered_tree(file.path(cd,"mcmctree.out"),template$tip.label,names(x))
 if(nrow(x)<50) stop("Fewer than 50 samples")
 for(i in sort(sample(seq_len(nrow(x)),50))) {
  ages<-c(rep(0,Ntip(tr)),as.numeric(x[i,paste0("t_n",tr$node.label)])*100)
  z<-tr; z$edge.length<-ages[z$edge[,1]]-ages[z$edge[,2]]
  if(any(!is.finite(z$edge.length)) || any(z$edge.length<0)) stop("Invalid node ages")
  z$node.label<-NULL; alltrees[[length(alltrees)+1]]<-z
 }
}
if(dataset=="p28") {m<-read.delim(a[3]); labels<-setNames(m$population,m$sampleID)} else {
 labels<-c(drum_BABY63="Dryas_drummondii",sp_POB10="Dryas_sp_Pobeda",cauc_ODP7438="Dryas_caucasica",oxyd_ODP43154="Dryas_oxyodonta",octo_END01="Dryas_subincisa",cren_ODP53249="Dryas_grandis",ajan_ODP48916="Dryas_ajanensis",alas_ODP43138="Dryas_alaskensis",hook_GLA03="Dryas_hookeriana",intg_ODP7653="Dryas_integrifolia",octo_JAT08="Dryas_octopetala")
}
if(anyDuplicated(labels)) stop("Duplicate display labels")
dryas<-lapply(alltrees,function(t) {t<-drop.tip(t,outgroups);t$tip.label<-unname(labels[t$tip.label]);if(anyNA(t$tip.label))stop("Missing label");t})
for(i in seq_along(alltrees)) {t<-alltrees[[i]]; idx<-t$tip.label %in% names(labels);t$tip.label[idx]<-unname(labels[t$tip.label[idx]]);alltrees[[i]]<-t}
class(alltrees)<-class(dryas)<-"multiPhylo"
dir.create(file.path(root,"exports"),showWarnings=FALSE)
write.nexus(alltrees,file=file.path(root,"exports/full.nex"))
write.nexus(dryas,file=file.path(root,"exports/dryas.nex"))
