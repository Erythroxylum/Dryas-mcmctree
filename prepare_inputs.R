#!/usr/bin/env Rscript
library(ape)
m<-read.delim("pop28_representatives.tsv",stringsAsFactors=FALSE)
stopifnot(all(c("sampleID","population") %in% names(m)),nrow(m)==28,!anyDuplicated(m$sampleID))
if(!file.exists("clade_membership.tsv")) {
 write.table(data.frame(sampleID=m$sampleID,population=m$population,clade="UNASSIGNED"),"clade_membership.tsv",sep="\t",quote=FALSE,row.names=FALSE)
 message("Created clade_membership.tsv: assign every row drummondii, AIH, or other_core before p28 setup.")
}
tips15<-c("drum_BABY63","sp_POB10","cauc_ODP7438","oxyd_ODP43154","octo_END01","cren_ODP53249","ajan_ODP48916","alas_ODP43138","hook_GLA03","intg_ODP7653","octo_JAT08")
out<-c("Cermon_DW847","Chafol_1611","Pursta_DW848","Purtri_DW846b")
files<-list.files("chromosome_alignments_trimmed",pattern="\\.(phy|phylip|txt|aln)$",full.names=TRUE)
if(length(files)!=9)stop("Expected nine chromosome PHYLIP files")
for(f in files) {
 x<-tryCatch(read.dna(f,format="sequential",as.character=TRUE),error=function(e) NULL)
 if(is.null(x)) x<-read.dna(f,format="interleaved",as.character=TRUE)
 for(dataset in c("15tip","p28")) {
  wanted<-c(if(dataset=="15tip") tips15 else m$sampleID,out)
  missing<-setdiff(wanted,rownames(x))
  if(length(missing))stop(dataset," missing samples in ",f,": ",paste(missing,collapse=", "),". Supply the corresponding original alignments; rows cannot be reconstructed.")
  dest<-file.path("prepared",dataset);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  write.dna(x[wanted,,drop=FALSE],file=file.path(dest,basename(f)),format="sequential",nbcol=-1,colsep="")
 }
}
