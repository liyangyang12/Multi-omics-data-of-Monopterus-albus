ibrary(DiffBind)
tamoxifen <- dba(sampleSheet="SampleSheet7.csv")
tamoxifen <- dba.count(tamoxifen,bUseSummarizeOverlaps=TRUE,bParallel = TRUE)
tamoxifen <- dba.normalize(tamoxifen)
tamoxifen <- dba.contrast(tamoxifen,categories=DBA_TREATMENT,minMembers = 2)
tamoxifen <- dba.analyze(tamoxifen, method=DBA_DESEQ2)
tamoxifen.DB <- dba.report(tamoxifen)
show <-dba.show(tamoxifen,bContrasts=T)
show <- as.data.frame(show)
for (i in 1:nrow(show)){
	tamoxifen.DB <- dba.report(tamoxifen,contrast = i)
        result <- as.data.frame(tamoxifen.DB)
	write.table(result, file=paste(show[i,2],"_vs_",show[i,4],".txt",sep=""),sep="\t",quote=F,row.names=FALSE)}



