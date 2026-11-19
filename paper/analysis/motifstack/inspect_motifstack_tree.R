suppressPackageStartupMessages(library(motifStack))
suppressPackageStartupMessages(library(jsonlite))
root <- "/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019"
hc <- readRDS(file.path(root,"analysis/motifStack_hclust.rds"))
cat("motifSignature defaults\n"); print(args(motifSignature))
cat("height quantiles\n");print(quantile(hc$height,c(0,.5,.75,.9,.95,.975,1)))
for(k in c(8,10,12,15,17,20)){
  cat("k=",k,"sizes:",sort(table(cutree(hc,k=k)),decreasing=TRUE),"\n")
}
phy <- ade4::hclust2phylog(hc)
write_json(list(labels=hc$labels,merge=unname(hc$merge),heights=hc$height,
                order=hc$order,leaves=names(phy$leaves),droot=as.list(phy$droot),
                parts=phy$parts),file.path(root,"analysis/tree_for_reference_layout.json"),
           auto_unbox=TRUE,pretty=FALSE)
cat("motifSignature source\n"); print(motifSignature)
