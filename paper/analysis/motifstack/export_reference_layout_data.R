suppressPackageStartupMessages(library(motifStack))
suppressPackageStartupMessages(library(jsonlite))
root <- "/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019"
x <- readRDS(file.path(root,"input_matrices.rds"))
hc <- readRDS(file.path(root,"analysis/motifStack_hclust.rds"))
# Cut affects only which major branch gets an outer label. Keep the full tree.
groups <- cutree(hc,k=15)
sizes <- sort(table(groups),decreasing=TRUE)
cat("Display groups and native-tree order\n")
patterns <- lapply(seq_along(x$ids),function(i){
  p<-x$data$patterns[[i]]
  list(id=x$ids[i],key=p$key,group=unname(groups[x$ids[i]]),
       positive=p$sign=="pos_patterns",matched=length(p$matches)>0,
       ppm=p$ppm)
})
for(g in names(sizes)){
  members <- which(groups[x$ids]==as.integer(g))
  positions <- match(members,hc$order)
  cat(g,length(members),"order",min(positions),max(positions),
      "matched",sum(vapply(x$data$patterns[members],function(p)length(p$matches)>0,TRUE)),"\n")
}
write_json(list(patterns=patterns,labels=hc$labels,merge=unname(hc$merge),
               heights=hc$height,order=hc$order,display_k=15,
               matrix_method="motifStack defaults, unchanged"),
           file.path(root,"analysis/reference_layout_data.json"),auto_unbox=TRUE,pretty=FALSE)
