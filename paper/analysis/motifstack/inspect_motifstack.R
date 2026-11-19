suppressPackageStartupMessages(library(motifStack))
cat("R:", R.version.string, "\n")
cat("motifStack:", as.character(packageVersion("motifStack")), "\n")
for (name in c("clusterMotifs", "matalign", "motifHclust", "DNAmotifAlignment",
               "plotMotifStackWithRadialPhylog", "plotMotifLogo", "plotMotifLogoStack")) {
  cat("\nFUNCTION", name, "\n")
  print(args(get(name, asNamespace("motifStack"))))
  if (name %in% c("clusterMotifs", "motifHclust")) print(get(name, asNamespace("motifStack")))
}
cat("\nGRID DEFAULT COLOUR\n")
x <- new("pfm", name="test", mat=matrix(.25, 4, 10, dimnames=list(c("A","C","G","T"),NULL)))
print(x@color)
cat("\nFONTS\n")
system("fc-match Arial")
root <- "/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019"
hc <- readRDS(file.path(root,"analysis/motifStack_hclust.rds"))
phylog <- ade4::hclust2phylog(hc)
cat("\nHCLUST labels\n"); print(head(hc$labels))
cat("\nPHYLOG labels\n"); print(head(names(phylog$leaves)))
