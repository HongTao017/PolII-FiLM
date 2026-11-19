suppressPackageStartupMessages(library(motifStack))
suppressPackageStartupMessages(library(jsonlite))
root <- "/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019"
source <- "/mnt/HDD8TB/hongtao/Data_Seq/task/motif_pattern_families_20261018/source_data/patterns.json.gz"
data <- fromJSON(gzfile(source, "rb"), simplifyVector=FALSE)
regions <- c(tss_sequence="TSS", first_donor_sequence="Donor", first_acceptor_sequence="Acceptor", pre_tes_sequence="preTES")
ids <- vapply(data$patterns, function(p) {
  number <- as.integer(sub(".*pattern_", "", p$pattern_id))
  sprintf("F%02d-%s-%s%02d", p$fold, regions[p$region], ifelse(p$sign=="pos_patterns","P","N"), number)
}, "")
pfms <- setNames(lapply(seq_along(data$patterns), function(i) {
  mat <- do.call(cbind, lapply(data$patterns[[i]]$ppm, unlist))
  rownames(mat) <- c("A", "C", "G", "T")
  stopifnot(all(dim(mat)==c(4L,50L)), max(abs(colSums(mat)-1))<1e-5)
  new("pfm", name=ids[i], mat=mat)
}), ids)
saveRDS(list(pfms=pfms, data=data, ids=ids), file.path(root,"input_matrices.rds"))
cat("Benchmark 8 original 50bp PFMs, defaults\n")
timing <- system.time(alignment <- matalign(pfms[1:8]))
print(timing)
print(head(alignment))
cat("Expected all-pair seconds:", unname(timing["elapsed"])*choose(314,2)/choose(8,2), "\n")
saveRDS(alignment, file.path(root,"benchmark_alignment.rds"))
