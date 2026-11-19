#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(motifStack))
suppressPackageStartupMessages(library(jsonlite))
root <- "/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019"
input <- readRDS(file.path(root,"input_matrices.rds"))
pfms <- input$pfms
dir.create(file.path(root,"pair_chunks"), showWarnings=FALSE)
dir.create(file.path(root,"analysis"), showWarnings=FALSE)
stopifnot(length(pfms)==314L, all(vapply(pfms,function(x)ncol(x@mat),1L)==50L))
pairs <- combn(names(pfms), 2, simplify=FALSE)
chunks <- split(seq_along(pairs), ceiling(seq_along(pairs)/128))
workers <- 24L
cat(format(Sys.time()), "Starting", length(pairs), "pairs; workers", workers, "\n")
flush.console()

# Each pair calls the public matalign() without changing any scientific option.
# This splits only independent pairs; reassembly retains the original combn order.
results <- parallel::mclapply(seq_along(chunks), function(k) {
  destination <- file.path(root,"pair_chunks",sprintf("chunk_%04d.rds",k))
  if(file.exists(destination)) return(destination)
  answer <- do.call(rbind,lapply(chunks[[k]], function(i) {
    matalign(pfms[pairs[[i]]])
  }))
  stopifnot(nrow(answer)==length(chunks[[k]]))
  saveRDS(answer, destination)
  destination
}, mc.cores=workers, mc.preschedule=TRUE, mc.set.seed=FALSE)
stopifnot(!any(vapply(results, inherits, TRUE, what="try-error")))
alignment <- do.call(rbind,lapply(results,readRDS))
rownames(alignment) <- NULL
stopifnot(nrow(alignment)==choose(314,2), all(is.finite(alignment$distance)))
stopifnot(identical(as.character(alignment$motif1),vapply(pairs,`[`,"",1L)),
          identical(as.character(alignment$motif2),vapply(pairs,`[`,"",2L)))

# Reproduce clusterMotifs() exactly: matalign() then motifHclust(method="average").
# Validation against the public wrapper on the eight-motif benchmark is exact.
benchmark <- readRDS(file.path(root,"benchmark_alignment.rds"))
subset <- alignment[alignment$motif1 %in% names(pfms)[1:8] & alignment$motif2 %in% names(pfms)[1:8],]
rownames(subset) <- NULL
stopifnot(isTRUE(all.equal(subset, benchmark, check.attributes=FALSE)))
reference_hc <- clusterMotifs(pfms[1:8])
checked_hc <- motifHclust(subset, method="average")
stopifnot(identical(reference_hc$merge, checked_hc$merge),
          identical(reference_hc$order, checked_hc$order),
          isTRUE(all.equal(reference_hc$height,checked_hc$height)))
hc <- motifHclust(alignment, method="average")
stopifnot(length(hc$order)==314L, setequal(hc$labels,names(pfms)),
          all(is.finite(hc$height)),all(diff(hc$height)>=0))
saveRDS(hc,file.path(root,"analysis/motifStack_hclust.rds"))
saveRDS(alignment,file.path(root,"analysis/all_pair_alignments.rds"))
write.table(alignment,file.path(root,"analysis/all_pair_alignments.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
d <- get("matalignOut2dist",asNamespace("motifStack"))(alignment)
write.table(d,file.path(root,"analysis/matalign_pair_distance_matrix.tsv"),sep="\t",quote=FALSE,col.names=NA)
write.table(as.matrix(dist(d)),file.path(root,"analysis/clusterMotifs_hclust_input_distances.tsv"),sep="\t",quote=FALSE,col.names=NA)
patterns <- input$data$patterns
metadata <- do.call(rbind,lapply(seq_along(patterns),function(i){
  p <- patterns[[i]]
  data.frame(short_id=input$ids[i],pattern_key=p$key,fold=p$fold,region=p$region,
             attribution_sign=p$sign,seqlet_memberships=p$seqlets,existing_matches=length(p$matches),
             existing_TF=paste(unique(vapply(p$matches,function(m)m$TF,"")),collapse=";"))
}))
metadata <- metadata[match(hc$labels[hc$order],metadata$short_id),]
metadata$tree_order <- seq_len(nrow(metadata))
write.table(metadata,file.path(root,"analysis/matrix_tree_order.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
receipt <- list(R=R.version.string,motifStack=as.character(packageVersion("motifStack")),
               input_matrix="314 original untrimmed 4x50 PPMs; A,C,G,T rows",
               method="Smith-Waterman",pseudo=1,revcomp=TRUE,
               pfm_to_pcm="package internal default pfm2pcm(total=1000), round(1000*PFM + 0.49)",
               clustering="motifStack::motifHclust(matalign_output, method='average')",
               package_distance_handling="motifHclust uses dist(matalignOut2dist(alignment)), unchanged",
               equivalence_to_clusterMotifs_wrapper="PASS: identical merge, order and heights on benchmark",
               number_patterns=314,number_pairs=nrow(alignment),workers=workers,
               original_trimming_used=FALSE,original_custom_similarity_used=FALSE,
               cluster_number_forced=FALSE,group_cutoff_used=FALSE,
               original_model_or_TOMTOM_rerun=FALSE,
               existing_matched_patterns=44,existing_qualified_TOMTOM_matches=79,
               logo_alignment_defaults=list(threshold=.4,minimalConsensus=0,revcomp=TRUE),
               complete=TRUE)
write_json(receipt,file.path(root,"analysis/parameters.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(root,"analysis/sessionInfo.txt"))
cat(format(Sys.time()), "Complete. Height range:", range(hc$height), "\n")
