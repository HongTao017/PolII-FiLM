#!/usr/bin/env Rscript
# Native motifStack rendering. Only visual layout, font and colour are styled.
# All computation uses untrimmed original PPMs and package defaults.
suppressPackageStartupMessages(library(motifStack))
suppressPackageStartupMessages(library(Cairo))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(grid))
root <- "/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019"
out <- file.path(root,"figures")
dir.create(out,showWarnings=FALSE)
input <- readRDS(file.path(root,"input_matrices.rds"))
palette <- c(A="#68947E",C="#658CA9",G="#B29B67",T="#B5787C")
pfms <- lapply(input$pfms,function(p){p@color <- palette; p})
patterns <- setNames(input$data$patterns,input$ids)
references <- input$data$reference_motifs
CairoFonts(regular="Arial:style=Regular",bold="Arial:style=Bold",
           italic="Arial:style=Italic",bolditalic="Arial:style=Bold Italic")
open_pdf <- function(name) {
  CairoPDF(file=file.path(out,name),width=5,height=4,pointsize=12,family="Arial",onefile=TRUE)
}
header <- function(title,subtitle=NULL){
  grid.text(title,x=.5,y=.963,gp=gpar(fontfamily="Arial",fontsize=12,fontface="bold",col="black"))
  if(!is.null(subtitle))grid.text(subtitle,x=.5,y=.902,gp=gpar(fontfamily="Arial",fontsize=12,col="black"))
}
footer <- function(text)grid.text(text,x=.5,y=.035,gp=gpar(fontfamily="Arial",fontsize=12,col="black"))
ref_pfm <- function(m) {
  mat <- do.call(cbind,lapply(m$ppm,unlist))
  rownames(mat) <- c("A","C","G","T")
  new("pfm",mat=mat,name=paste0("REF-",m$matrix_id),color=palette)
}
save_alignment <- function(aligned, filename){
  saveRDS(aligned,file.path(root,"analysis",filename))
}

if("--smoke" %in% commandArgs(TRUE)){
  cat("Optional native-logo backends:\n")
  for(p in c("grImport2","grImport","Cairo"))cat(p,requireNamespace(p,quietly=TRUE),"\n")
  open_pdf("motifStack_render_check.pdf")
  sample <- DNAmotifAlignment(pfms[1:3])
  plotMotifLogoStack(sample,font="Arial",ncex=1,xaxis=FALSE,yaxis=FALSE,ylab=NA)
  dev.off()
  quit(status=0)
}

hc <- readRDS(file.path(root,"analysis/motifStack_hclust.rds"))
stopifnot(setequal(hc$labels,names(pfms)))
phylog <- ade4::hclust2phylog(hc)
leaves <- names(phylog$leaves)
leaf_ids <- input$ids[match(leaves,gsub("-","_",input$ids,fixed=TRUE))]
stopifnot(!anyNA(leaf_ids),!anyDuplicated(leaf_ids),setequal(leaf_ids,input$ids))
hit <- vapply(patterns[leaf_ids],function(p)length(p$matches)>0,TRUE)
stopifnot(sum(hit)==44L)

# Native radial renderer, with each real matched leaf highlighted grey.
# No cut height, number of groups, previous C label, or TF label determines the tree.
if(!("--skip-radial" %in% commandArgs(TRUE))){
open_pdf("motifStack_radial_314_default.pdf")
par(family="Arial",ps=12,cex=1,col="black",fg="black",lwd=1)
radial_pfms <- pfms[leaf_ids]
for(i in seq_along(radial_pfms))radial_pfms[[i]]@name <- leaves[i]
plotMotifStackWithRadialPhylog(phylog,pfms=radial_pfms,circle=1.22,
    labels.leaves=rep("",length(leaves)),clabel.leaves=1,cleaves=0,
    col.leaves=rep("black",length(leaves)),
    col.bg=ifelse(hit,"#E1E1E1","white"),
    circle.motif=1.29,col.inner.label.circle=ifelse(hit,"#777777","#D8DDE1"),
    inner.label.circle.width=.045,font="Arial")
text(0,2.24,"314 pattern PFMs",cex=1,font=2)
text(0,-2.25,"Grey: existing JASPAR matches",cex=1)
dev.off()
}

# All 314 PPMs in the global default-tree order. Alignment is local to each page,
# not a claim that four motifs are a statistically defined biological family.
ordered <- hc$labels[hc$order]
groups <- split(ordered,ceiling(seq_along(ordered)/4))
open_pdf("motifStack_all_314_PFM_matrix_atlas.pdf")
atlas_alignments <- list()
for(page in seq_along(groups)) {
  ids <- groups[[page]]
  q <- pfms[ids]
  aligned <- if(length(q)>1L)DNAmotifAlignment(q) else q
  atlas_alignments[[as.character(page)]] <- aligned
  for(j in seq_along(aligned)){
    original <- sub("\\(RC\\)$","",aligned[[j]]@name)
    reversed <- endsWith(aligned[[j]]@name,"(RC)")
    p <- patterns[[original]]
    aligned[[j]]@name <- paste0(original,if(reversed)" (RC)"else"", "  n=",p$seqlets)
  }
  grid.newpage()
  pushViewport(viewport(gp=gpar(fontfamily="Arial",fontsize=12,col="black")))
  header(sprintf("Pattern PFMs: page %d/%d",page,length(groups)))
  pushViewport(viewport(x=.5,y=.493,width=.98,height=.80))
  grid.draw(plotMotifLogoStack(aligned,font="Arial",ncex=1,xaxis=FALSE,yaxis=FALSE,
                    ylab=NA,draw=FALSE))
  popViewport()
  footer("Global tree order; default page alignment")
  popViewport()
}
dev.off()
save_alignment(atlas_alignments,"all_314_page_alignments.rds")

# For every existing reference, show up to three best qualified query matches.
# Native motifStack aligns the PFMs here; the displayed q remains from old TOMTOM.
# The reference enters only this display's small tree, never the 314-pattern tree.
match_records <- do.call(rbind,lapply(names(patterns),function(id){
  p <- patterns[[id]]
  if(!length(p$matches))return(NULL)
  do.call(rbind,lapply(p$matches,function(m){
    data.frame(query=id,TF=m$TF,matrix_id=m$matrix_id,q=m$q_value,
               p=m$p_value,overlap=m$overlap_bp,stringsAsFactors=FALSE)
  }))
}))
stopifnot(nrow(match_records)==79L, length(unique(match_records$query))==44L)
write.table(match_records,file.path(root,"analysis/original_TOMTOM_matches.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
open_pdf("motifStack_JASPAR_reference_matrix_panels.pdf")
shown <- list()
for(matrix_id in names(references)){
  reference <- references[[matrix_id]]
  matches <- match_records[match_records$matrix_id==matrix_id,]
  matches <- matches[order(matches$q,matches$p,matches$query),]
  matches <- head(matches[!duplicated(matches$query),],3)
  q <- c(pfms[matches$query],setNames(list(ref_pfm(reference)),paste0("REF-",matrix_id)))
  # motifStack's public wrapper handles two motifs as a plain aligned stack.
  # Its profile-distance hclust step requires at least three motifs.
  small_hc <- if(length(q)>2L)clusterMotifs(q) else NULL
  aligned <- DNAmotifAlignment(if(is.null(small_hc))q else q[small_hc$order])
  shown[[matrix_id]] <- list(matches=matches,hc=small_hc,aligned=aligned)
  for(j in seq_along(aligned)){
    id <- sub("\\(RC\\)$","",aligned[[j]]@name)
    rc <- endsWith(aligned[[j]]@name,"(RC)")
    if(startsWith(id,"REF-")){
      aligned[[j]]@name <- paste(reference$TF,matrix_id,if(rc)"(RC)"else"")
    }else{
      m <- matches[matches$query==id,][1,]
      aligned[[j]]@name <- paste0(id,if(rc)" RC"else"","  q=",formatC(m$q,format="g",digits=2))
    }
  }
  grid.newpage()
  pushViewport(viewport(gp=gpar(fontfamily="Arial",fontsize=12,col="black")))
  header(paste(reference$TF,matrix_id),"Query / JASPAR PFM comparison")
  pushViewport(viewport(x=.5,y=.468,width=.98,height=.735))
  if(is.null(small_hc)){
    grid.draw(plotMotifLogoStack(aligned,font="Arial",ncex=1,xaxis=FALSE,yaxis=FALSE,
                      ylab=NA,draw=FALSE))
  }else{
    grid.draw(plotMotifLogoStackWithTree(aligned,hc=small_hc,font="Arial",ncex=1,
         xaxis=FALSE,yaxis=FALSE,ylab=NA,draw=FALSE))
  }
  popViewport()
  footer("Default alignment; q from original TOMTOM")
  popViewport()
}
dev.off()
save_alignment(shown,"reference_panel_selection_and_alignment.rds")
write_json(list(output="Native motifStack PDF figures",canvas_inches=c(5,4),
                font="Arial",text_fontsize_pt=12,
                outputs=list(radial="motifStack_radial_314_default.pdf",
                matrix_atlas="motifStack_all_314_PFM_matrix_atlas.pdf",
                reference_panels="motifStack_JASPAR_reference_matrix_panels.pdf"),
                atlas_pages=length(groups),reference_pages=length(references),
                nucleotide_colours=as.list(palette),
                display_adjustments="page layout, label suppression in dense radial tree, matched-leaf grey highlight, font and colours only",
                no_positional_histogram="seqlet-coordinate distributions were not supplied to this redraw; none fabricated",
                sequence_logos="PPM information logos, not signed CWM contribution logos",
                previous_25_cluster_labels_reused=FALSE),
           file.path(root,"analysis/plot_receipt.json"),auto_unbox=TRUE,pretty=TRUE)
cat("Completed native motifStack figures\n")
