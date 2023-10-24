#' @title cellAnnot
#' @param seurat_obj A Seurat object
#' @param db.type scType, CellMarker, CellTaxonomy, or snRandom (internal).
#' @param Species Human or Mouse.
#' @param Tissue 
#' @param Cancer TRUE or FALSE
#' 
#' @import Seurat
#' @importFrom dplyr `%>%` group_by top_n
#' @importFrom plyr mapvalues
#' 
#' @export
#' 
cellAnnot <- function(seurat_obj, 
                      db.type='CellMarker',
                      Tissue,
                      Species='Human', 
                      Cancer=TRUE) {
  
  if (db.type == 'scType') {
    gs <- prepareGeneset(db.type = db.type,Tissue = Tissue)
  } else if (db.type == 'CellMarker') {
    gs <- prepareGeneset(
      db.type = db.type, Species = Species,
      Tissue = Tissue, Cancer = Cancer
    )
  } else if (db.type == 'CellTaxonomy') {
    gs <- prepareGeneset(db.type = db.type, Species = Species,Tissue = Tissue)
  } else if (db.type == 'snRandom') {
    gs <- prepareGeneset(db.type = db.type)
  }
  
  DefaultAssay(seurat_obj) <- 'RNA'
  seurat_obj <- ScaleData(seurat_obj)
  es <- Scoring(expr = seurat_obj@assays$RNA@scale.data, gs = gs)
  result <- do.call('rbind',lapply(unique(seurat_obj$seurat_clusters),function(x){
    es.max.cl <- sort(rowSums(es[,rownames(seurat_obj@meta.data[seurat_obj$seurat_clusters == x,])]),decreasing = TRUE)
    head(data.frame(cluster = x, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj$seurat_clusters==x)),10)
  }))
  score <- result %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 1, wt = scores)
  
  if (db.type == 'scType') {
    seurat_obj$scTypeAnnot <- plyr::mapvalues(seurat_obj$seurat_clusters, from = score$cluster, to = score$type)
    annot.result <- seurat_obj@meta.data[,c('seurat_clusters','scTypeAnnot')]
  } else if (db.type == 'CellMarker') {
    seurat_obj$CellMarkerAnnot <- plyr::mapvalues(seurat_obj$seurat_clusters, from = score$cluster, to = score$type)
    annot.result <- seurat_obj@meta.data[, c('seurat_clusters','CellMarkerAnnot')]
  } else if (db.type == 'CellTaxonomy') {
    seurat_obj$CellTaxonomyAnnot <- plyr::mapvalues(seurat_obj$seurat_clusters, from = score$cluster, to = score$type)
    annot.result <- seurat_obj@meta.data[, c('seurat_clusters','CellTaxonomyAnnot')]
  } else if (db.type == 'snRandom') {
    seurat_obj$snRandomAnnot <- plyr::mapvalues(seurat_obj$seurat_clusters, from = score$cluster, to = score$type)
    annot.result <- seurat_obj@meta.data[, c('seurat_clusters','snRandomAnnot')]
  }
  
  return(annot.result)
}
