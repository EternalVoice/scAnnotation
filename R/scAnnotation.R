#' @title scAnnotation
#' @param seurat_obj
#' @param annot.method dataset or marker
#' @param dataset_ref If annot.method = 'dataset', you should
#' provide a SingleR reference dataset in `.rds` format.
#' @param db.type scType, CellMarker, CellTaxonomy, snRandom.
#' Only CellMarker and CellTaxonomy support Mouse.
#' @param Tissue
#' @param Species
#' @param Cancer
#' 
#' @importFrom SingleR SingleR
#' @import Seurat
#' 
#' @examples 
#' 
#' @export
#' 
scAnnotation <- function(seurat_obj,
                         annot.method = 'marker',
                         dataset_ref = NULL,
                         db.type = 'CellMarker',
                         Tissue,
                         Species = 'Human',
                         Cancer = TRUE) {
  if (annot.method == 'dataset') {
    REF <- readRDS(dataset_ref)
    pred <- SingleR::SingleR(
      test = seurat_obj@assays$RNA@data,
      ref = REF,
      labels = REF$label.fine,
      clusters = seurat_obj$seurat_clusters
    )
    pred_df <- data.frame(pred)
    cell.annot <- data.frame(ClusterID = rownames(pred_df), SingleRAnnot = pred_df$labels)
    return(cell.annot)
  } else if (annot.method == 'marker') {
    cell.annot <- cellAnnot(
      seurat_obj = seurat_obj, db.type = db.type, Tissue = Tissue, 
      Species = Species, Cancer = Cancer
    )
  }
  return(cell.annot)
}
