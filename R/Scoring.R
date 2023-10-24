#' @title Scoring
#' @param expr input scRNA-seq matrix (rounames - genes, column names - cells)
#' @param gs list of gene sets positively expressed in the cell type
#' 
#' @importFrom scales rescale
#' 
#' @export
#' 
Scoring <- function(expr,gs){
  # marker sensitivity
  marker_stat <- sort(table(unlist(gs)), decreasing = T)
  marker_sensitivity <- data.frame(
    score_marker_sensitivity = scales::rescale(as.numeric(marker_stat),to = c(0,1),from = c(length(gs),1)),
    gene = names(marker_stat),
    stringsAsFactors = FALSE
  )
  # subselect genes only found in dataset
  names_gs_cp <- names(gs)
  gs <- lapply(1:length(gs), function(i){
    GeneIndToKeep <- rownames(expr) %in% as.character(gs[[i]])
    rownames(expr)[GeneIndToKeep]
  })
  names(gs) <- names_gs_cp
  cell_markers_genescore <- marker_sensitivity[marker_sensitivity$gene %in% unique(unlist(gs)),]
  
  # multiple by marker sensitivity
  for (j in 1:nrow(cell_markers_genescore)) {
    expr[cell_markers_genescore[j,'gene'],] <- expr[cell_markers_genescore[j,'gene'],] * cell_markers_genescore[j,'score_marker_sensitivity']
  }
  # subselect only with marker genes
  expr <- expr[unique(c(unlist(gs))),]
  # combine scores
  es <- do.call('rbind', lapply(names(gs), function(i){
    sapply(1:ncol(expr), function(j){
      gs_z <- expr[gs[[i]], j]
      sum_t1 <- (sum(gs_z)/sqrt(length(gs_z)))
      sum_t1
    })
  }))
  dimnames(es) <- list(names(gs), colnames(expr))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),]
  es.max
}
