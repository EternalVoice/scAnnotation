#' @title prepareGeneset
#' @param db reference cell marker
#' @param db.type scType, CellMarker, CellTaxonomy, or snRandom-seq
#' @param Species Human or Mouse
#' @param Tissue
#' @param Cancer TRUE or FALSE
#' 
#' @export
#' 
prepareGeneset <- function(db.type,Tissue,Species='Human',Cancer=TRUE) {
  
  if (db.type == 'scType') {
    message('[[', Sys.time(), ']]: --- scTypeDB currently support Species: Human')
    message('[[', Sys.time(), ']]: --- scTypeDB currently support tissue: ')
    print(unique(ScTypeDB_full$tissueType))
    print(paste0('Select Tissue: ', Tissue))
    cellmarker <- ScTypeDB_full[grep(Tissue, ScTypeDB_full$tissueType, ignore.case = TRUE),]
    cellmarker$geneSymbolmore1 <- gsub(' ','',cellmarker$geneSymbolmore1)
    cellmarker$geneSymbolmore1 <- sapply(1:nrow(cellmarker), function(i){
      markers_all <- gsub(' ', '', unlist(strsplit(cellmarker$geneSymbolmore1[i],',')))
      markers_all <- markers_all[markers_all != 'NA' & markers_all != '']
      if (length(markers_all) > 0) {
        markers_all <- unique(na.omit(markers_all))
        paste0(markers_all, collapse = ',')
      } else {
        ''
      }
    })
    cellmarker$geneSymbolmore1 <- gsub('///',',',cellmarker$geneSymbolmore1)
    cellmarker$geneSymbolmore1 <- gsub(' ', '', cellmarker$geneSymbolmore1)
    gs <- lapply(1:nrow(cellmarker), function(j){
      gsub(' ','',unlist(strsplit(toString(cellmarker$geneSymbolmore1[j]),',')))
    })
    names(gs) <- cellmarker$cellName
    
  } else if (db.type == 'CellMarker') {
    message('[[', Sys.time(), ']]: --- CellMarker currently support Species: ')
    print(unique(Cell_marker_All$species))
    print(paste0('Select Species: ', Species))
    cellmarker <- subset(Cell_marker_All, species == Species)
    message('[[', Sys.time(), ']]: --- CellMarker currently support tissue: ')
    print(unique(cellmarker$tissue_type))
    print(paste0('Select Tissue: ', Tissue))
    cellmarker <- cellmarker[grep(Tissue,cellmarker$tissue_type, ignore.case = T),]
    message('[[', Sys.time(), ']]: --- CellMarker currently support Cells: ')
    print(unique(cellmarker$cell_type))
    if (Cancer == TRUE) {
      print(paste0('Select Cells: Cancer cell'))
      cellmarker <- subset(cellmarker, cell_type == 'Cancer cell')
    } else {
      print(paste0('Select Cells: Normal cell'))
      cellmarker <- subset(cellmarker, cell_type == 'Normal cell')
    }
    celltype <- unique(cellmarker$cell_name)
    gs <- lapply(1:length(celltype), function(i){
      markers_all <- cellmarker[cellmarker$cell_name == celltype[i],]$Symbol
      markers_all <- as.character(na.omit(markers_all))
      markers_all <- sort(markers_all)
      if (length(markers_all) > 0) {
        markers_all <- unique(na.omit(markers_all))
        paste0(markers_all, collapse = ',')
      } else {
        ''
      }
    })
    gs <- lapply(1:length(gs), function(j){
      gsub(' ', '', unlist(strsplit(toString(gs[j]),',')))
    })
    names(gs) <- celltype
    
  } else if (db.type == 'CellTaxonomy') {
    message('[[', Sys.time(), ']]: --- CellTaxonomy currently support Species: ')
    print(unique(Cell_Taxonomy$species))
    print(paste0('Select Species: ', Species))
    cellmarkers <- subset(Cell_Taxonomy, species == Species)
    message('[[', Sys.time(), ']]: --- CellTaxonomy currently support tissue: ')
    print(unique(cellmarkers$Tissue_standard))
    print(paste0('Select Tissue: ', Tissue))
    cellmarkers <- cellmarkers[grep(Tissue, cellmarkers$Tissue_standard, ignore.case = T),]
    
    celltype <- unique(cellmarkers$Cell_standard)
    gs <- lapply(1:length(celltype), function(i){
      markers_all <- cellmarkers[cellmarkers$Cell_standard == celltype[i],]$Cell_Marker
      markers_all <- as.character(na.omit(markers_all))
      markers_all <- sort(markers_all)
      if (length(markers_all) > 0) {
        markers_all <- unique(na.omit(markers_all))
        paste0(markers_all, collapse = ',')
      } else {
        ''
      }
    })
    gs <- lapply(1:length(gs), function(j){
      gsub(' ', '', unlist(strsplit(toString(gs[j]),',')))
    })
    names(gs) <- celltype
    
  } else if (db.type == 'snRandom') {
    celltype <- unique(snRandom_markers$CelltypeShortName)
    gs <- lapply(1:length(celltype), function(i){
      markers_all <- snRandom_markers[snRandom_markers$CelltypeShortName == celltype[i],]$Markers
      markers_all <- as.character(na.omit(markers_all))
      markers_all <- sort(markers_all)
      markers_all <- unique(na.omit(markers_all))
      paste0(markers_all, collapse = ',')
    })
    gs <- lapply(1:length(gs), function(j){
      gsub(' ', '', unlist(strsplit(toString(gs[j]),',')))
    })
    names(gs) <- celltype
  }
  
  return(gs)
}
