##' Scale the matrix
#'
#' Will convert the angles from cell cycle genes into 10 features to compute
#' the dppc classifier.
#'
#' @param m: matrix to scale, no default
#' @return scaled matrix
#' @export
.scale = function(m) {
  x = as.matrix(m)
  
  mu = rowMeans(m)
  sigma = matrixStats::rowSds(x)
  
  scaled = sweep(m, 1, mu, "-")
  scaled = sweep(scaled, 1, sigma, "/")
  
  return(scaled)
}


##' Get the expression data from a Seurat object
#' 
#' This function will collect the gene expression data from a Seurat object.
#'
#' @param seurat_obj: the Seurat object.
#' @param assay: which assay in the Seurat object should be exported, default is 'SCT'.
#' @param layer: which layer in the Seurat object assay should be exported, default is 'data'.
#' @return Seurat expression data.
#' @export
#' 
GetSeuratExpr = function(seurat_obj, assay = 'SCT', layer = 'scale.data') {
  if (packageVersion('SeuratObject') >= '5.0.0') {
     # Seurat v5
     return(LayerData(
         object = seurat_obj,
         assay = assay,
         layer = layer
     ))
  } else {
     # Seurat v4
     return(GetAssayData(
         object = seurat_obj,
         assay = assay,
         slot = layer
     ))
  }
}


##' Predict GBM developmentally partitioned proliferative compartments
#'
#' This function predicts the cell cycle state for each cell in the object
#' using the dppc cell cycle classifier. The possible cell cycle states that
#' dppc can predict are: Neural G0, G1, G1/other, Late G1, S, S/G2, G2/M,
#' and M/Early G1.
#'
#' The dppc predicts the cell cycle state for each cell in the object by
#' selecting the cell cycle state for each cell with the maximum cell cycle
#' state probability. If the cell cycle state probability for a cell does not
#' meet the probability threshold, the cell will receive an 'Unknown' cell cycle
#' state prediction. dppc cell cycle state predictions and probabilities for
#' each cell in the object will be stored in the object .obs after classification.
#'
#' @param seurat0: a seurat object must be supplied to classify, no default
#' @param do_sctransform: whether to do SCTransform before classifying, default is TRUE
#' @param assay: which seurat_obj assay to use for classification, helpful if data is prenormalized, default is 'SCT'
#' @param gene_id: what type of gene ID is used, either 'ensembl' or 'symbol', defaults to 'ensembl'
#' @param spatial: whether the data is spatial, defaults to FALSE
#' @return Seurat object with dppc calls and probabilities for each cell cycle state
#' @export
PredictDPPC = function(seurat_obj, do_sctransform=TRUE, assay='SCT', layer = 'scale.data', gene_id='ensembl', spatial = FALSE) {
  cat('Running gbmDPPC:\n')
  # Make a copy of object
  seurat1 = seurat_obj
  
  # Load model and marker genes
  classes = read.csv(system.file('extdata', 'dppc_classes.txt', package='gbmDPPC'), header=FALSE)$V1
  marker_genes = read.csv(system.file('extdata', 'dppc_genes.csv', package='gbmDPPC'), header=TRUE, row.names=1)[,paste0('human_',gene_id)] # make file
  
  
  # Run SCTransform on data being sure to include the marker_genes
  if(assay=='SCT' & do_sctransform) {
    cat('  Redoing SCTransform to ensure maximum overlap with classifier genes...\n')
    if(!spatial) {
      seurat1 = SCTransform(seurat1, return.only.var.genes=FALSE, verbose=FALSE)
    } else {
      seurat1 = SCTransform(seurat1, assay = 'Spatial', return.only.var.genes=FALSE, verbose=FALSE)
    }
  }
  
  # Subset data marker genes to marker genes included in classification
  common_genes = intersect(row.names(seurat1),marker_genes)
  
  # Find missing genes and assign 0s to each cell
  cat(paste0('  Total possible marker genes for this classifier: ', length(marker_genes),'\n'))
  input_mat = GetSeuratExpr(seurat1, assay = assay, layer = layer)
  input_mat = input_mat[common_genes,]
  #if(assay=='SCT') {
  #    input_mat = seurat1@assays$SCT@scale.data[common_genes,]
  #} else {
  #    input_mat = seurat1@assays$RNA@data[common_genes,]
  #}
  
  missing_genes = setdiff(marker_genes, rownames(input_mat))
  
  cat(paste0('    Marker genes present in this dataset: ', nrow(input_mat),'\n'))
  cat(paste0('    Missing marker genes in this dataset: ', length(missing_genes),'\n'))
  
  
  input_mat_scaled = .scale(as.matrix(input_mat))
  
  # Create the input and output arrays (oup_preds) here with 
  # names and dimensions of marker_genes x samples (cols in seurat object 1)
  nscaled_data = matrix(min(input_mat_scaled,na.rm=T),
                        nrow = length(marker_genes),
                        ncol = ncol(seurat1), 
                        dimnames = list(marker_genes,colnames(seurat1)))
  
  # Add in the nromalized expression data from the seurat data set 
  nscaled_data[common_genes, ] = input_mat_scaled[common_genes, ]           
  nscaled_data[!is.finite(nscaled_data)] = 0
  
  
  cat(paste0('  Predicting dppc state probabilities...\n'))
  
  # Apply classifier to normalized and scaled data.  Then give the resulting output (oup) row names
  oup_preds = apply(nscaled_data, 2, dppc_classifier)
  rownames(oup_preds) = classes   
  
  # Organize the predictions into a dataframe and return the found dppc states
  # We need the dataframe with rows as samples hence the transpose here
  df1 = data.frame(t(oup_preds))
  cat(paste0('  Choosing dppc state...\n'))
    
  DPPCState = data.frame(factor(rownames(oup_preds)[apply(oup_preds,2,which.max)], levels=c('dppc1','dppc2','dppc3','dppc4','dppc5','dppc6','dppc7','dppc8')), row.names = colnames(oup_preds))
    
  colnames(DPPCState) = 'dppc'
  df1[,'dppc'] = DPPCState$dppc
  
  cat('  Adding probabilities and predictions to metadata\n')
  seurat_obj = AddMetaData(object = seurat_obj, metadata = df1)
  
  cat('Done\n')
  return(seurat_obj)
}


#' DimPlot of dppc predictions with standard colors
#'
#' This function plots the cell cycle onto a DimPlot for single cell or nuclei data.
#'
#' @param Seurat object that should have dppc cell cycle states predicted.
#' @return A DimPlot object that can be plotted.
#' @export
DimPlot.dppc = function(seurat_obj, ...) {
  dp1 = DimPlot(seurat_obj, group.by='dppc', cols = c('dppc1' = '#E41A1C', 'dppc2' = '#377EB8', 'dppc3' = '#4DAF4A','dppc4' = '#984EA3', 'dppc5' = '#FF7F00', 'dppc6' = '#E6AB02', 'dppc7' = '#A65628', 'dppc8' = '#F781BF'), ...)
  return(dp1)
}


#' SpatialDimPlot of dppc predictions with standard colors
#'
#' This function plots the cell cycle onto a DimPlot for spatial data.
#'
#' @param Seurat object with dppc cell cycle states predicted.
#' @return A DimPlot object that can be plotted.
#' @export
SpatialDimPlot.dppc = function(seurat_obj, ...) {
  dp1 = SpatialDimPlot(seurat_obj, group.by='dppc', cols = c('dppc1' = '#E41A1C', 'dppc2' = '#377EB8', 'dppc3' = '#4DAF4A','dppc4' = '#984EA3', 'dppc5' = '#FF7F00', 'dppc6' = '#E6AB02', 'dppc7' = '#A65628', 'dppc8' = '#F781BF'), ...)
  return(dp1)
}


#' dppc classifier function
#' 
#' This function calls the C code C_dppc.c to run the neural network classifier. 
#' dppc_classifier returns a 1 x 8 double precision vector of predictions 
#' 
#' @param  norm_expVec a double precision 1 x 750 vector of normalized gene expression values
#' @return  dppc_classifier returns a 1 x 8 double precision vector of class probabilites
#'
#' @export
dppc_classifier = function(norm_expVec) {
  .Call("C_dppc", norm_expVec)
}
