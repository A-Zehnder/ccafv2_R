##########################################################
## DPPC classifier inference for Seurat (ccAFv2-style)
##
##########################################################

.strip_ensembl_version = function(x) {
  sub('\\.[0-9]+$', '', as.character(x))
}

#' Pull an expression matrix from a Seurat v4 or v5 object.
#'
#' @param seurat_obj Seurat object
#' @param assay assay name
#' @param layer layer / slot name ("data", "scale.data", "counts")
#' @return genes x cells matrix-like object
#' @export
GetSeuratExpr = function(seurat_obj, assay = 'SCT', layer = 'data') {
  if (!assay %in% names(seurat_obj@assays)) {
    stop('Assay "', assay, '" is not in the object. Available: ',
         paste(names(seurat_obj@assays), collapse = ', '))
  }
  if (utils::packageVersion('SeuratObject') >= '5.0.0' &&
      exists('LayerData', where = asNamespace('SeuratObject'), mode = 'function')) {
    SeuratObject::LayerData(seurat_obj, assay = assay, layer = layer)
  } else if (utils::packageVersion('Seurat') >= '5.0.0') {
    Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  } else {
    Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer)
  }
}

#' Scale genes (rows) to z-scores across cells (columns).
#'
#' Same transform as Python `_scale`: (x - rowMean) / rowSd.
#'
#' @param m matrix, genes x cells
#' @return scaled matrix
#' @export
.scale = function(m) {
  x = as.matrix(m)
  mu = rowMeans(x, na.rm = TRUE)
  sigma = matrixStats::rowSds(x, na.rm = TRUE)
  sigma[!is.finite(sigma) | sigma == 0] = 1
  scaled = sweep(x, 1, mu, '-')
  scaled = sweep(scaled, 1, sigma, '/')
  scaled
}

.pkg_extdata = function(fname) {
  p = system.file('extdata', fname, package = 'ccAFv2')
  if (p == '') {
    p = system.file('extdata', fname, package = 'dppc')
  }
  p
}

.find_table = function(explicit_path, package_names) {
  if (!is.null(explicit_path) && file.exists(explicit_path)) {
    return(explicit_path)
  }
  for (nm in package_names) {
    p = .pkg_extdata(nm)
    if (nzchar(p) && file.exists(p)) return(p)
  }
  stop('Could not find required table. Looked for explicit path "',
       explicit_path, '" and package extdata files: ',
       paste(package_names, collapse = ', '),
       '. Copy dppc_marker_genes_unique.csv / dppc_classes.txt into ',
       'inst/extdata or pass genes_file / classes_file.')
}

#' Load DPPC class names in model output order.
#'
#' @param classes_file optional path to a headerless one-column file
#' @return character vector
#' @export
LoadDPPCClasses = function(classes_file = NULL) {
  path = .find_table(classes_file, c('dppc_classes.txt', 'ccAFv2_classes.txt'))
  classes = read.csv(path, header = FALSE, stringsAsFactors = FALSE)$V1
  classes = as.character(classes)
  classes = classes[!is.na(classes) & classes != '']
  if (length(classes) == 0) {
    classes = paste0('dppc', 1:8)
  }
  classes
}

#' Load DPPC marker genes in the exact order the network was trained on.
#'
#' Training did `sorted(unique(gene))` on the Ensembl column of
#' `dppc_marker_genes_unique.csv`, then saved `marker_genes_order.csv`.
#' If an order file is present it wins; otherwise genes are sorted
#' the same way.
#'
#' @param genes_file unique-marker CSV (gene, symbol)
#' @param order_file optional marker_genes_order.csv
#' @param gene_id 'ensembl' or 'symbol'
#' @return character vector of feature names in classifier order
#' @export
LoadDPPCMarkerGenes = function(genes_file = NULL, order_file = NULL,
                               gene_id = 'ensembl') {
  genes_path = .find_table(
    genes_file,
    c('dppc_marker_genes_unique.csv', 'marker_genes_order.csv',
      'ccAFv2_genes.csv')
  )
  mg = utils::read.csv(genes_path, stringsAsFactors = FALSE, check.names = FALSE)
  names(mg) = tolower(names(mg))

  if ('gene' %in% names(mg)) {
    ensembl = .strip_ensembl_version(mg$gene)
  } else {
    ensembl = .strip_ensembl_version(mg[[1]])
  }

  if (tolower(gene_id) == 'symbol') {
    if (!'symbol' %in% names(mg)) {
      stop('gene_id = "symbol" but the marker table has no "symbol" column: ',
           genes_path)
    }
    ids = as.character(mg$symbol)
    # Fall back to Ensembl when the symbol is missing so feature count
    # stays aligned with the network weights.
    bad = is.na(ids) | trimws(ids) == '' | tolower(ids) %in% c('nan', 'na', 'none', 'null', '.')
    ids[bad] = ensembl[bad]
  } else {
    ids = ensembl
  }

  order_path = NULL
  if (!is.null(order_file) && file.exists(order_file)) {
    order_path = order_file
  } else {
    cand = .pkg_extdata('marker_genes_order.csv')
    if (nzchar(cand) && file.exists(cand)) order_path = cand
  }

  if (!is.null(order_path)) {
    ord = utils::read.csv(order_path, stringsAsFactors = FALSE, check.names = FALSE)
    names(ord) = tolower(names(ord))
    col = if ('gene' %in% names(ord)) 'gene' else names(ord)[1]
    order_ids = unique(.strip_ensembl_version(ord[[col]]))
    # If the caller asked for symbols, map order Ensembl -> symbol.
    if (tolower(gene_id) == 'symbol' && 'gene' %in% names(mg) && 'symbol' %in% names(mg)) {
      map = setNames(ids, ensembl)
      mapped = unname(map[order_ids])
      mapped[is.na(mapped) | mapped == ''] = order_ids[is.na(mapped) | mapped == '']
      ids = mapped
    } else {
      ids = order_ids
    }
  } else {
    ids = sort(unique(ids[ids != '' & !is.na(ids)]))
  }

  if (length(ids) == 0) {
    stop('No marker genes loaded from ', genes_path)
  }
  ids
}

#' Map object rownames onto classifier feature IDs (version-stripped).
#'
#' @param rownames_obj character
#' @param features classifier feature IDs
#' @return named character vector: feature -> object rowname (or NA)
.match_features = function(rownames_obj, features) {
  obj_ids = .strip_ensembl_version(rownames_obj)
  map_obj = rownames_obj
  names(map_obj) = obj_ids
  # First exact, then version-stripped.
  out = setNames(rep(NA_character_, length(features)), features)
  exact = features %in% rownames_obj
  out[exact] = features[exact]
  rest = !exact
  hit = features[rest] %in% names(map_obj)
  out[rest][hit] = unname(map_obj[features[rest][hit]])
  out
}

#' Build the genes x cells classifier matrix.
#'
#' Present genes are z-scored across cells. Missing genes and non-finite
#' values are filled with min(scaled present), matching
#' 1_Final_model_loop_v2_dppc_shared_concat.py.
#'
#' @param expr genes x cells matrix
#' @param features classifier gene order
#' @return genes x cells numeric matrix
#' @export
BuildDPPCInput = function(expr, features) {
  expr = as.matrix(expr)
  rownames(expr) = as.character(rownames(expr))
  matched = .match_features(rownames(expr), features)
  present_feat = names(matched)[!is.na(matched)]
  present_rows = unname(matched[!is.na(matched)])

  if (length(present_feat) == 0) {
    stop('None of the DPPC marker genes are present in the expression matrix.')
  }

  present_mat = expr[present_rows, , drop = FALSE]
  rownames(present_mat) = present_feat
  present_scaled = .scale(present_mat)

  fill_val = suppressWarnings(min(present_scaled, na.rm = TRUE))
  if (!is.finite(fill_val)) fill_val = 0

  out = matrix(fill_val,
               nrow = length(features),
               ncol = ncol(expr),
               dimnames = list(features, colnames(expr)))
  out[present_feat, ] = present_scaled[present_feat, , drop = FALSE]
  out[!is.finite(out)] = fill_val
  out
}

#' Predict DPPC state for every cell in a Seurat object.
#'
#' Applies the compiled C network (`C_ccAFv2`) to all cells by default.
#' Training used S-phase / DPPC-labeled cells only; inference is not
#' restricted to S unless `restrict_to = "S"` and a `ccAFv2` (or similar)
#' column is already on the object.
#'
#' Predictions and per-class probabilities are written to metadata:
#'   dppc, dppc1, dppc2, ..., dppc8
#'
#' @param seurat_obj Seurat object
#' @param threshold max-probability cutoff; below this -> "Unknown"
#' @param do_sctransform rerun SCTransform when assay == "SCT"
#' @param assay assay to read
#' @param layer layer to read; use "data" to match the training export
#' @param species kept for API compatibility; DPPC markers are human GSC
#' @param gene_id "ensembl" or "symbol"
#' @param spatial if TRUE, SCTransform uses the Spatial assay
#' @param restrict_to NULL for all cells, or a named list / character
#'   filter. Character "S" keeps cells whose ccAFv2 (or cell_cycle_col)
#'   value is "S". A named list filters metadata equality, e.g.
#'   list(ccAFv2 = "S").
#' @param cell_cycle_col metadata column used when restrict_to = "S"
#' @param genes_file optional path to dppc_marker_genes_unique.csv
#' @param order_file optional path to marker_genes_order.csv
#' @param classes_file optional path to dppc_classes.txt
#' @param prediction_col metadata column for the hard call
#' @return Seurat object with DPPC probabilities and calls
#' @export
PredictDPPC = function(seurat_obj,
                       threshold = 0,
                       do_sctransform = FALSE,
                       assay = 'SCT',
                       layer = 'data',
                       species = 'human',
                       gene_id = 'ensembl',
                       spatial = FALSE,
                       restrict_to = NULL,
                       cell_cycle_col = 'ccAFv2',
                       genes_file = NULL,
                       order_file = NULL,
                       classes_file = NULL,
                       prediction_col = 'dppc') {
  cat('Running DPPC classifier:\n')
  seurat1 = seurat_obj

  classes = LoadDPPCClasses(classes_file)
  marker_genes = LoadDPPCMarkerGenes(genes_file, order_file, gene_id)
  n_features = length(marker_genes)
  cat(paste0('  Classes (', length(classes), '): ',
             paste(classes, collapse = ', '), '\n'))
  cat(paste0('  Classifier features: ', n_features, '\n'))

  if (assay == 'SCT' && isTRUE(do_sctransform)) {
    cat('  Redoing SCTransform to maximize overlap with classifier genes...\n')
    if (!spatial) {
      seurat1 = Seurat::SCTransform(seurat1, return.only.var.genes = FALSE,
                                    verbose = FALSE)
    } else {
      seurat1 = Seurat::SCTransform(seurat1, assay = 'Spatial',
                                    return.only.var.genes = FALSE,
                                    verbose = FALSE)
    }
  }

  if (!assay %in% names(seurat1@assays)) {
    assay = Seurat::DefaultAssay(seurat1)
    cat(paste0('  Requested assay missing; using DefaultAssay: ', assay, '\n'))
  }

  input_mat = GetSeuratExpr(seurat1, assay = assay, layer = layer)
  rownames(input_mat) = as.character(rownames(input_mat))

  cells_use = colnames(input_mat)
  if (!is.null(restrict_to)) {
    md = seurat1[[]]
    if (is.character(restrict_to) && length(restrict_to) == 1 &&
        is.null(names(restrict_to)) && restrict_to %in% c('S', 's')) {
      if (!cell_cycle_col %in% colnames(md)) {
        stop('restrict_to = "S" needs metadata column ', cell_cycle_col)
      }
      keep = which(as.character(md[[cell_cycle_col]]) == 'S')
    } else if (is.list(restrict_to)) {
      keep = rep(TRUE, nrow(md))
      for (nm in names(restrict_to)) {
        keep = keep & as.character(md[[nm]]) %in% as.character(restrict_to[[nm]])
      }
      keep = which(keep)
    } else if (is.character(restrict_to) && cell_cycle_col %in% colnames(md)) {
      keep = which(as.character(md[[cell_cycle_col]]) %in% restrict_to)
    } else {
      stop('restrict_to must be NULL, "S", a character vector of ',
           cell_cycle_col, ' values, or a named list of metadata filters.')
    }
    cells_use = rownames(md)[keep]
    cells_use = intersect(cells_use, colnames(input_mat))
    cat(paste0('  Restricted to ', length(cells_use), ' / ',
               ncol(input_mat), ' cells\n'))
    if (length(cells_use) == 0) {
      warning('No cells left after restrict_to; returning object unchanged.')
      return(seurat_obj)
    }
    input_mat = input_mat[, cells_use, drop = FALSE]
  } else {
    cat(paste0('  Scoring all ', length(cells_use), ' cells (all phases)\n'))
  }

  matched = .match_features(rownames(input_mat), marker_genes)
  n_present = sum(!is.na(matched))
  n_missing = n_features - n_present
  cat(paste0('    Marker genes present: ', n_present, '\n'))
  cat(paste0('    Marker genes missing: ', n_missing, '\n'))
  if (n_present / n_features < 0.8) {
    warning('Overlap below 80% (', n_present, '/', n_features,
            '). Try do_sctransform = TRUE or check gene_id.')
  }

  nscaled_data = BuildDPPCInput(input_mat, marker_genes)
  if (nrow(nscaled_data) != n_features) {
    stop('Built input has ', nrow(nscaled_data),
         ' rows but the network expects ', n_features,
         '. Check marker gene order files.')
  }

  cat('  Predicting DPPC state probabilities...\n')
  oup_preds = apply(nscaled_data, 2, ccAFv2_classifier)
  oup_preds = as.matrix(oup_preds)
  if (nrow(oup_preds) != length(classes)) {
    if (nrow(oup_preds) == 1 && length(oup_preds) == length(classes) * ncol(nscaled_data)) {
      oup_preds = matrix(oup_preds, nrow = length(classes))
    } else if (ncol(oup_preds) == length(classes) && nrow(oup_preds) == ncol(nscaled_data)) {
      oup_preds = t(oup_preds)
    } else {
      stop('C classifier returned ', paste(dim(oup_preds), collapse = ' x '),
           ' but ', length(classes), ' class rows were expected. ',
           'Regenerate C_ccAFv2 so the wrapper output length is ',
           length(classes), ' (dppc1-dppc8).')
    }
  }
  rownames(oup_preds) = classes

  cat('  Choosing DPPC state...\n')
  max_state = rownames(oup_preds)[apply(oup_preds, 2, which.max)]
  df1 = data.frame(t(oup_preds), check.names = FALSE)
  df1[[prediction_col]] = factor(max_state, levels = c(classes, 'Unknown'))
  df1[apply(oup_preds, 2, max) < threshold, prediction_col] = 'Unknown'
  rownames(df1) = colnames(oup_preds)

  # Make.names-safe copies so AdjustDPPCThreshold can find columns later
  # even if metadata names get sanitized.
  for (cl in classes) {
    df1[[make.names(cl)]] = df1[[cl]]
  }

  cat('  Adding probabilities and predictions to metadata\n')
  # Cells not scored stay NA for the new columns.
  seurat_obj = Seurat::AddMetaData(object = seurat_obj, metadata = df1)
  cat('Done\n')
  seurat_obj
}

#' Drop-in name used by the U5 package, now routing to PredictDPPC.
#'
#' Kept so existing scripts that call PredictCellCycle still run after
#' you swap in the DPPC C library. U5-only arguments (`include_g0`) are
#' accepted and ignored.
#'
#' @inheritParams PredictDPPC
#' @export
PredictCellCycle = function(seurat_obj, threshold = 0, include_g0 = FALSE,
                            do_sctransform = FALSE, assay = 'SCT',
                            layer = 'data', species = 'human',
                            gene_id = 'ensembl', spatial = FALSE, ...) {
  if (!missing(include_g0)) {
    message('PredictCellCycle is the DPPC wrapper; include_g0 is ignored.')
  }
  PredictDPPC(seurat_obj,
             threshold = threshold,
             do_sctransform = do_sctransform,
             assay = assay,
             layer = layer,
             species = species,
             gene_id = gene_id,
             spatial = spatial,
             ...)
}

#' Re-apply a probability threshold to stored DPPC calls.
#'
#' @param seurat_obj Seurat object already run through PredictDPPC
#' @param threshold new cutoff
#' @param classes_file optional class file
#' @param prediction_col metadata column to overwrite
#' @return Seurat object
#' @export
AdjustDPPCThreshold = function(seurat_obj, threshold = 0.5,
                               classes_file = NULL,
                               prediction_col = 'dppc') {
  cat('Adjusting DPPC threshold:\n')
  classes = LoadDPPCClasses(classes_file)
  md = seurat_obj[[]]
  pred_cols = intersect(c(classes, make.names(classes)), colnames(md))
  pred_cols = pred_cols[!duplicated(make.names(pred_cols))]
  if (length(pred_cols) < 2) {
    stop('No DPPC probability columns found in metadata. Run PredictDPPC first.')
  }
  predictions1 = md[, pred_cols, drop = FALSE]
  # Restore original class names if make.names changed them.
  map = setNames(classes, make.names(classes))
  colnames(predictions1) = ifelse(colnames(predictions1) %in% names(map),
                                  unname(map[colnames(predictions1)]),
                                  colnames(predictions1))
  max_state = colnames(predictions1)[apply(predictions1, 1, which.max)]
  call = as.character(max_state)
  call[apply(predictions1, 1, max) < threshold] = 'Unknown'
  seurat_obj[[prediction_col]] = factor(call, levels = c(classes, 'Unknown'))
  cat('Done\n')
  seurat_obj
}

#' @rdname AdjustDPPCThreshold
#' @export
AdjustCellCycleThreshold = function(seurat_obj, threshold = 0.5,
                                    include_g0 = FALSE, ...) {
  AdjustDPPCThreshold(seurat_obj, threshold = threshold, ...)
}

#' Add per-DPPC module scores for later regression.
#'
#' Reads `dppc_filtered_markers_by_dppc.csv` (dppc, gene, symbol) produced
#' by GBM_dppc_marker_filter_csv_dppc2_union.R. Falls back to the unique
#' marker table (no per-state split) if the by-dppc file is absent.
#'
#' @param seurat_obj Seurat object
#' @param assay assay for AddModuleScore
#' @param gene_id "ensembl" or "symbol"
#' @param markers_by_dppc_file optional path
#' @return Seurat object with module-score columns
#' @export
PrepareForDPPCRegression = function(seurat_obj, assay = 'SCT',
                                    gene_id = 'ensembl',
                                    markers_by_dppc_file = NULL) {
  path = tryCatch(
    .find_table(markers_by_dppc_file,
                c('dppc_filtered_markers_by_dppc.csv',
                  'dppc_marker_genes_unique.csv')),
    error = function(e) NULL
  )
  if (is.null(path)) {
    stop('Need dppc_filtered_markers_by_dppc.csv to build per-state gene lists.')
  }
  mg = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(mg) = tolower(names(mg))
  id_col = if (tolower(gene_id) == 'symbol' && 'symbol' %in% names(mg)) 'symbol' else 'gene'
  if (!'dppc' %in% names(mg)) {
    stop(path, ' has no "dppc" column; cannot split genes by state.')
  }
  cluster_genes = lapply(split(mg[[id_col]], mg$dppc), function(x) {
    unique(as.character(x[!is.na(x) & x != '']))
  })
  cluster_genes = cluster_genes[lengths(cluster_genes) > 0]
  Seurat::AddModuleScore(seurat_obj, features = cluster_genes, assay = assay,
                         name = paste0(names(cluster_genes), '_exprs'))
}

dppc_colors = c(
  dppc1 = '#a6cee3',
  dppc2 = '#1f78b4',
  dppc3 = '#b2df8a',
  dppc4 = '#33a02c',
  dppc5 = '#fb9a99',
  dppc6 = '#e31a1c',
  dppc7 = '#fdbf6f',
  dppc8 = '#ff7f00',
  Unknown = '#cccccc'
)

#' DimPlot of DPPC predictions with the training-script palette.
#'
#' @param seurat_obj Seurat object with a dppc column
#' @param group.by metadata column
#' @param ... passed to DimPlot
#' @export
DimPlot.dppc = function(seurat_obj, group.by = 'dppc', ...) {
  Seurat::DimPlot(seurat_obj, group.by = group.by, cols = dppc_colors, ...)
}

#' @rdname DimPlot.dppc
#' @export
DimPlot.ccAFv2 = function(seurat_obj, group.by = 'dppc', ...) {
  DimPlot.dppc(seurat_obj, group.by = group.by, ...)
}

#' SpatialDimPlot of DPPC predictions.
#'
#' @param seurat_obj Seurat object
#' @param group.by metadata column
#' @param ... passed to SpatialDimPlot
#' @export
SpatialDimPlot.dppc = function(seurat_obj, group.by = 'dppc', ...) {
  Seurat::SpatialDimPlot(seurat_obj, group.by = group.by, cols = dppc_colors, ...)
}

#' @rdname SpatialDimPlot.dppc
#' @export
SpatialDimPlot.ccAFv2 = function(seurat_obj, group.by = 'dppc', ...) {
  SpatialDimPlot.dppc(seurat_obj, group.by = group.by, ...)
}

#' Stacked barplot of DPPC calls across thresholds.
#'
#' @param seurat_obj Seurat object with DPPC probability columns
#' @param classes_file optional class file
#' @param ... unused, kept for compatibility
#' @export
ThresholdPlot = function(seurat_obj, classes_file = NULL, ...) {
  classes = LoadDPPCClasses(classes_file)
  md = seurat_obj[[]]
  pred_cols = intersect(c(classes, make.names(classes)), colnames(md))
  pred_cols = pred_cols[!duplicated(make.names(pred_cols))]
  predictions1 = md[, pred_cols, drop = FALSE]
  map = setNames(classes, make.names(classes))
  colnames(predictions1) = ifelse(colnames(predictions1) %in% names(map),
                                  unname(map[colnames(predictions1)]),
                                  colnames(predictions1))

  make_calls = function(th) {
    call = colnames(predictions1)[apply(predictions1, 1, which.max)]
    call[apply(predictions1, 1, max) < th] = 'Unknown'
    factor(call, levels = c(classes, 'Unknown'))
  }

  dfall = data.frame(table(make_calls(0)) / nrow(predictions1))
  names(dfall)[1:2] = c('dppc', 'Freq')
  dfall$Threshold = '0'
  for (threshold in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) {
    df1 = data.frame(table(make_calls(threshold)) / nrow(predictions1))
    names(df1)[1:2] = c('dppc', 'Freq')
    df1$Threshold = as.character(threshold)
    dfall = rbind(dfall, df1)
  }
  ggplot2::ggplot(dfall) +
    ggplot2::geom_bar(ggplot2::aes(x = Threshold, y = Freq, fill = dppc),
                      position = 'stack', stat = 'identity') +
    ggplot2::scale_fill_manual(values = dppc_colors) +
    ggplot2::theme_minimal()
}



#'
#' @param norm_expVec numeric vector, length = 750
#' @return numeric vector of class probabilities, length = 8
#' @export
ccAFv2_classifier = function(norm_expVec) {
  .Call('C_ccAFv2', as.numeric(norm_expVec))
}
