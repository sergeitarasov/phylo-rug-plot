library(ape)
library(phytools)
library(phangorn)
library(tibble)


transfer_node_labels <- function(target, ref) {
  # target: phylo tree whose node.label you want to (re)assign
  # ref   : phylo tree with correct node.label
  
  if (!inherits(target, "phylo") || !inherits(ref, "phylo")) {
    stop("Both 'target' and 'ref' must be objects of class 'phylo'.")
  }
  
  # Check that the trees have the same set of taxa (order can differ)
  if (!setequal(target$tip.label, ref$tip.label)) {
    stop("Trees do not have the same set of tip labels.")
  }
  
  # Helper: given a tree+node, return a key that uniquely identifies the clade
  get_clade_key <- function(tree, node) {
    tips <- Descendants(tree, node, type = "tips")[[1]]
    clade_tips <- sort(tree$tip.label[tips])
    paste(clade_tips, collapse = "|")
  }
  
  ## ---------- Build map: clade -> node label from ref ----------
  ntip_ref   <- length(ref$tip.label)
  nodes_ref  <- (ntip_ref + 1):(ntip_ref + ref$Nnode)
  # ref
  # lapply(c(1:317), function(x) get_clade_key(ref, x) %>% length) %>% unlist 
  # keys for each internal node in ref
  ref_keys <- vapply(nodes_ref, function(n) get_clade_key(ref, n), character(1))
  
  # ref$node.label is in the same order as nodes_ref
  ref_labels <- ref$node.label
  if (length(ref_labels) != length(nodes_ref)) {
    stop("Length of ref$node.label does not match ref$Nnode.")
  }
  names(ref_labels) <- ref_keys
  
  ## ---------- Assign labels to target by matching clades ----------
  ntip_target  <- length(target$tip.label)
  nodes_target <- (ntip_target + 1):(ntip_target + target$Nnode)
  
  # ensure node.label vector exists in target
  if (is.null(target$node.label)) {
    target$node.label <- rep(NA_character_, target$Nnode)
  } else if (length(target$node.label) != target$Nnode) {
    # be defensive
    target$node.label <- rep(NA_character_, target$Nnode)
  }
  
  # compute clade keys for target internal nodes
  target_keys <- vapply(nodes_target, function(n) get_clade_key(target, n), character(1))
  
  # for each target node, if its clade exists in ref, copy the label
  for (i in seq_along(nodes_target)) {
    key <- target_keys[i]
    if (!is.na(key) && key %in% names(ref_labels)) {
      target$node.label[i] <- ref_labels[[key]]
    } else {
      # leave as NA (no matching clade in ref)
      #stop()
      target$node.label[i] <- target$node.label[i]
    }
  }
  
  return(target)
}


# get_node_support(tree_list$`70p_partition_entropy.tre`, round = NULL)
get_node_support <- function(tree, round = NULL) {
  lbl <- tree$node.label
  
  # If no labels at all, return all NA
  if (is.null(lbl)) {
    return(data.frame(
      support_1 = rep(NA_real_, tree$Nnode),
      support_2 = rep(NA_real_, tree$Nnode)
    ))
  }
  
  # Ensure length == Nnode, pad with NA if needed
  if (length(lbl) < tree$Nnode) {
    lbl <- c(lbl, rep(NA_character_, tree$Nnode - length(lbl)))
  }
  
  spl <- strsplit(lbl, "/")
  
  supp_mat <- t(vapply(spl, function(x) {
    # If not of the form "val1/val2", return NA, NA
    if (length(x) < 2) {
      return(c(NA_real_, NA_real_))
    }
    as.numeric(x[1:2])
  }, numeric(2)))
  
  supp <- as.data.frame(supp_mat)
  names(supp) <- c("support_1", "support_2")
  
  if (!is.null(round)) {
    supp$support_1 <- round(supp$support_1, round)
    supp$support_2 <- round(supp$support_2, round)
  }
  
  supp
}




read_trees_from_dir <- function(dir,
                                ext = "tre",
                                format = c("newick", "nexus")) {
  
  format <- match.arg(format)
  
  # Get files with this extension
  files <- list.files(dir, pattern = paste0("\\.", ext, "$"), full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No files with extension .", ext, " found in directory ", dir)
  }
  
  # Reader function based on format
  reader <- switch(format,
                   newick = read.tree,
                   nexus  = read.nexus)
  
  # Read all trees
  trees <- lapply(files, function(f) {
    reader(f)
  })
  
  names(trees) <- basename(files)
  
  return(trees)
}

# rug_mt <- node_presence_matrix(backbone, tree_list)
# trees <- tree_list
node_presence_matrix <- function(backbone, trees) {
  # backbone: a single phylo
  # trees: a list of phylo (e.g. from read_trees_from_dir + unlist)
  
  # helper: get clade (tip labels) for a given node
  get_clade <- function(tree, node) {
    tips <- Descendants(tree, node, type = "tips")[[1]]
    sort(tree$tip.label[tips])
  }
  
  ntip <- Ntip(backbone)
  nnodes <- backbone$Nnode
  bb_nodes <- (ntip + 1):(ntip + nnodes)
  
  # clades for each backbone node
  bb_clades <- lapply(bb_nodes, get_clade, tree = backbone)
  
  # for each tree, compute clades and compare
  M <- sapply(seq_along(trees), function(i) {
    tr <- trees[[i]]
    tr_nodes <- (Ntip(tr) + 1):(Ntip(tr) + tr$Nnode)
    tr_clades <- lapply(tr_nodes, get_clade, tree = tr)
    
    # for each backbone clade: is it present in this tree?
    vapply(bb_clades, function(cl) {
      any(vapply(tr_clades, function(cl2) identical(cl2, cl), logical(1)))
    }, logical(1))
  })
  
  # rows = backbone nodes, cols = trees
  # M <- t(M)
  # colnames(M) <- paste0("node_", seq_along(bb_nodes))
  # rownames(M) <- paste0("tree_", seq_along(trees))
  colnames(M) <- names(trees)
  M <- cbind(node_id=bb_nodes, M)
  #as_tibble(M)
  M
}




node_presence_matrix2 <- function(backbone,
                                 trees,
                                 use_support   = FALSE,
                                 support_col   = 1,
                                 round_support = NULL) {
  # backbone: a single phylo
  # trees   : a named list of phylo
  # use_support: if FALSE -> 0/1 presence;
  #              if TRUE  -> support value from each tree (or 0 if absent)
  # support_col: 1 or 2 -> which column from get_node_support() to use
  
  # helper: get clade (tip labels) for a given node
  get_clade <- function(tree, node) {
    tips <- Descendants(tree, node, type = "tips")[[1]]
    sort(tree$tip.label[tips])
  }
  
  if (use_support && !support_col %in% c(1, 2)) {
    stop("support_col must be 1 or 2.")
  }
  
  ntip    <- Ntip(backbone)
  nnodes  <- backbone$Nnode
  bb_nodes <- (ntip + 1):(ntip + nnodes)
  
  # clades for each backbone node
  bb_clades <- lapply(bb_nodes, get_clade, tree = backbone)
  
  # For each tree, compute values per backbone node
  M <- sapply(seq_along(trees), function(i) {
    tr <- trees[[i]]
    
    # internal nodes of this tree
    tr_nodes  <- (Ntip(tr) + 1):(Ntip(tr) + tr$Nnode)
    tr_clades <- lapply(tr_nodes, get_clade, tree = tr)
    
    # if using support, parse support for this tree
    if (use_support) {
      supp <- get_node_support(tr, round = round_support)
      if (nrow(supp) != length(tr_nodes)) {
        stop("Number of node labels in tree ", i,
             " does not match Nnode(tree).")
      }
      support_vec <- as.numeric(supp[[support_col]])
    }
    
    # For each backbone clade: either 0/1 or that tree's support
    sapply(bb_clades, function(cl) {
      idx <- which(vapply(tr_clades,
                          function(cl2) identical(cl2, cl),
                          logical(1)))
      if (length(idx) == 0) {
        # backbone node absent in this tree
        return(0)
      } else {
        # present: either 1 or that tree's node support
        idx <- idx[1]  # unique match assumed
        if (!use_support) {
          return(1)
        } else {
          val <- support_vec[idx]
          ifelse(is.na(val), 0, val)  # NA support → treat as 0
        }
      }
    })
  })
  
  # At this point, M has rows = backbone nodes, cols = trees
  M <- as.matrix(M)
  colnames(M) <- if (!is.null(names(trees))) names(trees) else paste0("tree_", seq_along(trees))
  
  # Add backbone node IDs as first column
  M <- cbind(node_id = bb_nodes, M)
  
  M
}


# tree_list <- trees70
# verbose = TRUE
check_same_taxa <- function(tree_list, verbose = TRUE) {
  # Extract taxa from each tree
  taxa_list <- lapply(tree_list, function(tr) sort(tr$tip.label))
  
  # Use the first tree as reference
  ref <- taxa_list[[1]]
  
  # Compare all others against reference
  same <- vapply(taxa_list, function(x) identical(x, ref), logical(1))
  
  if (verbose) {
    if (all(same)) {
      message("✅ All trees contain the same set of taxa.")
    } else {
      #message("❌ Some trees differ in their taxa:.Comp a  ")
      message(paste0("❌  Some trees differ in their taxa.Compraing all trees to the Reference tree 1: ", names(ref)))
      bad <- which(!same)
      for (i in bad) {
        diff1 <- setdiff(ref, taxa_list[[i]])
        diff2 <- setdiff(taxa_list[[i]], ref)
        message(paste0(" - Tree ", i, paste0(': ', names(bad[i])), ": missing {", 
                       paste(diff1, collapse=","),
                       "} ; extra {",
                       paste(diff2, collapse=","), "}"))
      }
    }
  }
  
  return(all(same))
}





plot_node_rug <- function(tree, rug_mt,
                          cell_h, cell_w,
                          x_offset = 0.02,
                          y_offset = 0,
                          map_to_color, pal_info) {
  #   # tree:   phylo backbone used in current plot()
  #   # rug_mt: matrix / data.frame, first col = node_id, next 6 cols = support values
  #   # cell_h, cell_w: precomputed cell height/width (in data units, see above)
  #   # x_offset: fraction of max x to shift the rug to the right of the node
  #   # map_to_color: function(val, pal_info) -> color string
  #   # pal_info: whatever your map_to_color needs
  
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  dx_offset <- max(lastPP$xx) * x_offset
  dy_offset <- max(lastPP$yy) * y_offset
  
  n_nodes <- nrow(rug_mt)
  n_cells <- ncol(rug_mt) - 1
  
  #stopifnot(n_cells == 6)
  
  # n_cols <- 3
  # n_rows <- 2
  n_cols <- 2
  n_rows <- 2
  
  total_w <- n_cols * cell_w
  total_h <- n_rows * cell_h
  
  for (i in seq_len(n_nodes)) {
    node_id <- rug_mt[i, 1]
    vals <- as.numeric(rug_mt[i, -1])
    
    # node position with offsets
    x_center <- lastPP$xx[node_id] + dx_offset
    y_center <- lastPP$yy[node_id] + dy_offset
    
    # top-left corner of full 2x3 block
    x0 <- x_center - total_w / 2
    y0 <- y_center + total_h / 2
    
    for (k in seq_along(vals)) {
      val <- vals[k]
      
      row_idx <- ceiling(k / n_cols)
      col_idx <- ((k - 1) %% n_cols) + 1
      
      xleft   <- x0 + (col_idx - 1) * cell_w
      xright  <- xleft + cell_w
      ytop    <- y0 - (row_idx - 1) * cell_h
      ybottom <- ytop - cell_h
      
      # col <- if (is.na(val) || val == 0) {
      #   NA
      # } else {
      #   map_to_color(val, pal_info)
      # }
      col <-map_to_color(val, pal_info)
      
      rect(xleft, ybottom, xright, ytop,
           col = col,
           border = "black", lwd = 0.4)
    }
  }
}

rug_layout_map <- function(rug_mt, n_rows = 2, n_cols = 3) {
  # rug_mt: matrix or data.frame
  #   first column = node_id
  #   remaining columns = tree names
  #
  # returns a long data.frame:
  #   node_id | tree_name | row | col | cell_index
  
  tree_names <- colnames(rug_mt)[-1]
  n_cells <- length(tree_names)
  
  if (n_cells > n_rows * n_cols) {
    stop("More trees than available cells in the rug grid.")
  }
  
  layout_df <- do.call(rbind, lapply(seq_len(nrow(rug_mt)), function(i) {
    node_id <- rug_mt[i, 1]
    
    do.call(rbind, lapply(seq_along(tree_names), function(k) {
      data.frame(
        node_id    = node_id,
        tree_name  = tree_names[k],
        row        = ceiling(k / n_cols),           # 1..n_rows (top to bottom)
        col        = ((k - 1) %% n_cols) + 1,        # 1..n_cols (left to right)
        cell_index = k,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  layout_df
}
