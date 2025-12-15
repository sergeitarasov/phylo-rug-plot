# - There five analysis for ech 50p, and 70p (the backbone)
# - astral - partition
# - astral  - UCE
# - ghost -partition
# - iqtree - partition
# - iqtree - UCE
#
# iqtree_70p_UCEs.tre is a BACKBONE

library(readxl)
source('utils/utils.R')
source('utils/utils-sunken.R')

# tr <- read.nexus("trees/70p/iqtree_70p_UCEs.tre")
# read.nexus("trees/70p/astral_70p_partitions.tre")



#-------------------- 70p -------------------------------------------------------

# this file is used fro translations
biogeo_data <- read_excel('data/biogeo.xlsx', sheet = "BioGeo")

trees70 <- read_trees_from_dir(dir="raw_trees/70p", ext = "tre", format = "newick")
trees70 <- trees70
names(trees70)
check_same_taxa(trees70)


file <- file.path("plots", "70p-all-trees.pdf")
pdf(file = file, width = 8.27*2, height = 11.69 *2)
# i=1
for (i in 1:length(trees70)){
  file.name <- names(trees70[i])
  tr <- trees70[[i]]
  tr <- root(tr, outgroup=c('NicorbUCE', 'NicvesUCE'), resolve.root = T)
  tr <- ladderize(tr)
  #tr <- transfer_node_labels(target=tr, ref=trees70[[i]]) 
  #tr <- reorder(tr, order = "cladewise")
  # plot.phylo(tr, show.tip.label = F)
  #
  #Translate tree tips
  tr <- translate_tree_tips(phy = tr,
                            data = biogeo_data,
                            from_col = "from",  # your actual column name
                            to_col = "to")    # your actual column name
  
  file <- file.path("trees_processed", "316-tips", file.name)
  write.tree(tr, file)
  #
  # #------------ Extract Scarabaeinae + Aphodiinae
  node <- findMRCA(tr, c('Heteronitis_ragazzii_STL10032', 'Aphodius_immundus_ST003'))
  tr <- extract.clade(tr, node)
  #plot.phylo(tr, show.tip.label = F)
  file <- file.path("trees_processed", "289-tips", file.name)
  write.tree(tr, file)
  trees70[[i]] <- tr
  # plot to pdf
  plot.phylo(tr, show.tip.label = T, align.tip.label = F, 
             cex=0.4, label.offset=.005, no.margin = F, 
             main=file.name)
  nodelabels(tr$node.label, frame = 'none', cex=0.3, col='red')
}
dev.off()

#-------------------- -------------------------------------------------------
#--------------------  50p -------------------------------------------------------

library(readxl)
source('utils/utils.R')
source('utils/utils-sunken.R')

biogeo_data <- read_excel('data/biogeo.xlsx', sheet = "BioGeo")
trees50 <- read_trees_from_dir(dir="raw_trees/50p", ext = "tre", format = "newick")
names(trees50)
check_same_taxa(trees50)


file <- file.path("plots", "50p-all-trees.pdf")
pdf(file = file, width = 8.27*2, height = 11.69 *2)
# i=3
for (i in 1:length(trees50)){
  file.name <- names(trees50[i])
  tr <- trees50[[i]]
  tr <- root(tr, outgroup=c('NicorbUCE', 'NicvesUCE'), resolve.root = T)
  tr <- ladderize(tr)
  #Translate tree tips
  tr <- translate_tree_tips(phy = tr,
                            data = biogeo_data,
                            from_col = "from",  # your actual column name
                            to_col = "to")    # your actual column name
  
  file <- file.path("trees_processed", "316-tips", file.name)
  write.tree(tr, file)
  #
  # #------------ Extract Scarabaeinae + Aphodiinae
  node <- findMRCA(tr, c('Kurtops_signatus _STL10043', 'Aphodius_immundus_ST003'))
  #node <- findMRCA(tr, c('Heteronitis_ragazzii_STL10032', 'Aphodius_immundus_ST003'))
  tr <- extract.clade(tr, node)
  file <- file.path("trees_processed", "289-tips", file.name)
  write.tree(tr, file)
  trees50[[i]] <- tr
  # plot to pdf
  plot.phylo(tr, show.tip.label = T, align.tip.label = F, 
             cex=0.4, label.offset=.005, no.margin = F, 
             main=file.name)
  nodelabels(tr$node.label, frame = 'none', cex=0.3, col='red')
}
dev.off()

#-------------------- -------------------------------------------------------
#--------------------  get matrix  -------------------------------------------------------
# get_node_support(tree_list$`70p_partition_entropy.tre`, round = NULL)
#rug_mt <- node_presence_matrix(backbone, tree_list)

library(readxl)
source('utils/utils.R')
source('utils/utils-sunken.R')

trees <- read_trees_from_dir(dir="trees_processed/289-tips", ext = "tre", format = "newick")
names(trees)
check_same_taxa(trees)

# making fake node labels as two columns for ASTRAL, so it can be parsed
lb <- trees$`70p_ASTRAL_uce.tre`$node.label
trees$`70p_ASTRAL_uce.tre`$node.label <- paste0(lb, '/', lb)

lb <- trees$`70p_ASTRAL_partition_entropy.tre`$node.label
trees$`70p_ASTRAL_partition_entropy.tre`$node.label <- paste0(lb, '/', lb)
#------

backbone <- trees$`70p_uce.tre`

tree_list <- trees[c("70p_partition_entropy.tre", "70p_ASTRAL_partition_entropy.tre", "70p_ASTRAL_uce.tre",              
  "70p_ghost.tre",  "50p_uce.tre", "50p_partition_entropy.tre")]                          

names(tree_list)

#---- RF distance


#trees_unr <- lapply(trees, function(x) unroot(x))
class(trees) <- 'multiPhylo'
rf <- dist.topo(trees, method = "PH85")
rf
rf.mt <- as.matrix(rf)
sort(rf)
#---


rug_mt <- node_presence_matrix2(backbone,
                                tree_list,
                                use_support   = T,
                                support_col   = 1,    # 1 or 2
                                round_support = 2)

rug_mt
rug_mt[,'70p_ASTRAL_uce.tre'] <- rug_mt[,'70p_ASTRAL_uce.tre']*100
rug_mt[,'70p_ASTRAL_partition_entropy.tre'] <- rug_mt[,'70p_ASTRAL_partition_entropy.tre']*100
rug_mt

plot.phylo(trees$`50p_uce.tre`, show.tip.label = F)
nodelabels()

#-------------------- -------------------------------------------------------
names(trees)

### Make colors

# pal_info <- list()
# pal_info$pal <- gplots::rich.colors(64, palette="temperature", alpha=1.0, rgb=FALSE, plot=F)
# pal_info$pal[1:8] <- pal_info$pal[9]
# breaks<- BAMMtools::assignColorBreaks(seq(0,100, 1), spex=NULL, NCOLORS = 64, method="jenks")
# pal_info$breaks <- unique(breaks)

pal_info <- list()

# Number of colors in gradient
ncol <- 64

# Light grey to black
pal_info$pal <- gray.colors(
  ncol,
  start = 0.95,  # light grey
  end   = 0.00   # black
)

# Breaks from 0 to 1
pal_info$breaks <- seq(0, 100, length.out=ncol)

# Fun: vals -> colors
map_to_color <- function(vals, pal_info) {
  cols <- pal_info$pal[cut(vals, breaks = pal_info$breaks, include.lowest = TRUE)]
  #cols[vals == log(eps)] <- "#ffffff"  # force exact zeros to white
  #cols[vals == 0] <- "#ffffff"
  #cols[vals == 0] <- "lightgrey"
  cols[vals == 0] <- 'white'  #"aliceblue"  #"skyblue1" "#deebf7"  # this is a special color for 0s
  cols
}
map_to_color(70, pal_info)


#--------------------  RUG PLOT  -------------------------------------------------------
layout_df <- rug_layout_map(rug_mt)
head(layout_df)

# filter out 100 support
rug_mt_filt <- rug_mt[
  !apply(rug_mt[, -1, drop = FALSE], 1, function(x) all(x == 100)),
  ,
  drop = FALSE
]

rug_mt_100 <- rug_mt[
  apply(rug_mt[, -1, drop = FALSE], 1, function(x) all(x == 100)),
  ,
  drop = FALSE
]



file_base <- "plots/rug_plot_70p_uce"
# A4 fromat: width = 8.27, height = 11.69
pdf(
  file = paste0(file_base, ".pdf"),
  width = 8.27*2, height = 11.69 *2
)

# important patamters:
# cell_h <- dy * 1: cells width
# x_offset = -.014,       # tweak horizontal distance from node
# y_offset = 0.005,

plot.phylo(backbone, show.tip.label = T, align.tip.label = F, cex=0.4, label.offset=.005, no.margin = T)
#nodelabels(node=rug_mt_100[,1])
 
Ntip <- Ntip(backbone)

# nodelabels(text='X', node=577)
# rug_mt[,1]

lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# 100% nodes as points
# Vector of node IDs to mark
node_ids <- rug_mt_100[,1]  # example
# Plot small dots
points(
  lastPP$xx[node_ids],
  lastPP$yy[node_ids],
  pch = 16,        # solid circle
  cex = 0.8,       # size
  col = "black"    # or any color
)

# Compute square cell sizes (as above)
pin <- par("pin")
xrange <- diff(lastPP$x.lim)
yrange <- diff(lastPP$y.lim)
x_per_inch <- xrange / pin[1]
y_per_inch <- yrange / pin[2]
dy <- median(diff(sort(lastPP$yy[1:Ntip])))

cell_h <- dy * 1
cell_w <- cell_h * (x_per_inch / y_per_inch)

# Draw node rugs
plot_node_rug(tree, rug_mt_filt,
              cell_h = cell_h,
              cell_w = cell_w,
              x_offset = -.014,       # tweak horizontal distance from node
              y_offset = 0.005,
              map_to_color = map_to_color,
              pal_info = pal_info)

dev.off()

# # ==========================
# # Navajo rug on NODES (2x3 grid per node)
# # ==========================
# 
# lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
# 
# # Node coordinates (ape's internal numbering)
# node_x <- lastPP$xx
# node_y <- lastPP$yy
# 
# # Rug data:
# # first column: node_id (ape node index)
# # remaining columns: support for each analysis/tree
# node_ids   <- rug_mt[, 1]
# support_mt <- as.matrix(rug_mt[, -1, drop = FALSE])
# n_nodes    <- nrow(support_mt)
# n_cells    <- ncol(support_mt)  # number of trees / cells
# 
# # We want a 2 (rows) x 3 (cols) layout
# n_cols <- 3
# n_rows <- 2
# if (n_cells > n_cols * n_rows) {
#   stop("rug_mt has more columns than 2x3 cells can display.")
# }
# 
# # Size of each cell (adjust these multipliers to taste)
# box_w <- max(lastPP$xx) * 0.012
# box_h <- max(lastPP$yy) * 0.012
# # box_w <-.1
# # box_h <- .1
# some_offset <- max(lastPP$xx) * 0.01
# 
# 
# i=1
# for (i in seq_len(n_nodes)) {
#   nid <- node_ids[i]
#   x_center <- node_x[nid]
#   y_center <- node_y[nid]
#   
#   # total rug size
#   total_w <- n_cols * box_w
#   total_h <- n_rows * box_h
#   
#   # top-left corner of the whole 2x3 rug
#   # x_start <- x_center - total_w / 2
#   # y_start <- y_center + total_h / 2
#   x_start <- x_center - 0.01
#   y_start <- y_center + 5
#   
#   # draw each cell
#   # k=1
#   for (k in seq_len(n_cells)) {
#     val <- support_mt[i, k]
#     
#     # # skip if NA or 0 support
#     # if (is.na(val) || val == 0) {
#     #   fill_col <- NA
#     # } else {
#     #   fill_col <- map_to_color(val, pal_info)
#     # }
#     fill_col <- map_to_color(val, pal_info)
#     
#     # map cell index k -> grid (row, col)
#     # k = 1..6; fill row-wise from top-left
#     cell_row <- (k - 1) %/% n_cols  # 0 (top) or 1 (bottom)
#     cell_col <- (k - 1) %%  n_cols  # 0,1,2
#     
#     xleft   <- x_start + cell_col * box_w
#     xright  <- xleft + box_w
#     ytop    <- y_start - cell_row * box_h
#     ybottom <- ytop - box_h
#     
#     rect(xleft, ybottom, xright, ytop,
#          col    = fill_col,
#          border = "black", lwd = 0.6)
#   }
# }
