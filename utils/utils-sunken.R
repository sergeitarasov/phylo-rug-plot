translate_tree_tips <- function(phy, data, from_col = "from", to_col = "to") {
  # Check inputs
  if (!inherits(phy, "phylo")) {
    stop("'phy' must be a phylogenetic tree of class 'phylo'")
  }

  if (!all(c(from_col, to_col) %in% colnames(data))) {
    stop("Specified columns not found in data")
  }

  # Create translation dictionary
  trans_dict <- setNames(data[[to_col]], data[[from_col]])

  # Find which tips need translation
  tips_to_translate <- phy$tip.label %in% names(trans_dict)

  # Perform translation
  phy$tip.label[tips_to_translate] <- trans_dict[phy$tip.label[tips_to_translate]]

  # Print summary
  cat(sprintf("Total tips: %d\n", length(phy$tip.label)))
  cat(sprintf("Tips translated: %d\n", sum(tips_to_translate)))
  cat(sprintf("Tips unchanged: %d\n", sum(!tips_to_translate)))

  return(phy)
}

#
# get_geo_regions <- function(phy, biogeo_data, geo_col = "geo", species_col = "to") {
#   # Input validation
#   if (!inherits(phy, "phylo")) {
#     stop("'phy' must be a phylogenetic tree of class 'phylo'")
#   }
#
#   if (!all(c(species_col, geo_col) %in% colnames(biogeo_data))) {
#     stop(sprintf("biogeo_data must contain '%s' and '%s' columns",
#                  species_col, geo_col))
#   }
#
#   # Filter biogeo_data to include only species in tree
#   biogeo_filtered <- biogeo_data %>%
#     filter(.data[[species_col]] %in% phy$tip.label)
#
#   # Create named vector of regions
#   regions <- biogeo_filtered[[geo_col]]
#   names(regions) <- biogeo_filtered[[species_col]]
#
#   # Match tree tips to regions
#   matched_regions <- regions[phy$tip.label]
#
#   # Print summary
#   cat(sprintf("Total tips in tree: %d\n", length(phy$tip.label)))
#   cat(sprintf("Tips matched with regions: %d\n", sum(!is.na(matched_regions))))
#   cat(sprintf("Tips without regions: %d\n", sum(is.na(matched_regions))))
#
#   return(matched_regions)
# }

get_geo_regions <- function(phy, biogeo_data, geo_col = "geo", species_col = "to") {
  # Input validation
  if (!inherits(phy, "phylo")) {
    stop("'phy' must be a phylogenetic tree of class 'phylo'")
  }

  if (!all(c(species_col, geo_col) %in% colnames(biogeo_data))) {
    stop(sprintf("biogeo_data must contain '%s' and '%s' columns",
                 species_col, geo_col))
  }

  # Create regions vector in tree tips order
  matched_regions <- character(length(phy$tip.label))
  names(matched_regions) <- phy$tip.label

  # Fill in regions for matching species
  matched_indices <- match(names(matched_regions), biogeo_data[[species_col]])
  matched_regions[!is.na(matched_indices)] <- biogeo_data[[geo_col]][matched_indices[!is.na(matched_indices)]]

  # Print summary and order check
  cat(sprintf("Total tips in tree: %d\n", length(phy$tip.label)))
  cat(sprintf("Tips matched with regions: %d\n", sum(!is.na(matched_regions))))
  cat(sprintf("Tips without regions: %d\n", sum(is.na(matched_regions))))
  cat("\nOrder check:\n")
  cat("First 5 tree tips:", paste(head(phy$tip.label, 5), collapse=", "), "\n")
  cat("First 5 region names:", paste(head(names(matched_regions), 5), collapse=", "), "\n")

  return(matched_regions)
}


convert_to_binary_matrix <- function(geo_regions) {
  # Get unique regions
  regions <- sort(unique(geo_regions[!is.na(geo_regions)]))

  # Create empty tibble with species names
  library(tidyverse)
  binary_mat <- tibble(species = names(geo_regions))

  # Add columns for each region
  for(region in regions) {
    binary_mat[[region]] <- as.numeric(geo_regions == region)
  }

  return(binary_mat)
}
