# phylo-rug-plot

## Directory overview

### `main.R`
R script to generate **phylogenetic rug plots** and related visualizations.

---

### `plots/`
Final figures for publication, exported in **PDF format**.

---

### `raw_trees/`
Input directory containing raw phylogenetic trees.

Two datasets are included:
- **50p**
- **70p**

For each dataset, five phylogenetic analyses were performed:

- **ASTRAL (partitioned)**
- **ASTRAL (UCE)**
- **GHOST (partitioned)**
- **IQ-TREE (partitioned)**
- **IQ-TREE (UCE)**

The tree **`70p_uce.tre`** is treated as the **primary result** and used as the **backbone topology** in downstream analyses.

---

### `trees_processed/`
Processed tree files used in figures and publications.

Trees are provided in two versions:
- **With outgroup** (316 tips)
- **Without outgroup** (289 tips)
