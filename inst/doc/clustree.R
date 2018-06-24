## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load-iris-----------------------------------------------------------
library(clustree)
data("iris_clusts")

head(iris_clusts)

## ----iris-plot-----------------------------------------------------------
clustree(iris_clusts, prefix = "K")

## ----iris-aes-static-----------------------------------------------------
clustree(iris_clusts, prefix = "K", node_colour = "purple", node_size = 10,
         node_alpha = 0.8)

## ----iris-aes------------------------------------------------------------
clustree(iris_clusts, prefix = "K", node_colour = "Sepal.Width",
         node_colour_aggr = "mean")

## ----iris-layout---------------------------------------------------------
clustree(iris_clusts, prefix = "K", layout = "sugiyama")

## ----iris-layout-nocore--------------------------------------------------
clustree(iris_clusts, prefix = "K", layout = "sugiyama", use_core_edges = FALSE)

## ----sim_sc3-------------------------------------------------------------
library("SingleCellExperiment")

data("sim_sc3")
sim_sc3

## ----sim_sc3-colData-----------------------------------------------------
head(colData(sim_sc3))

## ----sim_sc3-plot--------------------------------------------------------
clustree(sim_sc3, prefix = "sc3_", suffix = "_clusters")

## ----sim_seurat----------------------------------------------------------
library("Seurat")

data("sim_seurat")
sim_seurat

## ----sim_seurat-meta-----------------------------------------------------
head(sim_seurat@meta.data)

## ----sim_seurat-plot-----------------------------------------------------
clustree(sim_seurat)

## ----plot-gene-----------------------------------------------------------
clustree(sim_seurat, node_colour = "Gene730", node_colour_aggr = "median")

## ----iris-overlay--------------------------------------------------------
clustree_overlay(iris_clusts, prefix = "K", x_value = "PC1", y_value = "PC2")

## ----iris-overlay-colour-------------------------------------------------
clustree_overlay(iris_clusts, prefix = "K", x_value = "PC1", y_value = "PC2",
                 use_colour = "points", alt_colour = "blue")

## ----iris-overlay-lables-------------------------------------------------
clustree_overlay(iris_clusts, prefix = "K", x_value = "PC1", y_value = "PC2",
                 label_nodes = TRUE)

## ----iris-overlay-sides--------------------------------------------------
overlay_list <- clustree_overlay(iris_clusts, prefix = "K", x_value = "PC1",
                                 y_value = "PC2", plot_sides = TRUE)

names(overlay_list)

overlay_list$x_side
overlay_list$y_side

## ----modify--------------------------------------------------------------
clustree(iris_clusts, prefix = "K") +
    scale_color_brewer(palette = "Set1") +
    scale_edge_color_continuous(low = "blue", high = "red")

## ----legends-------------------------------------------------------------
clustree(iris_clusts, prefix = "K") +
    guides(edge_colour = FALSE, edge_alpha = FALSE) +
    theme(legend.position = "bottom")

## ----citation------------------------------------------------------------
citation("clustree")

