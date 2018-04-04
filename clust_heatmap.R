#---------------------------------------------------------
# Data Prep
#---------------------------------------------------------

# load libraries
library("gplots")  #advanced plotter - better heat maps
library("RColorBrewer")   #edit color of heat maps
library("ggplot2")   #advanced plotter for PCA
library("ggfortify")   #plotting for statistical tools (PCA)
library("cluster")    #other statistical tools (kmeans)

setwd("~/Desktop/Insight")
TF <- read.csv(file = "TF_samples.csv", header = TRUE, sep = ',', row.names = 1)
# heatmaps need input as matrix
TF_mat <- as.matrix(TF)

#---------------------------------------------------------
# Choose clustering method
#---------------------------------------------------------

### clustering standard is Pearson, use Spearman if non-parametric / ranking
## Complete (max pair-wise distance), single (min pair-wise distance) or average linkage possible
# uses cor for correlation function

# TF dendrogram
hr <- hclust(as.dist(1-cor(t(TF_mat), method="spearman")), method="complete"); 
# Sample dendrogram
hc <- hclust(as.dist(1-cor(TF_mat, method="spearman")), method="complete"); 

# Dendrogram of sample clustering
plot(hr, hang = .1)

#---------------------------------------------------------
# Heatmap
#---------------------------------------------------------

# Standard gene/TF heatmaps use red, black to green
color_scheme <- colorRampPalette(c("green", "black", "red"))(n = 100)

# Use heatmap.2 from gplots instead of R's heatmap - more customizable

heatmap.2(TF_mat,
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          # Set color
          col = color_scheme,
          # Keep scale according to provided data, use "none"; for scale to row, "row"
          scale="none",
          density.info="none",
          # change size of plot - trial and error
          margins = c(5,9),
          # Change font sizes
          cexRow = .75,
          cexCol = .35,
          # Remove trace lines
          trace="none")

#---------------------------------------------------------------
# Principal Component Analysis (via single value decomposition) and Clustering
#---------------------------------------------------------------

# ggplot use for PCA
autoplot(prcomp(TF_mat))

## Clustering using PAM (Partitioning Around Medoids), a more robust method of K-means clustering
# cluster package
pca_samples <- pam(TF_mat, 2)
autoplot(pca_samples, 
         frame = TRUE, 
         frame.type = 'norm', 
         label = TRUE, 
         label.size = 2.5,)
