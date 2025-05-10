# R_Script_Hierarchical_Gene_Cluster_Analysis
# Author_Dr_Asad_Prodhan_prodhan82@gmail.com
# 
## Load required libraries
## DEG
library(DESeq2)   # DEG analysis
library(ggplot2)  # Making plots
library(genefilter) # Allows 'DESeq2 result' to filter genes based on adjusted-p   
library(apeglm)   # Allows DESeq2 to produce shrunken log2 fold changes to reduce the noise
library(ashr)     # Alternative to ashr
library(ggrepel)  # EnhancedVolcano requires this 
library(EnhancedVolcano) # allows making Volcano plot
library(magrittr) # provides with a set of operators like pipe %>%
library(dplyr)    # Data organisation and sub-setting # needed for cluster plots
library(tidyr)    # Data manipulation # needed for cluster plots
library(stringr)  # String match and filter
library(tibble)   # Making table
## Venn diagram
library(grid)           # VennDiagram requires this
library(futile.logger)  # VennDiagram requires this
library(VennDiagram)    # Venn diagram
## Combining plots
library(ggpubr)         # Combining plots
library(gridExtra)      # Combining plots
library(cowplot)        # Combining plots  
## PCA
library(usethis)        # Required by devtools
library(devtools)
library(plyr)
library(scales)
library(ggbiplot)       # Requires ggplot2, plyr, scales, grid
library(ggforce)        # Drawing ellipse with 3 samples
                        # Stat_ellipe (which doesn't require ggforce) for greater samples
## Heatmap
#library(pheatmap)      # pheatmap conflicts with ComplexHeatmap
## Colouring the row dendogram
library(ComplexHeatmap)
library(grid)           # For gpar function
library(dendextend)     # For dendrogram manipulation
library(wesanderson)    # For heatmap color
library(circlize)       # For heatmap color
#
## Read the data 
#
# Load gene(/transcript) count matrix and labels
countData <- as.matrix(read.csv("gene_count_matrix_hakea_young_roots_v2.csv", row.names="Gene_name"))
countData
metaData <- read.csv("metadata_hakea_young_roots.csv", row.names = 1)
metaData
#
# Check all sample IDs in metaData are also in CountData and match their orders
all(rownames(metaData) %in% colnames(countData))
countData <- countData[, rownames(metaData)]
countData
all(rownames(metaData) == colnames(countData))
#
# Check the countData for non-integers and missing values, otherwise there will be the following error
# Error in DESeqDataSet(se, design = design, ignoreRank) : some values in assay are not integers
#
str(countData)
summary(countData)
#
# Convert Non-Integers to Integers
# If the counts are floating-point numbers (e.g., decimals), you can round them:
countData <- round(countData)
countData
#
# Check for NA or Non-Numeric Values
any(is.na(countData))
any(!is.numeric(countData))
any(countData == 0)
sum(countData == 0)
#
# Verify Input Data Structure
# Confirm that countData is a matrix or a data.frame of numeric (integer) values
is.matrix(countData) 
is.data.frame(countData)
#
# Make treatments as a factor
metaData$Treatments <- as.factor(metaData$Treatments)
# check if the Treatments column is a character or factor using
str(metaData)  # Check data structure
table(metaData$Treatments)  # See levels of the factor
#
## Create a DESeqDataSet from count matrix and labels
dds_yr <- DESeqDataSetFromMatrix(countData = countData, colData = metaData,
                                 design = ~Treatments)
dds_yr
nrow(dds_yr)
#
# For Hierarchial Gene Cluster analysis, DESeq2(dds_yr) is not required
# DESeq2(dds_yr) for analysing and determining the DEGs
#
## Transform the dataset using vst or rlog. Transformation required for stabilizing variance
# variance-stabilizing transformation (vst) is similar to r-log transformation but suitable for large dataset
#
vsd <- vst(dds_yr, blind = FALSE)
vsd
colData(vsd)
#
## Extract the transformed data from DESeqDataSet using assay()
mat  <- assay(vsd)
head(mat)
#
#
any(mat == 0)
sum(mat == 0)
## Update the sample names
# Old colnames
colnames(mat) 
# New colnames
name_mapping <- c(
  "YR10" = "N0_R1",
  "YR10_2" = "N0_R2",
  "YR10_3" = "N0_R3",
  "YR_mix_new" = "N5_R1",
  "YR18" = "N5_R2",
  "YR19" = "N5_R3",
  "YR9" = "N25_R1",
  "YR15" = "N25_R2",
  "YR21" = "N25_R3"
)
# Replace column names in mat
colnames(mat) <- name_mapping[colnames(mat)]
colnames(mat)
#
# Calculate the z-scores (z-score normalization)
z_scores <- t(scale(t(mat)))
z_scores
head(z_scores)
#
any(z_scores == 0)
sum(z_scores == 0)
# Check for missing numbers
any(is.na(z_scores))
any(is.nan(z_scores))
any(is.infinite(z_scores))
sum(is.na(z_scores))
# If any TRUE is returned, inspect and clean those entries accordingly
#
# Save the z_scores
write.csv(z_scores, file = "z_scores.csv", row.names = FALSE)
#
# Calculate the replicate averages
#
average_by_treatment <- function(df, treatment_groups) {
  df <- as.data.frame(df)
  averaged_df <- data.frame(row.names = rownames(df))  # Create a new dataframe with row names
  for (treatment in names(treatment_groups)) {
    averaged_df[[treatment]] <- rowMeans(df[, treatment_groups[[treatment]], drop = FALSE], na.rm = TRUE)
  }
  return(averaged_df)
}
# Define treatment groups
treatment_groups <- list(
  N0 = c("N0_R1", "N0_R2", "N0_R3"),
  N5 = c("N5_R1", "N5_R2", "N5_R3"),
  N25 = c("N25_R1", "N25_R2", "N25_R3")
)
# Compute averages and only keep the average columns
z_scores_avg <- average_by_treatment(z_scores, treatment_groups)
z_scores_avg
head(z_scores_avg)
#
# Check if the z_scores_avg is a matrix or dataframe 
is.matrix(z_scores_avg)
is.data.frame(z_scores_avg)
# if not a dataframe, then convert as follows:
#z_scores_avg <- as.data.frame(z_scores_avg)
#
# Does z_scores_avg has any missing value?
sum(is.na(z_scores_avg)) 
# Replace all NA (missing) values in the z_scores_avg object with 0
z_scores_avg[is.na(z_scores_avg)] <- 0
sum(is.na(z_scores_avg))  # Count NAs
#
# Now, cluster the row
row_dist <- dist(z_scores_avg)  # Compute distance matrix
head(row_dist)
row_clustering <- hclust(row_dist)  # Hierarchical clustering
plot(row_clustering, main = "Hierarchical Clustering of Genes", xlab = "Genes")
#
## Cut the dendrogram at a height of your interest
# Set a threshold
height_threshold <- 2  # Adjust this value as needed
# Add a horizontal line at the height threshold
abline(h = height_threshold, col = "red")
# Cut the tree
row_clusters_by_height <- cutree(row_clustering, h = height_threshold)
row_clusters_by_height
# Get the number of clusters based on the height threshold
num_clusters_by_height <- length(unique(row_clusters_by_height))
num_clusters_by_height
# View the number of genes in each cluster
table(row_clusters_by_height)  # This will show the number of genes in each cluster
# If you sum up all the clusters, then the total sum will be the total number of rows
#
# Convert the row clusters to dendrogram
row_dendrogram <- as.dendrogram(row_clustering)
row_dendrogram
#
#
## Make a HeatMap
#
# Define the sample order on the x-axis
colnames(z_scores_avg)
sample_order <- c("N0", "N5", "N25")
sample_order
# Match the column names of z_scores object to the above sample order
column_order <- match(sample_order, colnames(z_scores_avg))
column_order
#
ncol(z_scores_avg)
length(column_order)
print(colnames(z_scores_avg))  # Check column names
print(column_order)  # Check current order
#
# Get the order of samples from the dendrogram
dend_order <- order.dendrogram(row_dendrogram)
dend_order
# Reorder the cluster assignment to match dendrogram order
ordered_clusters <- row_clusters_by_height[dend_order]
ordered_clusters
#
is.data.frame(z_scores_avg)
#
# Convert z_scores_avg to matrix
z_scores_avg <- as.matrix(z_scores_avg)
print(class(z_scores_avg))  # Should print "matrix"
colnames(z_scores_avg)
#
# Define color mapping: blue for low, white for mid, red for high
col_fun <- colorRamp2(c(min(z_scores_avg), 0, max(z_scores_avg)), c("deepskyblue3", "black", "red"))
#
# Define annotation colors
cluster_colors <- c("1" = "red", "2" = "blue", "3" = "green", 
                    "4" = "purple", "5" = "coral", "6" = "cyan", "7" = "cadetblue")

# Define row annotation with a visible legend
left_annotation = rowAnnotation(
  Cluster = anno_simple(row_clusters_by_height, 
                        col = cluster_colors, 
                        width = unit(4, "mm")),
  annotation_name_gp = gpar(fontsize = 10, col = "Azure4")
)
#
# Create heatmap
ht <- Heatmap(z_scores_avg, 
              name = "Intensity",
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 10, col = "Azure4"),
                grid_width = unit(3, "mm"),
                labels_gp = gpar(fontsize = 8)),
              show_row_names = FALSE, 
              show_column_names = TRUE, 
              column_names_rot = 45, 
              column_names_gp = gpar(fontsize = 10, col = "black", alpha = 1), 
              cluster_columns = FALSE, 
              row_dend_reorder = TRUE, 
              column_dend_reorder = FALSE, 
              col = col_fun, 
              width = ncol(z_scores) * unit(9, "mm"),
              height = nrow(z_scores) * unit(0.04, "mm"),
              row_dend_side = "left", 
              row_dend_width = unit(25, "mm"),  
              column_order = column_order, 
              gap = unit(0, "mm"),
              cluster_rows = TRUE,
              row_dend_gp = gpar(col = "azure3", lwd = 1),
              left_annotation = left_annotation)

ht
#
# Manually define the legend
cluster_legend <- Legend(
  labels = names(cluster_colors), 
  title = "Cluster",
  legend_gp = gpar(fill = cluster_colors),
  title_gp = gpar(fontsize = 10, col = "Azure4"),
  grid_width = unit(3.5, "mm"),
  labels_gp = gpar(fontsize = 8)
)
#
# Positioning the legend
pushViewport(viewport(
  x = 0.753,  # 0 =left, 1 = right, every time clear the plot first
  y = 0.70,   # 0 = bottom, 1 = top
  width = 0.15,  
  height = 0.4,  
  just = c("left", "center") 
))
# Adjust the legend position within the viewport
draw(cluster_legend, x = unit(0, "npc"), y = unit(0.5, "npc"))
#
# Exit the viewport
popViewport()
#

#######
# Individual co-expression groups
#
# Create a named vector of cluster colors based on the cluster assignments
cluster_color_mapping <- cluster_colors[as.character(row_clusters_by_height)]
#
# Cluster red
#
# Extract rows that belong to the "red" cluster
red_cluster_rows <- which(cluster_color_mapping == "red")

# Subset the z_scores_avg matrix to keep only the red cluster rows
red_cluster_data <- z_scores_avg[red_cluster_rows, ]

# Optionally, view the subset of the data for the red cluster
print(red_cluster_data)

# Convert the red_cluster_data matrix to a data frame
red_cluster_data_df <- as.data.frame(red_cluster_data)

# Add Gene column (row names) as a separate column
red_cluster_data_df$Gene <- rownames(red_cluster_data_df)

# Reshape the data into long format using pivot_longer
red_cluster_data_long <- red_cluster_data_df %>%
  pivot_longer(cols = -Gene,      # Reshape to long format (all columns except 'Gene')
               names_to = "Treatments",   # Column names (N0, N5, N25) become values in a new column 'Treatments'
               values_to = "Gene_Count")  # The values from N0, N5, N25 are placed in a new column 'Gene_Count'

# Check the reshaped data
head(red_cluster_data_long)
#
# Ensure Treatment is a factor with the specified order
red_cluster_data_long <- red_cluster_data_long %>%
  mutate(Treatments = factor(Treatments, levels = sample_order))
red_cluster_data_long
red_DEGs <- (nrow(red_cluster_data_long)/3)
red_DEGs
# Save the red cluster data to a CSV file
write.csv(red_cluster_data_long, file = "red_cluster_data.csv", row.names = FALSE)
#
# Plot red cluster gene counts across treatments with uniform line color
red_cluster_plot <- ggplot(red_cluster_data_long, aes(x = Treatments, y = Gene_Count, group = Gene)) +
  geom_line(color = "red", linewidth = 0.1) +   # Set all lines to 'red'
  #geom_point(color = "red") +  # Optional: Uncomment if you want points along the lines
  labs(title = paste0("Cluster 1 (", red_DEGs, ")")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5,color = "red",size = 10),
    #plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.text.x = element_text(angle = 45, hjust = 1, color = "Azure4", size = 10),
    axis.text.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", color = "grey")
  )

# Display the plot
red_cluster_plot
###
# Cluster blue
#
# Extract rows that belong to the "blue" cluster
blue_cluster_rows <- which(cluster_color_mapping == "blue")

# Subset the z_scores_avg matrix to keep only the blue cluster rows
blue_cluster_data <- z_scores_avg[blue_cluster_rows, ]

# Optionally, view the subset of the data for the blue cluster
print(blue_cluster_data)

# Convert the blue_cluster_data matrix to a data frame
blue_cluster_data_df <- as.data.frame(blue_cluster_data)

# Add Gene column (row names) as a separate column
blue_cluster_data_df$Gene <- rownames(blue_cluster_data_df)

# Reshape the data into long format using pivot_longer
blue_cluster_data_long <- blue_cluster_data_df %>%
  pivot_longer(cols = -Gene,      # Reshape to long format (all columns except 'Gene')
               names_to = "Treatments",   # Column names (N0, N5, N25) become values in a new column 'Treatments'
               values_to = "Gene_Count")  # The values from N0, N5, N25 are placed in a new column 'Gene_Count'

# Check the reshaped data
head(blue_cluster_data_long)
#
# Ensure Treatment is a factor with the specified order
blue_cluster_data_long <- blue_cluster_data_long %>%
  mutate(Treatments = factor(Treatments, levels = sample_order))
blue_cluster_data_long
blue_DEGs <- (nrow(blue_cluster_data_long)/3)
blue_DEGs
# Save the blue cluster data to a CSV file
write.csv(blue_cluster_data_long, file = "blue_cluster_data.csv", row.names = FALSE)
#
# Plot blue cluster gene counts across treatments with uniform line color
blue_cluster_plot <- ggplot(blue_cluster_data_long, aes(x = Treatments, y = Gene_Count, group = Gene)) +
  geom_line(color = "blue", linewidth = 0.1) +   # Set all lines to 'blue'
  #geom_point(color = "blue") +  # Optional: Uncomment if you want points along the lines
  labs(title = paste0("Cluster 2 (", blue_DEGs, ")")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5,color = "blue",size = 10),
    #plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.text.x = element_text(angle = 45, hjust = 1, color = "Azure4", size = 10),
    axis.text.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", color = "grey")
  )

# Display the plot
blue_cluster_plot
###
# Cluster green
#
# Extract rows that belong to the "green" cluster
green_cluster_rows <- which(cluster_color_mapping == "green")

# Subset the z_scores_avg matrix to keep only the green cluster rows
green_cluster_data <- z_scores_avg[green_cluster_rows, ]

# Optionally, view the subset of the data for the green cluster
print(green_cluster_data)

# Convert the green_cluster_data matrix to a data frame
green_cluster_data_df <- as.data.frame(green_cluster_data)

# Add Gene column (row names) as a separate column
green_cluster_data_df$Gene <- rownames(green_cluster_data_df)

# Reshape the data into long format using pivot_longer
green_cluster_data_long <- green_cluster_data_df %>%
  pivot_longer(cols = -Gene,      # Reshape to long format (all columns except 'Gene')
               names_to = "Treatments",   # Column names (N0, N5, N25) become values in a new column 'Treatments'
               values_to = "Gene_Count")  # The values from N0, N5, N25 are placed in a new column 'Gene_Count'

# Check the reshaped data
head(green_cluster_data_long)
#
# Ensure Treatment is a factor with the specified order
green_cluster_data_long <- green_cluster_data_long %>%
  mutate(Treatments = factor(Treatments, levels = sample_order))
green_cluster_data_long
green_DEGs <- (nrow(green_cluster_data_long)/3)
green_DEGs
# Save the green cluster data to a CSV file
write.csv(green_cluster_data_long, file = "green_cluster_data.csv", row.names = FALSE)
#
# Plot green cluster gene counts across treatments with uniform line color
green_cluster_plot <- ggplot(green_cluster_data_long, aes(x = Treatments, y = Gene_Count, group = Gene)) +
  geom_line(color = "green", linewidth = 0.1) +   # Set all lines to 'green'
  #geom_point(color = "green") +  # Optional: Uncomment if you want points along the lines
  labs(title = paste0("Cluster 3 (", green_DEGs, ")")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5,color = "green",size = 10),
    #plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.text.x = element_text(angle = 45, hjust = 1, color = "Azure4", size = 10),
    axis.text.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", color = "grey")
  )

# Display the plot
green_cluster_plot
###
# Cluster purple
#
# Extract rows that belong to the "purple" cluster
purple_cluster_rows <- which(cluster_color_mapping == "purple")

# Subset the z_scores_avg matrix to keep only the purple cluster rows
purple_cluster_data <- z_scores_avg[purple_cluster_rows, ]

# Optionally, view the subset of the data for the purple cluster
print(purple_cluster_data)

# Convert the purple_cluster_data matrix to a data frame
purple_cluster_data_df <- as.data.frame(purple_cluster_data)

# Add Gene column (row names) as a separate column
purple_cluster_data_df$Gene <- rownames(purple_cluster_data_df)

# Reshape the data into long format using pivot_longer
purple_cluster_data_long <- purple_cluster_data_df %>%
  pivot_longer(cols = -Gene,      # Reshape to long format (all columns except 'Gene')
               names_to = "Treatments",   # Column names (N0, N5, N25) become values in a new column 'Treatments'
               values_to = "Gene_Count")  # The values from N0, N5, N25 are placed in a new column 'Gene_Count'

# Check the reshaped data
head(purple_cluster_data_long)
#
# Ensure Treatment is a factor with the specified order
purple_cluster_data_long <- purple_cluster_data_long %>%
  mutate(Treatments = factor(Treatments, levels = sample_order))
purple_cluster_data_long
purple_DEGs <- (nrow(purple_cluster_data_long)/3)
purple_DEGs
# Save the purple cluster data to a CSV file
write.csv(purple_cluster_data_long, file = "purple_cluster_data.csv", row.names = FALSE)
#
# Plot purple cluster gene counts across treatments with uniform line color
purple_cluster_plot <- ggplot(purple_cluster_data_long, aes(x = Treatments, y = Gene_Count, group = Gene)) +
  geom_line(color = "purple", linewidth = 0.1) +   # Set all lines to 'purple'
  #geom_point(color = "purple") +  # Optional: Uncomment if you want points along the lines
  labs(title = paste0("Cluster 4 (", purple_DEGs, ")")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5,color = "purple",size = 10),
    #plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.text.x = element_text(angle = 45, hjust = 1, color = "Azure4", size = 10),
    axis.text.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", color = "grey")
  )

# Display the plot
purple_cluster_plot
###
# Cluster coral
#
# Extract rows that belong to the "coral" cluster
coral_cluster_rows <- which(cluster_color_mapping == "coral")

# Subset the z_scores_avg matrix to keep only the coral cluster rows
coral_cluster_data <- z_scores_avg[coral_cluster_rows, ]

# Optionally, view the subset of the data for the coral cluster
print(coral_cluster_data)

# Convert the coral_cluster_data matrix to a data frame
coral_cluster_data_df <- as.data.frame(coral_cluster_data)

# Add Gene column (row names) as a separate column
coral_cluster_data_df$Gene <- rownames(coral_cluster_data_df)

# Reshape the data into long format using pivot_longer
coral_cluster_data_long <- coral_cluster_data_df %>%
  pivot_longer(cols = -Gene,      # Reshape to long format (all columns except 'Gene')
               names_to = "Treatments",   # Column names (N0, N5, N25) become values in a new column 'Treatments'
               values_to = "Gene_Count")  # The values from N0, N5, N25 are placed in a new column 'Gene_Count'

# Check the reshaped data
head(coral_cluster_data_long)
#
# Ensure Treatment is a factor with the specified order
coral_cluster_data_long <- coral_cluster_data_long %>%
  mutate(Treatments = factor(Treatments, levels = sample_order))
coral_cluster_data_long
coral_DEGs <- (nrow(coral_cluster_data_long)/3)
coral_DEGs
# Save the coral cluster data to a CSV file
write.csv(coral_cluster_data_long, file = "coral_cluster_data.csv", row.names = FALSE)
#
# Plot coral cluster gene counts across treatments with uniform line color
coral_cluster_plot <- ggplot(coral_cluster_data_long, aes(x = Treatments, y = Gene_Count, group = Gene)) +
  geom_line(color = "coral", linewidth = 0.1) +   # Set all lines to 'coral'
  #geom_point(color = "coral") +  # Optional: Uncomment if you want points along the lines
  labs(title = paste0("Cluster 5 (", coral_DEGs, ")")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5,color = "coral",size = 10),
    #plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.text.x = element_text(angle = 45, hjust = 1, color = "Azure4", size = 10),
    axis.text.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", color = "grey")
  )

# Display the plot
coral_cluster_plot
###
# Cluster cyan
#
# Extract rows that belong to the "cyan" cluster
cyan_cluster_rows <- which(cluster_color_mapping == "cyan")

# Subset the z_scores_avg matrix to keep only the cyan cluster rows
cyan_cluster_data <- z_scores_avg[cyan_cluster_rows, ]

# Optionally, view the subset of the data for the cyan cluster
print(cyan_cluster_data)

# Convert the cyan_cluster_data matrix to a data frame
cyan_cluster_data_df <- as.data.frame(cyan_cluster_data)

# Add Gene column (row names) as a separate column
cyan_cluster_data_df$Gene <- rownames(cyan_cluster_data_df)

# Reshape the data into long format using pivot_longer
cyan_cluster_data_long <- cyan_cluster_data_df %>%
  pivot_longer(cols = -Gene,      # Reshape to long format (all columns except 'Gene')
               names_to = "Treatments",   # Column names (N0, N5, N25) become values in a new column 'Treatments'
               values_to = "Gene_Count")  # The values from N0, N5, N25 are placed in a new column 'Gene_Count'

# Check the reshaped data
head(cyan_cluster_data_long)
#
# Ensure Treatment is a factor with the specified order
cyan_cluster_data_long <- cyan_cluster_data_long %>%
  mutate(Treatments = factor(Treatments, levels = sample_order))
cyan_cluster_data_long
cyan_DEGs <- (nrow(cyan_cluster_data_long)/3)
cyan_DEGs
# Save the cyan cluster data to a CSV file
write.csv(cyan_cluster_data_long, file = "cyan_cluster_data.csv", row.names = FALSE)
#
# Plot cyan cluster gene counts across treatments with uniform line color
cyan_cluster_plot <- ggplot(cyan_cluster_data_long, aes(x = Treatments, y = Gene_Count, group = Gene)) +
  geom_line(color = "cyan", linewidth = 0.1) +   # Set all lines to 'cyan'
  #geom_point(color = "cyan") +  # Optional: Uncomment if you want points along the lines
  labs(title = paste0("Cluster 6 (", cyan_DEGs, ")")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5,color = "cyan",size = 10),
    #plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "Azure4", size = 10),
    #axis.text.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", color = "grey")
  )

# Display the plot
cyan_cluster_plot
###
# Cluster cadetblue
#
# Extract rows that belong to the "cadetblue" cluster
cadetblue_cluster_rows <- which(cluster_color_mapping == "cadetblue")

# Subset the z_scores_avg matrix to keep only the cadetblue cluster rows
cadetblue_cluster_data <- z_scores_avg[cadetblue_cluster_rows, ]

# Optionally, view the subset of the data for the cadetblue cluster
print(cadetblue_cluster_data)

# Convert the cadetblue_cluster_data matrix to a data frame
cadetblue_cluster_data_df <- as.data.frame(cadetblue_cluster_data)

# Add Gene column (row names) as a separate column
cadetblue_cluster_data_df$Gene <- rownames(cadetblue_cluster_data_df)

# Reshape the data into long format using pivot_longer
cadetblue_cluster_data_long <- cadetblue_cluster_data_df %>%
  pivot_longer(cols = -Gene,      # Reshape to long format (all columns except 'Gene')
               names_to = "Treatments",   # Column names (N0, N5, N25) become values in a new column 'Treatments'
               values_to = "Gene_Count")  # The values from N0, N5, N25 are placed in a new column 'Gene_Count'

# Check the reshaped data
head(cadetblue_cluster_data_long)
#
# Ensure Treatment is a factor with the specified order
cadetblue_cluster_data_long <- cadetblue_cluster_data_long %>%
  mutate(Treatments = factor(Treatments, levels = sample_order))
cadetblue_cluster_data_long
cadetblue_DEGs <- (nrow(cadetblue_cluster_data_long)/3)
cadetblue_DEGs
# Save the cadetblue cluster data to a CSV file
write.csv(cadetblue_cluster_data_long, file = "cadetblue_cluster_data.csv", row.names = FALSE)
#
# Plot cadetblue cluster gene counts across treatments with uniform line color
cadetblue_cluster_plot <- ggplot(cadetblue_cluster_data_long, aes(x = Treatments, y = Gene_Count, group = Gene)) +
  geom_line(color = "cadetblue", linewidth = 0.1) +   # Set all lines to 'cadetblue'
  #geom_point(color = "cadetblue") +  # Optional: Uncomment if you want points along the lines
  labs(title = paste0("Cluster 7 (", cadetblue_DEGs, ")")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5,color = "cadetblue",size = 10),
    #plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "Azure4", size = 10),
    #axis.text.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", color = "grey")
  )

# Display the plot
cadetblue_cluster_plot
###
#
###
### Combine the plots
# Combine the two plots side-by-side with spacing for Y-axis label
combined_plot <- plot_grid(
  red_cluster_plot,
  blue_cluster_plot,
  green_cluster_plot,
  purple_cluster_plot,
  coral_cluster_plot,
  cyan_cluster_plot,
  cadetblue_cluster_plot,
  #labels = c("A", "B", "C", "D", "E", "F", "G"),   # Optional plot labels
  ncol = 2,              # Arrange vertically
  nrow = 4,
  align = "hv"           # Align both horizontally and vertically
)
combined_plot
#
# Add a common Y-axis title outside the plot labels 
final_plot <- ggdraw() +
  draw_label("Gene Count", x = 0.06, y = 0.5, angle = 90, vjust = 1.5, size = 12, color = "Azure4") + 
  draw_plot(combined_plot, x = 0.05, y = 0.05, width = 0.95, height = 0.90) +
  draw_label("Treatments", x = 0.5, y = 0, vjust = -1, size = 12, color = "Azure4")

# Display final plot
final_plot
###







