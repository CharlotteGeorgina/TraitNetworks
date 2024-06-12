# Load necessary libraries
library(dplyr)
library(tidyr)
library(Hmisc)
library(corrplot)
library(tidyverse)
library(ggplot2)
library(sqldf)
library(cooccur)
library(reshape2)
library(circleplot)
library(circlize)
library(graph4lg)
library(igraph)

# Read and preprocess data
CWM <- read.csv("Data/CWM_JP_Fish_Pres_Ab.csv") %>% rename(Site = X)
CWM <- CWM %>% relocate(Clusters, .before = Site) %>% column_to_rownames(., var = 'Site')
colnames(CWM) <- c("Cluster", "ML_L", "ML_M", "ML_S", "ML_VS", "PLD_L", "PLD_M", "PLD_S", "PLD_VL", "T_C", "T_D", "T_H", "T_O", "T_Pis", "T_Pla", "T_Pre", "P_B", "P_CA", "P_D", "P_EA", "P_P", "P_RP", "P_SA", "P_SB", "P_UB", "R_B", "R_D", "R_LB", "R_N", "R_S")
CWM$Cluster[CWM$Cluster == "ColdTropical"] <- "WarmSubtropical"
CWM$Cluster[CWM$Cluster == "Subtropical"] <- "ColdSubtropical"

# Create a data frame for regions
countR <- data.frame(Cluster = c("Tropical", "WarmSubtropical", "ColdSubtropical", "Temperate"))

# Initialize the result data frame
all_metrics <- list()
bootstrap_iterations <- 1000

# Bootstrapping and network metric calculation
for (i in 1:nrow(countR)) {
  Cluster <- countR[i, 1]
  CWM_i <- CWM[CWM$Cluster == Cluster,] %>% select(-Cluster)
  
  # Initialize data frames to store bootstrap results
  bootstrap_metrics <- data.frame(DegCents=numeric(), BetwCentr=numeric(), EigCentr=numeric(), EdgeDen=numeric(), NoModule=numeric(), NetModularity=numeric())
  
  for (j in 1:bootstrap_iterations) {
    # Bootstrapping
    CWM_bootstrap <- CWM_i[sample(nrow(CWM_i), replace = TRUE), ]
    
    # Calculate correlation matrix and significance
    CWM_Cor <- rcorr(as.matrix(CWM_bootstrap))
    CWM_prob <- ifelse(CWM_Cor$P < 0.05, 1, 0)
    CWM_eff <- CWM_Cor$r
    SignMatrix <- CWM_prob * CWM_eff
    SignMatrix[is.na(SignMatrix)] <- 0
    SignMatrix <- ifelse(SignMatrix > 0.1, 1, 0)
    SignMatrix2 <- SignMatrix[rowSums(SignMatrix == 0) != ncol(SignMatrix), colSums(SignMatrix == 0) != nrow(SignMatrix)]
    
    # Generate network graph
    TGraph <- graph_from_adjacency_matrix(SignMatrix2, mode = "undirected")
    WT_partition <- cluster_walktrap(TGraph, weights = NULL)
    
    # Calculate network-wide metrics
    DegCents <- centr_clo(TGraph, mode = "all")$centralization
    BetwCentr <- centr_betw(TGraph, directed = FALSE)$centralization
    EigCentr <- centr_eigen(TGraph, directed = FALSE)$centralization
    EdgeDen <- edge_density(TGraph, loops = FALSE)
    NoModule <- length(WT_partition)
    NetModularity <- modularity(TGraph, membership(WT_partition))
    
    # Store metrics
    bootstrap_metrics[j, ] <- c(DegCents, BetwCentr, EigCentr, EdgeDen, NoModule, NetModularity)
  }
  
  # Store all bootstrap metrics for the current cluster
  all_metrics[[Cluster]] <- bootstrap_metrics
}

# Calculate mean and standard deviation for each region and metric
mean_sd_metrics <- data.frame()
for (cluster in names(all_metrics)) {
  cluster_metrics <- all_metrics[[cluster]]
  means <- colMeans(cluster_metrics)
  sds <- apply(cluster_metrics, 2, sd)
  mean_sd_metrics <- rbind(mean_sd_metrics, cbind(Cluster=cluster, t(means), t(sds)))
}
colnames(mean_sd_metrics) <- c("Cluster", "Mean_DegCents", "Mean_BetwCentr", "Mean_EigCentr", "Mean_EdgeDen", "Mean_NoModule", "Mean_NetModularity", "SD_DegCents", "SD_BetwCentr", "SD_EigCentr", "SD_EdgeDen", "SD_NoModule", "SD_NetModularity")

# Save mean and SD results
write.csv(mean_sd_metrics, 'Results/Eocene_Oligocene_Mean_SD_22_03_23.csv')

# Perform t-tests between tropical and temperate groups for each metric
t_test_results <- data.frame(Metric=character(), p_value=numeric(), stringsAsFactors=FALSE)
metrics <- colnames(all_metrics[[1]])

for (metric in metrics) {
  tropical_values <- all_metrics[["Tropical"]][[metric]]
  temperate_values <- all_metrics[["Temperate"]][[metric]]
  t_test <- t.test(tropical_values, temperate_values)
  t_test_results <- rbind(t_test_results, data.frame(Metric=metric, p_value=t_test$p.value))
}

# Save t-test results
write.csv(t_test_results, 'Results/Eocene_Oligocene_T_Test_Results_22_03_23.csv')
