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
library(circlize)
library(graph4lg)
library(igraph)
library(psych)

#### Data #### 
CWM <- read.csv("Data/CWM_JPFish_Biomass_TropTemp.csv")
CWM <- CWM %>% rename(Site = X)

#organise and rename
CWM <- CWM %>% relocate(Cluster, .before = Site)
CWM <- CWM %>% column_to_rownames(., var = 'Site')
colnames(CWM) <- c("Cluster", "ML_Large", "ML_Medium", "ML_Small",
                   "ML_V.Small", "PLD_Long", "PLD_Medium", "PLD_Short",
                   "PLD_V.Long", "T_Corrallivore", "T_Detritivore",
                   "T_Herbivore", "T_Omnivore", "T_Piscivore",
                   "T_Planktivore", "T_Predator", "P_Benthic",
                   "P_AnthozoanA", "P_Demersal", 
                   "P_Pelagic", "P_ReefPelagic", "P_SandA", "P_SubBenthic",
                   "P_UpperBenthic", "R_Brooders", "R_Demersal",
                   "R_livebearers", "R_Nesters", "R_Scatterers")

##### Make a counter file for regions for looping:  ----
countR <- data.frame(Cluster = c("Tropical", "Temperate"))

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

# Combine metrics from all_metrics into a single data frame with Cluster indicator
combined_metrics <- NULL

for (cluster in names(all_metrics)) {
  cluster_metrics <- all_metrics[[cluster]]
  cluster_df <- as.data.frame(cluster_metrics)
  cluster_df$Cluster <- cluster
  combined_metrics <- rbind(combined_metrics, cluster_df)
}

combined_metrics <- combined_metrics %>% select(Cluster, everything())

# Save combined metrics to CSV
write.csv(combined_metrics, 'Results/JP_Fish_Combined_Metrics_Tropical_Temperate.csv', row.names = FALSE)

# Group by Cluster and calculate metrics
metrics_summary <- combined_metrics %>%
  group_by(Cluster) %>%
  summarise(across(everything(), list(
    Mean = mean, 
    SD = sd, 
    Median = median,
    IQR = IQR
  )))

# Flatten the grouped data frame for easier viewing
metrics_summary <- as.data.frame(metrics_summary)

# Save metrics summary to CSV
write.csv(metrics_summary, 'Results/JP_Fish_Metrics_Summary_Tropical_Temperate.csv', row.names = FALSE)

#Normality
shapiro.test(combined_metrics$DegCents)
shapiro.test(combined_metrics$BetwCentr)
shapiro.test(combined_metrics$EigCentr)
shapiro.test(combined_metrics$EdgeDen)
shapiro.test(combined_metrics$NoModule)
shapiro.test(combined_metrics$NetModularity)
#all non-noramal 


# Initialize an empty data frame to store results
wilcox_test_results <- data.frame(Metric = character(), W_statistic = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through each numeric metric
for (metric in colnames(combined_metrics)[sapply(combined_metrics, is.numeric)]) {
  # Perform Wilcoxon rank sum test
  wilcox_result <- wilcox.test(as.formula(paste(metric, "~ Cluster")), data = combined_metrics, exact = FALSE)
  
  # Extract W statistic, p-value, and metric name
  metric_name <- metric
  W_statistic <- wilcox_result$statistic
  p_value <- wilcox_result$p.value
  
  # Round p-value to 3 decimal places
  p_value <- round(p_value, digits = 3)
  
  # Append results to wilcox_test_results data frame
  wilcox_test_results <- rbind(wilcox_test_results, data.frame(Metric = metric_name, W_statistic = W_statistic, p_value = p_value))
}

# Save Wilcoxon test results to CSV
write.csv(wilcox_test_results, 'Results/JpFish_Wilcoxon_Test_Results.csv', row.names = FALSE)


