# Load necessary libraries
library(corrplot)
library(tidyverse)  # Includes dplyr and other packages
library(cooccur)
library(circlize)
library(igraph)
library(Hmisc)

#### Data #### 
Pre <- read.csv("Data/Speciesbytraitmatrix_Pre.csv")
Post <- read.csv("Data/Speciesbytraitmatrix_Post.csv")


#feeding_deposit is missing from: Post
Post$feeding_deposit <- rep(0, 37)

#add source column 
Pre$source <- "Pre"
Post$source <- "Post"

# Combine the data frames and retain the 'source' column
All <- bind_rows(Pre, Post)

#rename and reorder
All <- All %>% relocate(source, .before = X)
colnames(All) <- c("Time", "taxon", "F_D", "F_G", "F_M", "F_P", "F_S", "M_Fac",
                   "M_fast", "M_non", "M_slow", "S_G", "S_L", "S_S", "S_M", "S_T", "S_VL",
                   "T_E", "T_I", "T_P")
##### Make a counter file for times for looping #####

# Create a data frame to loop through times
countR <- data.frame(Time = c("Pre", "Post"))

##### Initialize a dataframe to store network metrics #####

Reg_metrix <- data.frame(
  Time = character(),
  DegCents = numeric(),
  BetwCentr = numeric(),
  EigCentr = numeric(),
  EdgeDen = numeric(),
  NoModule = numeric(),
  NetModularity = numeric(),
  stringsAsFactors = FALSE
)

##### Parameters for permutation testing #####

num_permutations <- 100  # Number of permutations

##### Loop through Regions, analyzing each time period #####

for (i in 1:nrow(countR)) {
  # Get the current time period
  Time <- countR[i, "Time"]
  
  # Subset the data by Time
  All_i <- All %>% filter(Time == Time) %>% select(-Time)  # Remove Time column
  
  # Convert to matrix for correlation analysis
  All_mat <- as.matrix(All_i[-1])  # Exclude taxon column
  
  # Calculate correlation matrix
  tra_Cor <- rcorr(All_mat)
  
  ## Probability of pairwise co-occurrence
  tra_prob <- ifelse(tra_Cor$P < 0.05, 1, 0)
  out_prob_s <- tra_prob
  write.csv(out_prob_s, file = paste0("Results/", Time, "/outprob_s", Time, ".csv"))
  
  ## Effect matrix
  effectM <- tra_Cor$r
  write.csv(effectM, file = paste0("Results/", Time, "/outeff_", Time, ".csv"))
  
  ## Generate a matrix where links are significant
  SignMatrix <- out_prob_s * effectM
  write.csv(SignMatrix, file = paste0("Results/", Time, "/SignMatrix_", Time, ".csv"))
  
  ## Binarize SignMatrix
  SignMatrix <- ifelse(SignMatrix > 0.1, 1, 0)
  SignMatrix2 <- SignMatrix[rowSums(SignMatrix == 0) != ncol(SignMatrix), colSums(SignMatrix == 0) != nrow(SignMatrix)]
  write.csv(SignMatrix2, file = paste0("Results/", Time, "/SignMatrix2_", Time, ".csv"))
  
  ### Graph-theoretic metrics using igraph ###
  
  # Original network
  TGraph <- graph_from_adjacency_matrix(SignMatrix2, mode = "undirected")
  
  # Initialize vectors to store metric results
  DegCents <- numeric(num_permutations + 1)
  BetwCentr <- numeric(num_permutations + 1)
  EigCentr <- numeric(num_permutations + 1)
  EdgeDen <- numeric(num_permutations + 1)
  NoModule <- numeric(num_permutations + 1)
  NetModularity <- numeric(num_permutations + 1)
  
  # Compute metrics for original network
  DegCents[1] <- centr_clo(TGraph, mode = "all")$centralization
  BetwCentr[1] <- centr_betw(TGraph, directed = FALSE)$centralization
  EigCentr[1] <- centr_eigen(TGraph, directed = FALSE)$centralization
  EdgeDen[1] <- edge_density(TGraph, loops = FALSE)
  WT_partition <- cluster_walktrap(TGraph, weights = NULL)
  NoModule[1] <- length(WT_partition)
  NetModularity[1] <- modularity(TGraph, membership(WT_partition))
  
  # Perform permutation testing
  for (perm in 1:num_permutations) {
    # Generate a null network (randomized network)
    null_network <- rewire(TGraph, with = keeping_degseq(niter = 10))
    
    # Compute metrics for null network
    DegCents[perm + 1] <- centr_clo(null_network, mode = "all")$centralization
    BetwCentr[perm + 1] <- centr_betw(null_network, directed = FALSE)$centralization
    EigCentr[perm + 1] <- centr_eigen(null_network, directed = FALSE)$centralization
    EdgeDen[perm + 1] <- edge_density(null_network, loops = FALSE)
    WT_partition_null <- cluster_walktrap(null_network, weights = NULL)
    NoModule[perm + 1] <- length(WT_partition_null)
    NetModularity[perm + 1] <- modularity(null_network, membership(WT_partition_null))
  }
  
  # Store average metrics across permutations
  Reg_metrix[nrow(Reg_metrix) + 1, ] <- c(Time,
                                          mean(DegCents),
                                          mean(BetwCentr),
                                          mean(EigCentr),
                                          mean(EdgeDen),
                                          mean(NoModule),
                                          mean(NetModularity))
  
  # Compare observed metrics to null distributions
  
  # Perform one-sided tests for metrics where higher values are of interest
  p_vals <- numeric(length(Reg_metrix) - 1)  # Initialize vector for p-values
  
  # Calculate p-values for each metric
  for (j in 2:ncol(Reg_metrix)) {  # Skip the 'Time' column
    observed <- Reg_metrix[nrow(Reg_metrix), j]
    null_distribution <- Reg_metrix[1:num_permutations, j]
    p_vals[j - 1] <- mean(null_distribution >= observed)
  }
  
  # Adjust for multiple comparisons using Bonferroni correction
  alpha <- 0.05  # Desired overall significance level
  num_tests <- length(p_vals)
  adjusted_alpha <- alpha / num_tests
  
  # Perform Bonferroni correction
  significant <- p_vals <= adjusted_alpha
  
  # Print or store adjusted p-values
  cat(paste0("Time period: ", Time, "\n"))
  cat("Observed vs Null p-values (Bonferroni corrected):\n")
  print(p_vals)
  
  # Print significant results
  cat("Significant metrics:\n")
  print(significant)
  
  # Save the metrics and results
  write.csv(Reg_metrix, file = 'Results/Reg_metrix_13_06_24.csv')
}

