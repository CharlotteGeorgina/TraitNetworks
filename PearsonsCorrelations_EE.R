# Load necessary libraries
library(corrplot)
library(tidyverse)  # Includes dplyr and other packages
library(cooccur)
library(circlize)
library(Hmisc)
library(igraph)

#### Data #### 
Pre <- read.csv("Data/Speciesbytraitmatrix_Pre.csv")
Post <- read.csv("Data/Speciesbytraitmatrix_Post.csv")


#feeding_deposit is missing from: Post, Early and Late
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
##### Make a counter file for times for looping:  ----
countR <- data.frame(Time = c("Pre", "Post"))


##### loop through Regions, making eveything ----

Reg_metrix = data.frame(Time=as.character(), DegCents = numeric(), BetwCentr = numeric(),
                        EigCentr  = numeric(), EdgeDen = numeric(), NoModule = numeric(), NetModularity = numeric(), stringsAsFactors=FALSE)

for (i in 1:nrow(countR)) { #circle through the  times
  # work out time
  Time = countR[i,1]
  
  # subsetting by Time
  All_i = All[All$Time == Time,]  # select the records for the Time we want
  All_i = All_i[,-1]
  All_ii <- All_i[-1]
  row.names(All_ii) <- All_i$taxon
  tra_Cor <- rcorr(as.matrix(All_ii)) #produces correlation matrix 

  ## probtable: write results for probability of pairwise co-occurrence to output table
  tra_prob <- tra_Cor$P
  out_prob <- tra_prob; csvfile<-paste0("Results/",Time,"/outprob_",Time,".csv"); write.csv(out_prob, csvfile) 
  
  #binarise the prob table, < 0.05 = 1, >0.05 = 0
  tra_prob <- ifelse(tra_prob < 0.05, 1, 0)
  out_prob_s <- tra_prob; csvfile<-paste0("Results/",Time,"/outprob_s",Time,".csv"); write.csv(out_prob_s, csvfile)
  
  #Get correlation matrix
  tra_eff <- tra_Cor$r
  effectM = as.matrix(tra_eff); csvfile<-paste0("Results/",Time,"/outeff_",Time,".csv"); write.csv(effectM, csvfile) #Get effects table
  
  ### Generate a matrix where links are significant
  SignMatrix <- out_prob_s * effectM
  csvfile=paste0("Results/",Time,"/SignMatrix_",Time,".csv"); write.csv(SignMatrix, csvfile)
  
  #change na and Nan to 0's 
  SignMatrix[is.na(SignMatrix)] <- 0
  SignMatrix <- ifelse(SignMatrix>0.1, 1, 0)    # Binarize 
  #cut out the zero columns 
  SignMatrix2 = SignMatrix[rowSums(SignMatrix == 0) != ncol(SignMatrix), colSums(SignMatrix == 0) != nrow(SignMatrix)]
  csvfile=paste0("Results/", Time, "/SignMatrix2_",Time,".csv"); write.csv(SignMatrix2, csvfile)
  
  ### graph-theoretic metrix code ### calls 'igraph'  
  TGraph = graph_from_adjacency_matrix(SignMatrix2, mode = "undirected")
  plot(TGraph) #it will ignore the weights at this point
  
  ### cluster any substructures (most techniques dont work because of neg relationships)
  WT_partition = cluster_walktrap(TGraph, weights = NULL)
  
  ### Node/ Trait specific metrics:
  t_Degree <- igraph::degree(TGraph, mode = "all", loops = FALSE)
  t_Degree <- t_Degree / max(t_Degree)  # Normalize the degree centrality manually
  t_DegCentr <- centr_clo(TGraph, mode = "all")$res
  tModule <- data.frame(WT_partition$membership)
  colnames(tModule) <- c("tModule")
  tModule$Trait <- V(TGraph)$name
  tModule <- rowid_to_column(tModule, "ID")
  
  # trait level modularity
  t_metrix = cbind(Trait = vertex_attr(TGraph, "name"),t_Degree,t_DegCentr,tModule)
  t_metrix[is.na(t_metrix)] <- 0
  csvfile=paste0("Results/",Time,"/t_metrix_2_",Time,".csv"); write.csv(t_metrix, csvfile)
  
  ### Network-wide metrics
  DegCents = centr_clo(TGraph, mode = "all")$centralization
  BetwCentr = centr_betw(TGraph, directed = FALSE)$centralization
  EigCentr = centr_eigen(TGraph, directed = FALSE)$centralization
  EdgeDen = edge_density(TGraph, loops = FALSE)
  NoModule = length(WT_partition)
  NetModularity = modularity(TGraph,membership(WT_partition))
  # Time level summary table (DegCents, BetwCentr, EigCentr, NoModule, NetModularity) 
  Reg_metrix[nrow(Reg_metrix)+1,] = c(Time, DegCents, BetwCentr, EigCentr, EdgeDen, NoModule, NetModularity)
  
  #Plot and save
  head(V(TGraph))
 
   V(TGraph)$color <- ifelse(V(TGraph)$name == "F_D", "dodgerblue4", ifelse(V(TGraph)$name == "F_G", "dodgerblue3", 
                                                                            ifelse(V(TGraph)$name == "F_M", "dodgerblue1", 
                                                                                   ifelse(V(TGraph)$name == "F_P", "deepskyblue2", 
                                                                                          ifelse(V(TGraph)$name == "F_S", "deepskyblue",
                                                                                                 ifelse(V(TGraph)$name == "M_Fac", "mediumpurple4",
                                                                                                        ifelse(V(TGraph)$name == "M_fast", "mediumpurple2",
                                                                                                               ifelse(V(TGraph)$name == "M_slow", "mediumpurple",
                                                                                                                      ifelse(V(TGraph)$name == "M_non", "mediumorchid",
                                                                                                                             ifelse(V(TGraph)$name == "S_G", "honeydew4",
                                                                                                                                    ifelse(V(TGraph)$name == "S_L", "honeydew2",
                                                                                                                                           ifelse(V(TGraph)$name == "S_S", "mistyrose2",
                                                                                                                                                  ifelse(V(TGraph)$name == "S_M", "mistyrose",
                                                                                                                                                         ifelse(V(TGraph)$name == "S_T", "antiquewhite3",
                                                                                                                                                                ifelse(V(TGraph)$name == "S_VL", "antiquewhite",
                                                                                                                                                                       ifelse(V(TGraph)$name == "T_E", "darkorange3",
                                                                                                                                                                              ifelse(V(TGraph)$name =="T_I", "darkorange1",
                                                                                                                                                                                     ifelse(V(TGraph)$name =="T_P", "darkgoldenrod1", "red"))))))))))))))))))
   

   
   #Assign group ideas = modules 
   group_ids <- lapply(tModule %>% split(.$tModule), function(grp) { grp$ID })
   group_ids
  # set colour
   group_colour <- c("gray", "gray", "gray", "gray", "gray")
   #same but add transparency
   group_colour_fill <- paste0(group_colour, '95')
   
   #set node sizes         
  V(TGraph)$size<-t_metrix$t_Degree*100
  #set layout
  layout = layout_with_kk(TGraph)

  ##to save plot 
  pdf(file=paste0("Results/Graphs/Network_",Time,"2.pdf"),
      width=8, height=6)
  plot(TGraph, 
       edge.width=2,
       edge.color="black",
       vertex.label.color="black",
       mark.groups = group_ids,
       mark.col = group_colour_fill,
       mark.border = group_colour,
       layout = layout, hovermode='closest')
  dev.off()



}


#save the Cooccurs and Metrix from the entire loop
write.csv(Reg_metrix, 'Results/Reg_metrix_14_06_24.csv')


#### Null Network Comparisons ####
# Function to calculate network metrics including node-specific metrics
# Function to calculate network metrics including node-specific metrics
calculate_metrics <- function(TGraph, num_permutations) {
  # Initialize vectors to store metric results
  DegCents <- numeric(num_permutations + 1)
  BetwCentr <- numeric(num_permutations + 1)
  EigCentr <- numeric(num_permutations + 1)
  EdgeDen <- numeric(num_permutations + 1)
  NoModule <- numeric(num_permutations + 1)
  NetModularity <- numeric(num_permutations + 1)
  
  # Compute metrics for original network
  DegCents[1] <- centr_clo(TGraph)$centralization
  BetwCentr[1] <- centr_betw(TGraph, directed = FALSE)$centralization
  EigCentr[1] <- centr_eigen(TGraph, directed = FALSE)$centralization
  EdgeDen[1] <- edge_density(TGraph, loops = FALSE)
  WT_partition <- cluster_walktrap(TGraph, weights = NULL)
  NoModule[1] <- length(WT_partition)
  NetModularity[1] <- modularity(TGraph, membership(WT_partition))
  
  # Node-specific metrics for original network
  t_Degree <- igraph::degree(TGraph, mode = "all", loops = FALSE, normalized = T)
  t_DegCentr <- centr_clo(TGraph)$res
  tModule <- data.frame(WT_partition$membership)
  colnames(tModule) <- c("tModule")
  tModule$Trait <- V(TGraph)$name
  tModule <- rowid_to_column(tModule, "ID")
  t_metrix <- cbind(Trait = vertex_attr(TGraph, "name"), t_Degree, t_DegCentr, tModule)
  t_metrix[is.na(t_metrix)] <- 0
  
  # Initialize lists to store null distributions for node-specific metrics
  null_t_Degree <- matrix(0, nrow = length(V(TGraph)), ncol = num_permutations)
  null_t_DegCentr <- matrix(0, nrow = length(V(TGraph)), ncol = num_permutations)
  
  # Perform permutation testing
  for (perm in 1:num_permutations) {
    # Generate a null network (ER random network)
    null_network <- erdos.renyi.game(vcount(TGraph), p = edge_density(TGraph), type = "gnp", directed = FALSE)
    
    # Compute metrics for null network
    DegCents[perm + 1] <- centr_clo(null_network)$centralization
    BetwCentr[perm + 1] <- centr_betw(null_network, directed = FALSE)$centralization
    EigCentr[perm + 1] <- centr_eigen(null_network, directed = FALSE)$centralization
    EdgeDen[perm + 1] <- edge_density(null_network, loops = FALSE)
    WT_partition_null <- cluster_walktrap(null_network, weights = NULL)
    NoModule[perm + 1] <- length(WT_partition_null)
    NetModularity[perm + 1] <- modularity(null_network, membership(WT_partition_null))
    
    # Compute metrics for null network
    null_t_Degree[, perm] <- igraph::degree(null_network, mode = "all", loops = FALSE, normalized = T) / max(igraph::degree(null_network, mode = "all", loops = FALSE, normalized = T))
    null_t_DegCentr[, perm] <- centr_clo(null_network)$res
  }
  
  # Compare observed node-specific metrics to null distributions
  p_vals_Degree <- rowMeans(null_t_Degree >= t_Degree)
  p_vals_DegCentr <- rowMeans(null_t_DegCentr >= t_DegCentr)
  
  return(list(DegCents = DegCents, BetwCentr = BetwCentr, EigCentr = EigCentr, EdgeDen = EdgeDen, NoModule = NoModule, NetModularity = NetModularity, t_metrix = t_metrix, p_vals_Degree = p_vals_Degree, p_vals_DegCentr = p_vals_DegCentr))
}


# Function to save node-specific p-values
save_node_specific_pvals <- function(Time, metrics) {
  node_pvals <- data.frame(
    Trait = metrics$t_metrix$Trait,
    p_vals_Degree = metrics$p_vals_Degree,
    p_vals_DegCentr = metrics$p_vals_DegCentr
  )
  csvfile <- paste0("Results/", Time, "/node_pvals_", Time, ".csv")
  write.csv(node_pvals, csvfile)
}

# Main analysis after your original code
num_permutations <- 100

# Create graphs from the saved SignMatrix2 files
SignMatrix2_Pre <- read.csv("Results/Pre/SignMatrix2_Pre.csv", row.names = 1)
SignMatrix2_Post <- read.csv("Results/Post/SignMatrix2_Post.csv", row.names = 1)

TGraph_Pre <- graph_from_adjacency_matrix(as.matrix(SignMatrix2_Pre), mode = "undirected")
TGraph_Post <- graph_from_adjacency_matrix(as.matrix(SignMatrix2_Post), mode = "undirected")

# Calculate metrics for each time period
metrics_Pre <- calculate_metrics(TGraph_Pre, num_permutations)
metrics_Post <- calculate_metrics(TGraph_Post, num_permutations)

#store Metrics and p-values for Pre
# Extract observed metrics from metrics object
# Extract the first 6 elements (lists) from metrics_Pre
metrics_Pre_list <- metrics_Pre[1:6]

# Convert each list to a data frame
metrics_Pre_df <- do.call(rbind, metrics_Pre_list)

observed_metrics <- metrics_Pre_df[, 1]  # First column is observed

# Extract null distributions from metrics object (remaining columns)
null_distributions <- metrics_Pre_df[, -1]  # Exclude the first column

# Convert NaNs to 0 in null distributions
null_distributions[is.nan(null_distributions)] <- 0

# Calculate means of null metrics
mean_null_metrics <- colMeans(null_distributions, na.rm = TRUE)

# Append null metrics (Pre_Null or Post_Null) as rows to Reg_metrix
Reg_metrix <<- rbind(Reg_metrix,
                     data.frame(
                       Time = paste("Pre_Null", sep = ""),
                       DegCents = mean_null_metrics[1],  # Adjust indices accordingly
                       BetwCentr = mean_null_metrics[2],
                       EigCentr = mean_null_metrics[3],
                       EdgeDen = mean_null_metrics[4],
                       NoModule = mean_null_metrics[5],
                       NetModularity = mean_null_metrics[6]
                       # Add more columns as needed
                     ))

# Calculate p-values for each metric
p_vals <- numeric(length(observed_metrics))  # Initialize vector for p-values

for (j in seq_along(observed_metrics)) {
  observed <- observed_metrics[j]
  null_distribution <- null_distributions[j, ]
  p_vals[j] <- mean(null_distribution >= observed)
}

# Append p-values as rows to Reg_metrix
Reg_metrix <<- rbind(Reg_metrix,
                     data.frame(
                       Time = paste("Pre_p", sep = ""),
                       DegCents = p_vals[1],  # Adjust indices accordingly
                       BetwCentr = p_vals[2],
                       EigCentr = p_vals[3],
                       EdgeDen = p_vals[4],
                       NoModule = p_vals[5],
                       NetModularity = p_vals[6]
                       # Add more columns as needed
                     ))

# Print p-values
cat(paste0("Time period: ", Time, "\n"))
cat("Observed vs Null p-values:\n")
print(p_vals)

#store Metrics and p-values for post
# Extract observed metrics from metrics object
# Extract the first 6 elements (lists) from metrics_Pre
metrics_Post_list <- metrics_Post[1:6]

# Convert each list to a data frame
metrics_Post_df <- do.call(rbind, metrics_Post_list)

observed_metrics <- metrics_Post_df[, 1]  # First column is observed

# Extract null distributions from metrics object (remaining columns)
null_distributions <- metrics_Post_df[, -1]  # Exclude the first column

# Convert NaNs to 0 in null distributions
null_distributions[is.nan(null_distributions)] <- 0

# Calculate means of null metrics
mean_null_metrics <- colMeans(null_distributions, na.rm = TRUE)

# Append null metrics (Pre_Null or Post_Null) as rows to Reg_metrix
Reg_metrix <<- rbind(Reg_metrix,
                     data.frame(
                       Time = paste("Post_Null", sep = ""),
                       DegCents = mean_null_metrics[1],  # Adjust indices accordingly
                       BetwCentr = mean_null_metrics[2],
                       EigCentr = mean_null_metrics[3],
                       EdgeDen = mean_null_metrics[4],
                       NoModule = mean_null_metrics[5],
                       NetModularity = mean_null_metrics[6]
                       # Add more columns as needed
                     ))

# Calculate p-values for each metric
p_vals <- numeric(length(observed_metrics))  # Initialize vector for p-values

for (j in seq_along(observed_metrics)) {
  observed <- observed_metrics[j]
  null_distribution <- null_distributions[j, ]
  p_vals[j] <- mean(null_distribution >= observed)
}

# Append p-values as rows to Reg_metrix
Reg_metrix <<- rbind(Reg_metrix,
                     data.frame(
                       Time = paste("Post_p", sep = ""),
                       DegCents = p_vals[1],  # Adjust indices accordingly
                       BetwCentr = p_vals[2],
                       EigCentr = p_vals[3],
                       EdgeDen = p_vals[4],
                       NoModule = p_vals[5],
                       NetModularity = p_vals[6]
                       # Add more columns as needed
                     ))

# Print p-values
cat(paste0("Time period: ", Time, "\n"))
cat("Observed vs Null p-values:\n")
print(p_vals)

# Save node-specific p-values for Pre and Post time periods
save_node_specific_pvals("Pre", metrics_Pre)
save_node_specific_pvals("Post", metrics_Post)

# Save the metrics and results
write.csv(Reg_metrix, file = 'Results/Reg_metrix_Null_P_14_06_24.csv')

