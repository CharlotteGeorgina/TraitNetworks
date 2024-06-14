# Load necessary libraries
library(corrplot)
library(tidyverse)  # Includes dplyr and other packages
library(cooccur)
library(circlize)
library(igraph)
library(Hmisc)

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


##### loop through Regions, making eveything ----

Reg_metrix = data.frame(Cluster=as.character(), DegCents = numeric(), BetwCentr = numeric(),
                        EigCentr  = numeric(), EdgeDen = numeric(), NoModule = numeric(), NetModularity = numeric(), stringsAsFactors=FALSE)

for (i in 1:nrow(countR)) { #circle through the Regions
  # work out region
  Cluster = countR[i,1]
  
  # subsetting by Region
  CWM_i = CWM[CWM$Cluster == Cluster,]  # select the records for the Region we want
  CWM_i = CWM_i[,-1]
  CWM_Cor <- rcorr(as.matrix(CWM_i)) #produces correlation matrix 

  ## probtable: write results for probability of pairwise co-occurrence to output table
  CWM_prob <- CWM_Cor$P
  out_prob <- CWM_prob; csvfile<-paste0("Results/",Cluster,"/outprob_",Cluster,".csv"); write.csv(out_prob, csvfile, row.names = F) 
  
  #binarise the prob table, < 0.05 = 1, >0.05 = 0
  CWM_prob <- ifelse(CWM_prob < 0.05, 1, 0)
  out_prob_s <- CWM_prob; csvfile<-paste0("Results/",Cluster,"/outprob_s",Cluster,".csv"); write.csv(out_prob_s, csvfile, row.names = F)
  
  #Get correlation matrix
  CWM_eff <- CWM_Cor$r
  effectM = as.matrix(CWM_eff); csvfile<-paste0("Results/",Cluster,"/outeff_",Cluster,".csv"); write.csv(effectM, csvfile, row.names = F) #Get effects table
  
  ### Generate a matrix where links are significant
  SignMatrix <- out_prob_s * effectM
  csvfile=paste0("Results/",Cluster,"/SignMatrix_",Cluster,".csv"); write.csv(SignMatrix, csvfile, row.names=F)
  
  #change na and Nan to 0's 
  SignMatrix[is.na(SignMatrix)] <- 0
  SignMatrix <- ifelse(SignMatrix>0.1, 1, 0)    # Binarize 
  
  #cut out the zero columns 
  SignMatrix2 = SignMatrix[rowSums(SignMatrix == 0) != ncol(SignMatrix), colSums(SignMatrix == 0) != nrow(SignMatrix)]
  #csvfile=paste0("SignMatrix2_",HemiRegion,".csv"); write.csv(SignMatrix2, csvfile)
  
  ### graph-theoretic metrix code ### calls 'igraph'  
  TGraph = graph_from_adjacency_matrix(SignMatrix2, mode = "undirected")
  plot(TGraph) #it will ignore the weights at this point
  
  ### cluster any substructures (most techniques dont work because of neg relationships)
  WT_partition = cluster_walktrap(TGraph, weights = NULL)
  
  ### Node/ Trait specific metrics:
  t_Degree = igraph::degree(TGraph, v = V(TGraph), mode = c("all"), loops = F, 
                      normalized = T) # no of connections of a Trait, normalised for comparison
  t_DegCentr = centr_clo(TGraph, mode = "all")$res  # Trait level centrality from No of connections
  tModule = data.frame(WT_partition$membership); colnames(tModule) = c("tModule")
  tModule$Trait = V(TGraph)$name
  tModule <- rowid_to_column(tModule, "ID")
  
  # Trait level modularity
  t_metrix = cbind(Trait = vertex_attr(TGraph, "name"),t_Degree,t_DegCentr,tModule)
  t_metrix[is.na(t_metrix)] <- 0
  csvfile=paste0("Results/",Cluster,"/t_metrix_2_",Cluster,".csv"); write.csv(t_metrix, csvfile, row.names = F)
  
  ### Network-wide metrics
  DegCents = centr_clo(TGraph, mode = "all")$centralization
  BetwCentr = centr_betw(TGraph, directed = FALSE)$centralization
  EigCentr = centr_eigen(TGraph, directed = FALSE)$centralization
  EdgeDen = edge_density(TGraph, loops = FALSE)
  NoModule = length(WT_partition)
  NetModularity = modularity(TGraph,membership(WT_partition))
  
  # Region level summary table (DegCents, BetwCentr, EigCentr, NoModule, NetModularity) 
  Reg_metrix[nrow(Reg_metrix)+1,] = c(Cluster, DegCents, BetwCentr, EigCentr, EdgeDen, NoModule, NetModularity)
  
  #Plot and save
  V(TGraph)$color <- ifelse(V(TGraph)$name == "ML_Large", "dodgerblue4", ifelse(V(TGraph)$name == "ML_Medium", "dodgerblue3", 
                                                                           ifelse(V(TGraph)$name == "ML_Small", "dodgerblue1", 
                                                                                  ifelse(V(TGraph)$name == "ML_V.Small", "deepskyblue2", 
                                                                                         ifelse(V(TGraph)$name == "PLD_Long", "aquamarine3",
                                                                                                ifelse(V(TGraph)$name == "PLD_Medium", "aquamarine4",
                                                                                                       ifelse(V(TGraph)$name == "PLD_Short", "aquamarine2",
                                                                                                              ifelse(V(TGraph)$name == "PLD_V.Long", "darkslategray3",
                                                                                                                     ifelse(V(TGraph)$name == "T_Corrallivore", "magenta4",
                                                                                                                            ifelse(V(TGraph)$name == "T_Detritivore", "magenta3",
                                                                                                                                   ifelse(V(TGraph)$name == "T_Herbivore", "mediumpurple",
                                                                                                                                          ifelse(V(TGraph)$name == "T_Omnivore", "mediumorchid",
                                                                                                                                                 ifelse(V(TGraph)$name == "T_Piscivore", "honeydew4",
                                                                                                                                                        ifelse(V(TGraph)$name == "T_Planktivore", "honeydew2",
                                                                                                                                                               ifelse(V(TGraph)$name == "T_Predator", "mistyrose2",
                                                                                                                                                                      ifelse(V(TGraph)$name == "P_Benthic", "mistyrose",
                                                                                                                                                                             ifelse(V(TGraph)$name =="P_AnthozoanA", "darkorange3",
                                                                                                                                                                                    ifelse(V(TGraph)$name =="P_Demersal", "darkorange",
                                                                                                                                                                                           ifelse(V(TGraph)$name == "P_Pelagic", "gold2",
                                                                                                                                                                                                  ifelse(V(TGraph)$name == "P_ReefPelagic", "goldenrod2",
                                                                                                                                                                                                         ifelse(V(TGraph)$name == "P_SandA", "darkgoldenrod3",
                                                                                                                                                                                                               ifelse(V(TGraph)$name == "P_SubBenthic", "darkgoldenrod4",
                                                                                                                                                                                                                     ifelse(V(TGraph)$name == "P_UpperBenthic", "gold4",
                                                                                                                                                                                                                            ifelse(V(TGraph)$name == "R_Brooders", "coral3",
                                                                                                                                                                                                                                   ifelse(V(TGraph)$name == "R_livebearers", "coral",
                                                                                                                                                                                                                                          ifelse(V(TGraph)$name == "R_Nesters", "coral4",
                                                                                                                                                                                                                                                 ifelse(V(TGraph)$name == "R_Scatterers", "violetred2",
                                                                                                                                                                                                                                                        ifelse(V(TGraph)$name == "R_Demersal", "violetred4","red"))))))))))))))))))))))))))))
  

  
  #Assign group ideas = modules 
  group_ids <- lapply(tModule %>% split(.$tModule), function(grp) { grp$ID })
  group_ids
  # set colour
  group_colour <- c("gray", "gray", "gray", "gray", "gray")
  #same but add transparency
  group_colour_fill <- paste0(group_colour, '95')
  
  #set node sizes
  V(TGraph)$size<-t_metrix$t_Degree*100
  #specify layout
  layout = layout_with_kk(TGraph)
  
  ##to save plot 
  pdf(file=paste0("Results/Graphs/",Cluster,".pdf"))
  plot(TGraph, 
       edge.width=2,
       edge.color="black",
       vertex.label.color="black",
       Vertex.label.dist = 3,
       #vertex.label.degree = pi/4,
       mark.groups = group_ids,
       mark.col = group_colour_fill,
       mark.border = group_colour,
       layout = layout, hovermode='closest')
  dev.off()



}


#save the Cooccurs and Metrix from the entire loop
write.csv(Reg_metrix, 'Results/JPFish_TropTemp_Reg_metrix_13_06_24.csv')
