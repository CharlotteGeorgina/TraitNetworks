#Packages 
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


#### Data #### 
CWM <- read.csv("Data/CWM_JP_Fish.csv")
CWM <- CWM %>% rename(Collection_name = X)

#organise and rename
CWM <- CWM %>% relocate(Cluster, .before = Collection_name)
CWM <- CWM %>% column_to_rownames(., var = 'Collection_name')
colnames(CWM) <- c("Cluster", "ML_Large", "ML_Medium", "ML_Small",
                   "ML_V.Small", "PLD_Long", "PLD_Medium", "PLD_Short",
                   "PLD_V.Long", "T_Corrallivore", "T_Detritivore",
                   "T_Herbivore", "T_Omnivore", "T_Piscivore",
                   "T_Planktivore", "T_Predator", "P_Benthic",
                   "P_CnidarianA", "P_Demersal", "P_EchinodermA",
                   "P_Pelagic", "P_ReefPelagic", "P_SandA", "P_SubBenthic",
                   "P_UpperBenthic", "R_Brooders", "P_Demersal",
                   "P_livebearers", "P_Nesters", "P_Scatterers")

##### Make a counter file for regions for looping:  ----
countR <- data.frame(Cluster = c("Tropical", "Warm Subtropical", "Cold Subtropical", 
                             "Temperate"))


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
  out_prob <- CWM_prob; csvfile<-paste0("Results/",Era,"/outprob_",Era,".csv"); write.csv(out_prob, csvfile, row.names = F) 
  
  #binarise the prob table, < 0.05 = 1, >0.05 = 0
  CWM_prob <- ifelse(CWM_prob < 0.05, 1, 0)
  out_prob_s <- CWM_prob; csvfile<-paste0("Results/",Era,"/outprob_s",Era,".csv"); write.csv(out_prob_s, csvfile, row.names = F)
  
  #Get correlation matrix
  CWM_eff <- CWM_Cor$r
  effectM = as.matrix(CWM_eff); csvfile<-paste0("Results/",Era,"/outeff_",Era,".csv"); write.csv(effectM, csvfile, row.names = F) #Get effects table
  
  ### Generate a matrix where links are significant
  SignMatrix <- out_prob_s * effectM
  csvfile=paste0("Results/",Era,"/SignMatrix_",Era,".csv"); write.csv(SignMatrix, csvfile, row.names=F)
  
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
  csvfile=paste0("Results/",Era,"/t_metrix_2_",Era,".csv"); write.csv(t_metrix, csvfile, row.names = F)
  
  ### Network-wide metrics
  DegCents = centr_clo(TGraph, mode = "all")$centralization
  BetwCentr = centr_betw(TGraph, directed = FALSE)$centralization
  EigCentr = centr_eigen(TGraph, directed = FALSE)$centralization
  EdgeDen = edge_density(TGraph, loops = FALSE)
  NoModule = length(WT_partition)
  NetModularity = modularity(TGraph,membership(WT_partition))
  
  # Region level summary table (DegCents, BetwCentr, EigCentr, NoModule, NetModularity) 
  Reg_metrix[nrow(Reg_metrix)+1,] = c(Era, DegCents, BetwCentr, EigCentr, EdgeDen, NoModule, NetModularity)
  
  #Plot and save
  V(TGraph)$color <- ifelse(V(TGraph)$name == "MaxLength_Large", "dodgerblue4", ifelse(V(TGraph)$name == "MaxLength_Medium", "dodgerblue3", 
                                                                           ifelse(V(TGraph)$name == "MaxLegnth_Small", "dodgerblue1", 
                                                                                  ifelse(V(TGraph)$name == "MaxLength_V.Small", "deepskyblue2", 
                                                                                         ifelse(V(TGraph)$name == "PLD_long", "aquamarine3",
                                                                                                ifelse(V(TGraph)$name == "PLD_Medium", "aquamarine4",
                                                                                                       ifelse(V(TGraph)$name == "PLD_Short", "aquamarine2",
                                                                                                              ifelse(V(TGraph)$name == "PLD_V.Long", "darkslategray3",
                                                                                                                     ifelse(V(TGraph)$name == "Trophic_Corrallivore", "magenta4",
                                                                                                                            ifelse(V(TGraph)$name == "Trophic_Detritivore", "magenta3",
                                                                                                                                   ifelse(V(TGraph)$name == "Trophic_Herbivore", "mediumpurple",
                                                                                                                                          ifelse(V(TGraph)$name == "Trophic_Omnivore", "mediumorchid",
                                                                                                                                                 ifelse(V(TGraph)$name == "Trophic_Piscivore", "honeydew4",
                                                                                                                                                        ifelse(V(TGraph)$name == "Trophic_Planktivore", "honeydew2",
                                                                                                                                                               ifelse(V(TGraph)$name == "Trophic_Predator", "mistyrose2",
                                                                                                                                                                      ifelse(V(TGraph)$name == "Position_Benthic", "mistyrose",
                                                                                                                                                                             ifelse(V(TGraph)$name =="Position_CnidarianAssociated", "darkorange3",
                                                                                                                                                                                    ifelse(V(TGraph)$name =="Position_Demersal", "darkorange",
                                                                                                                                                                                           ifelse(V(TGraph)$name == "B_intra", "gold2",
                                                                                                                                                                                                  ifelse(V(TGraph)$name == "B_none", "darkgoldenrod1","red"))))))))))))))))))))
  

  
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
  pdf(file=paste0("Results/Graphs/",Era,".pdf"))
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
write.csv(Reg_metrix, 'Results/Eocene_Oligocene_Reg_metrix_22_03_23.csv')