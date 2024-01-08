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
Pre <- read.csv("Data/Speciesbytraitmatrix_Pre.csv")
Post <- read.csv("Data/Speciesbytraitmatrix_Post.csv")
Early <- read.csv("Data/Speciesbytraitmatrix_Early.csv")
Late <- read.csv("Data/Speciesbytraitmatrix_Late.csv")

#feeding_deposit is missing from: Post, Early and Late
Post$feeding_deposit <- rep(0, 37)
Early$feeding_deposit <- rep(0, 62)
Late$feeding_deposit <- rep(0, 79)

#Combine into one
All <- combine(Pre, Post, Early, Late)
#rename and reorder
All <- All %>% relocate(source, .before = X)
colnames(All) <- c("Time", "taxon", "F_D", "F_G", "F_M", "F_P", "F_S", "M_Fac",
                   "M_fast", "M_non", "M_slow", "S_G", "S_L", "S_S", "S_M", "S_T", "S_VL",
                   "T_E", "T_I", "T_P")
##### Make a counter file for times for looping:  ----
countR <- data.frame(Time = c("Pre", "Post", "Early", "Late" ))


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
  t_Degree = degree(TGraph, v = V(TGraph), mode = c("all"), loops = F, 
                      normalized = T) # no of connections of a trait, normalised for comparison
  t_DegCentr = centr_clo(TGraph, mode = "all")$res  # trait level centrality from No of connections
  tModule = data.frame(WT_partition$membership); colnames(tModule) = c("tModule")
  tModule$Trait = V(TGraph)$name
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
write.csv(Reg_metrix, 'Results/Reg_metrix_16_12_22.csv')

