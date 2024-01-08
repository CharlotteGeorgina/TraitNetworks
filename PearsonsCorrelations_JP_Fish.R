#Packages 
library(dplyr)
library(tidyr)
library(Hmisc)
library(corrplot)
library(tidyverse)
library(ggplot2)
library(sqldf)
library (cooccur)
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

for (i in 1:nrow(countR)) { #circle through the Regions, Aus first, then Jpn
  # work out hemi and Region
  Cluster = countR[1,1]
  
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
  #cut out the zero columns (but not if I want all spp in a FG or All)
  SignMatrix2 = SignMatrix[rowSums(SignMatrix == 0) != ncol(SignMatrix), colSums(SignMatrix == 0) != nrow(SignMatrix)]
  #csvfile=paste0("SignMatrix2_",HemiRegion,".csv"); write.csv(SignMatrix2, csvfile)
  
  ### graph-theoretic metrix code ### calls 'igraph'  
  TGraph = graph_from_adjacency_matrix(SignMatrix2, mode = "undirected")
  plot(TGraph) #it will ignore the weights at this point
  
  ### cluster any substructures (most techniques dont work because of neg relationships)
  WT_partition = cluster_walktrap(TGraph, weights = NULL)
  
  ### Node/ species specific metrics:
  t_Degree = igraph::degree(TGraph, v = V(TGraph), mode = c("all"), loops = F, 
                      normalized = T) # no of connections of a Spp, normalised for comparison
  t_DegCentr = centr_clo(TGraph, mode = "all")$res  # Spp level centrality from No of connections
  tModule = data.frame(WT_partition$membership); colnames(tModule) = c("tModule")
  tModule$Trait = V(TGraph)$name
  tModule <- rowid_to_column(tModule, "ID")
  
  # spp level modularity
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
  #V(TGraph)$color <- ifelse(V(TGraph)$name == "MaxLength_Large", "dodgerblue4", ifelse(V(TGraph)$name == "MaxLength_Medium", "dodgerblue3", 
   #                                                                        ifelse(V(TGraph)$name == "MaxLegnth_Small", "dodgerblue1", 
    #                                                                              ifelse(V(TGraph)$name == "MaxLength_V.Small", "deepskyblue2", 
     #                                                                                    ifelse(V(TGraph)$name == "PLD_long", "aquamarine3",
      #                                                                                          ifelse(V(TGraph)$name == "PLD_Medium", "aquamarine4",
       #                                                                                                ifelse(V(TGraph)$name == "PLD_Short", "aquamarine2",
        #                                                                                                      ifelse(V(TGraph)$name == "PLD_V.Long", "darkslategray3",
         #                                                                                                            ifelse(V(TGraph)$name == "Trophic_Corrallivore", "magenta4",
          #                                                                                                                  ifelse(V(TGraph)$name == "Trophic_Detritivore", "magenta3",
           #                                                                                                                        ifelse(V(TGraph)$name == "Trophic_Herbivore", "mediumpurple",
            #                                                                                                                              ifelse(V(TGraph)$name == "Trophic_Omnivore", "mediumorchid",
             #                                                                                                                                    ifelse(V(TGraph)$name == "Trophic_Piscivore", "honeydew4",
              #                                                                                                                                          ifelse(V(TGraph)$name == "Trophic_Planktivore", "honeydew2",
               #                                                                                                                                                ifelse(V(TGraph)$name == "Trophic_Predator", "mistyrose2",
                #                                                                                                                                                      ifelse(V(TGraph)$name == "Position_Benthic", "mistyrose",
                 #                                                                                                                                                            ifelse(V(TGraph)$name =="Position_CnidarianAssociated", "darkorange3",
                  #                                                                                                                                                                  ifelse(V(TGraph)$name =="Position_Demersal", "darkorange",
                   #                                                                                                                                                                        ifelse(V(TGraph)$name == "B_intra", "gold2",
                    #                                                                                                                                                                              ifelse(V(TGraph)$name == "B_none", "darkgoldenrod1","red"))))))))))))))))))))
  
  #TGraph = set_vertex_attr(TGraph, "Module", value = tModule$tModule) # add module membership
  
  #V(TGraph)$color <- ifelse(V(TGraph)$Module == "1", "deep sky blue", ifelse(V(TGraph)$Module =="2", "LightSkyBlue1", 
  #                                                           ifelse(V(TGraph)$Module == "3","dark turquoise",
  #                                                                  ifelse(V(TGraph)$Module == "4","dodger blue", 
  #                                                                         ifelse(V(TGraph)$Module == "5","cornflower blue", "red")))))
  group_ids <- lapply(tModule %>% split(.$tModule), function(grp) { grp$ID })
  group_ids
  # set colour
  group_colour <- c("gray", "gray", "gray", "gray", "gray")
  #same but add transparency
  group_colour_fill <- paste0(group_colour, '95')
  
  V(TGraph)$size<-t_metrix$t_Degree*100
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


#### working out code ####
####Pearson's correlation ####
####Tropical####
Trop_Cor <- cor(Trop) #produces correlation matrix 
Trop_V <- rcorr(as.matrix(Trop)) #produces p-values for the matrix
##flattens it - so you can see the significance 
flattened_T <- flattenCorrMatrix(Trop_V$r, Trop_V$P)  
#remove correlations where the p value is higher than 0.05
flattened_sig <- flattened_T %>%
  filter(p < 0.05)
### CHECK THE FLATTENED FOR AN IDEA OF A GOOD R VALUE 
Trop_Cor <- as.matrix(Trop_Cor) #turn to matrix
correl_Trop <- ifelse(Trop_Cor>0.4, 1, 0)    # Binarize 
diag(correl_Trop)<-0 #remove the diagonal 1's
corr_Trop <- graph_from_adjacency_matrix(correl_Trop, mode = "undirected") #turn to graphical
plot(corr_Trop) #plot
##to save plot 
png(file="Results/Tropical_co-occurrence_r0.4.png",
    width=1000, height=800)
plot(corr_Trop)
dev.off()

#### Warm subtropical ####
WS_Cor <- cor(WarmSub)
WS_V <- rcorr(as.matrix(WarmSub)) #produces p-values for the matrix
##flattens it - so you can see the significance 
flattened_WS <- flattenCorrMatrix(WS_V$r, WS_V$P)  
#remove correlations where the p value is higher than 0.05
flattened_sig_WS <- flattened_WS %>%
  filter(p < 0.05)
### CHECK THE FLATTENED FOR AN IDEA OF A GOOD R VALUE 
WS_Cor <- as.matrix(WS_Cor)
correl_WS <- ifelse(WS_Cor>0.4, 1, 0)    # Binarize 
diag(correl_WS)<-0
corr_WS <- graph_from_adjacency_matrix(correl_WS, mode = "undirected")
plot(corr_WS)

##to save plot 
png(file="Results/WarmSubtropical_co-occurrence_r0.4.png",
    width=1000, height=800)
plot(corr_WS)
dev.off()

#### Cold Subtropical ####
CS_Cor <- cor(ColdSub)
CS_V <- rcorr(as.matrix(ColdSub)) #produces p-values for the matrix
##flattens it - so you can see the significance 
flattened_CS <- flattenCorrMatrix(CS_V$r, CS_V$P)  
#remove correlations where the p value is higher than 0.05
flattened_sig_CS <- flattened_CS %>%
  filter(p < 0.05)
### CHECK THE FLATTENED FOR AN IDEA OF A GOOD R VALUE 
CS_Cor  <- as.matrix(CS_Cor )
correl_CS <- ifelse(CS_Cor >0.4, 1, 0)    # Binarize 
diag(correl_CS)<-0
corr_CS <- graph_from_adjacency_matrix(correl_CS, mode = "undirected")
plot(corr_CS)
##to save plot 
png(file="Results/ColdSubtropical_co-occurrence_R0.4.png",
    width=1000, height=800)
plot(corr_CS)
dev.off()

Temp_Cor <- cor(Temp)
Temp_V <- rcorr(as.matrix(Temp)) #produces p-values for the matrix
##flattens it - so you can see the significance 
flattened_Temp <- flattenCorrMatrix(Temp_V$r, Temp_V$P)  
#remove correlations where the p value is higher than 0.05
flattened_sig_Temp <- flattened_Temp %>%
  filter(p < 0.05)
### CHECK THE FLATTENED FOR AN IDEA OF A GOOD R VALUE 
Temp_Cor <- as.matrix(Temp_Cor)
correl_Temp <- ifelse(Temp_Cor>0.4, 1, 0)    # Binarize 
diag(correl_Temp)<-0
corr_Temp <- graph_from_adjacency_matrix(correl_Temp, mode = "undirected")
plot(corr_Temp)
##to save plot 
png(file="Results/Temperate_co-occurrence_r0.2.png",
    width=1000, height=800)
plot(corr_Temp)
dev.off()


#Centrality measures using igraph
Trop_deg <- degree(corr_Trop)
Trop_bet <- betweenness(corr_Trop)
Trop_clos <- closeness(corr_Trop)

Trop_eig <- eigen_centrality(corr_Trop)$vector

Trop_cent_df <- data.frame(Trop_deg, Trop_bet, Trop_clos, Trop_eig)

WT_partition = cluster_walktrap(corr_Trop, weights = NULL)
NetModularity = modularity(corr_Trop,membership(WT_partition))

### cluster any substructures (most techniques dont work because of neg relationships)
WT_partition = cluster_walktrap(corr_Trop, weights = NULL)

### Node/ trait specific metrics:
Trait_Degree = degree(corr_Trop) # no of connections of a Spp, normalised for comparison
Trait_DegCentr = centr_clo(corr_Trop, mode = "all")$res  # Spp level centrality from No of connections
TraitModule = data.frame(WT_partition$membership); colnames(TraitModule) = c("TraitModule")
# spp level modularity
Trait_metrix = cbind(Trait = vertex_attr(corr_Trop, "name"),Trait_Degree,Trait_DegCentr,TraitModule)

### Network-wide metrics
Region = "Tropical"
DegCents = centr_clo(corr_Trop, mode = "all")$centralization
BetwCentr = centr_betw(corr_Trop, directed = FALSE)$centralization
EigCentr = centr_eigen(corr_Trop, directed = FALSE)$centralization
NoModule = length(WT_partition)
NetModularity = modularity(corr_Trop,membership(WT_partition))
# Region level summary table (DegCents, BetwCentr, EigCentr, NoModule, NetModularity) 
Reg_metrix = data.frame(Reg=as.character(), DegCents = numeric(), BetwCentr = numeric(),
                        EigCentr  = numeric(), NoModule = numeric(), NetModularity = numeric(), stringsAsFactors=FALSE)
Reg_metrix[nrow(Reg_metrix)+1,] = c(Region, DegCents, BetwCentr, EigCentr, NoModule, NetModularity)


### plot and save
# work out the no of species in (significant) links per FG
Sign_sp = melt(SignMatrix); Sign_sp = Sign_sp[Sign_sp$value != 0,] # 
Spp_i = as.data.frame(unique((Sign_sp$Var2))); colnames(Spp_i) = c("Fish")
Spp_i$FG = All_fishFG$FG[match(Spp_i$Fish, All_fishFG$Fish)]
Spp_i$Spp_Degree = Spp_metrix$Spp_Degree[match(Spp_i$Fish, Spp_metrix$Fish)]
Spp_i$Spp_DegCentr = Spp_metrix$Spp_DegCentr[match(Spp_i$Fish, Spp_metrix$Fish)]
Spp_i$SppModule = Spp_metrix$SppModule[match(Spp_i$Fish, Spp_metrix$Fish)]
csvfile=paste0("SignSpp_",HemiRegion,".csv"); write.csv(Spp_i, csvfile)
Spp_per_FG = data.frame(table(Spp_i$FG)); colnames(Spp_per_FG) = c("FG", "SPinFG")
Spp_i = Spp_i[order(Spp_i$FG),]
Spp_i$Fish = as.character(Spp_i$Fish)
SpeciesSignificance[nrow(SpeciesSignificance) +1, ] = c(HemiRegion, nrow(Spp_i), 
                                                        nrow(Spp_i)/nrow(SignMatrix)*100, length((unique(Spp_i$FG))),
                                                        length((unique(Spp_i$FG)))/19*100)


