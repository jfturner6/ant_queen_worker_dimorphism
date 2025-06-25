#### NEW WITH FIXED TREE 

# RUN 2 RATE ARD MODEL


#Juliet October 2024
.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))

#Load packages
library(dplyr)
library(ape)
library(phytools)
library(ggplot2)
library(corHMM)
library(doParallel)
library(gtools)

#Run in parallel using 50 cores
registerDoParallel(8)


#Load custom functions for running corHMM models in parallel with adjustable model settings
corHMM_multiphylo <- function(multiphylo, trait_data, model_type, rate_cat) {
  
  # Run corHMM ancestral state reconstruction for each phylogeny
  foreach(i = 1:length(multiphylo)) %dopar% {
    output <- corHMM(phy = multiphylo[[i]], data = trait_data, rate.cat = rate_cat, model = model_type)
    
    # Define dynamic output file name
    output_file <- paste0("corHMM_ASR_WRP_only_discrete_", rate_cat, "rat_", model_type, "_", i, ".rda")
    
    # Save output to a file
    saveRDS(output, file = output_file)
    cat("Saved output to", output_file, "\n")
  }
}

#########

#Read in data file
ASRdata <- read.csv("/drives/4tb/Juliet/since_Nov_24/WRP_4level_ASR/ASR_WRP_df.csv")
# all 389 species after pruning

ASRdata <- ASRdata %>%
  rename(WRP = Inherent_Constraints_Scale)

# note: may need to assign WRP as numeric

#Read in sample of 400 phylogenetic trees
ant_trees_pruned <- read.tree(file ="/drives/4tb/Juliet/since_Nov_24/WRP_4level_ASR/ASR_WRP_pruned_tree.tre")

#data is already pruned to tree and filtered for WRP data

### ESTIMATING ANCESTRAL STATES ###
setwd("/drives/4tb/Juliet/since_Nov_24/WRP_4level_ASR/rat2_ARD") #Assign file path for models to be saved
corHMM_multiphylo(multiphylo = ant_trees_pruned, trait_data = ASRdata, model_type = "ARD", rate_cat = 2)
