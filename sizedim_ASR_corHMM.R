#Juliet January 2025

# Ancestral state reconstruction of queen-worker size dimorphism


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
registerDoParallel(4)


#Load custom functions for running corHMM models in parallel with adjustable model settings
corHMM_multiphylo <- function(multiphylo, trait_data, model_type, rate_cat) {
  
  # Run corHMM ancestral state reconstruction for each phylogeny
  foreach(i = 1:length(multiphylo)) %dopar% {
    output <- corHMM(phy = multiphylo[[i]], data = trait_data, rate.cat = rate_cat, model = model_type)
    
    # Define dynamic output file name
    output_file <- paste0("corHMM_ASR_sizedim_only_", rate_cat, "rat_", model_type, "_", i, ".rda")
    
    # Save output to a file
    saveRDS(output, file = output_file)
    cat("Saved output to", output_file, "\n")
  }
}

#########

#Read in data file
ASRdata <- read.csv("/drives/4tb/Juliet/since_Nov_24/SD_ASRs/Sizedim_ASR/ASR_sizedim_df.csv")

ASRdata <- ASRdata %>%
  rename(SD = dimorphism.binary) #1 = HIGH, 0 = LOW

# note: may need to assign WRP as numeric

#Read in sample of 400 phylogenetic trees
ant_trees_pruned <- read.tree(file ="/drives/4tb/Juliet/since_Nov_24/Sizedim_ASR/ASR_sizedim_tree.tre")

#data is already pruned to tree and filtered for data

### ESTIMATING ANCESTRAL STATES ###
setwd("/drives/4tb/Juliet/since_Nov_24/Sizedim_ASR/sizedim_output") #Assign file path for models to be saved
corHMM_multiphylo(multiphylo = ant_trees_pruned, trait_data = ASRdata, model_type = "ER", rate_cat = 1)
