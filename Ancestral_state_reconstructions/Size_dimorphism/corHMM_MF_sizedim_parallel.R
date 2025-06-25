######################## ASR with corHMM - analysis across 400 trees######################
####Runs ASR for size dimorphism as a binary variable 
###Analysing all species for which data is available for sizedim as well as MF
#Runs models in parallel
#Juliet Oct 2024, adapted from Louis Bell-Roberts

.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))
setwd("/drives/4tb/Juliet/since_Nov_24/SD_ASRs/Sizedim_MF/corHMM_sizedim_MF")


#Load packages
library(dplyr)
library(ape)
library(phytools)
library(ggplot2)
library(nloptr,
        lib.loc = "/drives/4tb/Juliet/miniconda3/envs/my_r/lib/R/library/")
library(gridExtra,
        lib.loc = "/drives/4tb/Juliet/miniconda3/envs/my_r/lib/R/library/")
library(corHMM)
library(doParallel)

######### NEW BEST MODEL IS 1rat_ER

#Run in parallel
registerDoParallel(20) #use 25 and only one at a time or you will scare Ming ## reduced 2ndOct

#Load custom functions for running corHMM models
corHMM_multiphylo_1rat_ER <- function(multiphylo, trait_data) {
  
  # Run corHMM ancestral state reconstruction for each phylogeny
  foreach(i = 1:length(multiphylo)) %dopar% {
    output <- corHMM(phy = multiphylo[[i]], data = trait_data, rate.cat = 1, model = "ER")
    
    # Save output to a file
    saveRDS(output, file = paste0("corHMM_ASR_MF_sizedim_1rat_ER_", i, ".rda"))
    cat("Saved output to", paste0("corHMM_ASR_MF_sizedim_1rat_ER_", i, ".rda"), "\n")
  }
}

#########

#Read in data file
ant_data <- read.csv("MF_sizedim_ASR_df.csv") #dimorphism is in label format

#Set variables so that they're in the correct structure
#ant_data[ant_data == ""] <- NA #Replace blank by NA

#Read in sample of 400 phylogenetic trees
ant.trees <- read.tree(file ="MF_sizedim_ASR_tree.tre")

### ESTIMATING ANCESTRAL STATES ###
## Models to estimate queen-worker size dimorphism at each node on the phylogeny ##

setwd("/drives/4tb/Juliet/since_Nov_24/SD_ASRs/Sizedim_MF/corHMM_sizedim_MF/sizedim_MF_output") #Assign file path for models to be saved
corHMM_multiphylo_1rat_ER(ant.trees, ant_data)


########END OF SCRIPT##########