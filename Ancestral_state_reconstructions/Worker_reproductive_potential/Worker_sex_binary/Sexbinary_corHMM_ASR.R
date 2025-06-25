#### SEX BINARY ASR ####
# Juliet 16th January 2025


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
    output_file <- paste0("corHMM_ASR_sexbinary_", rate_cat, "rat_", model_type, "_", i, ".rda")
    
    # Save output to a file
    saveRDS(output, file = output_file)
    cat("Saved output to", output_file, "\n")
  }
}

#########

#Read in data file
ASRdata <- read.csv("ASR_WRP_df.csv") #389 species

ASRdata <- ASRdata %>%
  mutate(sex_binary = recode(Inherent_Constraints_Scale,
                                     `1` = "sex",
                                     `2` = "sex",
                                     `3` = "loss",
                                     `4` = "loss")) 
#code as binary 0/1
ASRdata$sex_binary <- ifelse(ASRdata$sex_binary == "sex", "0", "1") #X1 SEX, X2 NO SEX
ASRdata$sex_binary <- as.factor(ASRdata$sex_binary)
ASRdata2 <-dplyr::select(ASRdata, animal, sex_binary)



# note: may need to assign WRP as numeric

#Read in sample of 400 phylogenetic trees
ant_trees_pruned <- read.tree(file ="ASR_WRP_pruned_tree.tre")

#data is already pruned to tree and filtered for WRP data

### ESTIMATING ANCESTRAL STATES ###
setwd("/drives/4tb/Juliet/since_Nov_24/WRP_ASRs/Sex_binary_ASR/rat2_ARD") #Assign file path for models to be saved
corHMM_multiphylo(multiphylo = ant_trees_pruned, trait_data = ASRdata2, model_type = "ARD", rate_cat = 2)

