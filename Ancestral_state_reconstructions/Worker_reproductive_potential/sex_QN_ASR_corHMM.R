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
registerDoParallel(25)


# reconstruction for WRP only, using WRP subset and pruned tree. need to convert to binary sex/no sex here. will incorporate queen number in MCMCglmm step of analysis. 07/12/24


#Load custom functions for running corHMM models in parallel with adjustable model settings
corHMM_multiphylo <- function(multiphylo, trait_data, model_type, rate_cat) {
  
  # Run corHMM ancestral state reconstruction for each phylogeny
  foreach(i = 1:length(multiphylo)) %dopar% {
    output <- corHMM(phy = multiphylo[[i]], data = trait_data, rate.cat = rate_cat, model = model_type)
    
    # Define dynamic output file name
    output_file <- paste0("sex_QN_corHMM_ASR", rate_cat, "rat_", model_type, "_", i, ".rda")
    
    # Save output to a file
    saveRDS(output, file = output_file)
    cat("Saved output to", output_file, "\n")
  }
}

#########

#Read in data file
ASRdata <- read.csv("/drives/4tb/Juliet/since_Nov_24/Sex_QN_ASR/ASR_QN_WRP_df.csv")

ASRdata <- ASRdata %>%
  mutate(totipotency_binary = recode(Inherent_Constraints_Scale,
                                     `1` = "toti",
                                     `2` = "toti",
                                     `3` = "loss",
                                     `4` = "loss")) # now defining by loss of totipotency, whether capable of sexual reproduction
#code as binary 0/1
ASRdata$totipotency_binary <- ifelse(ASRdata$totipotency_binary == "toti", "0", "1") # 1 = no sexual capacity (SWAPPED 12/09/2024)
ASRdata$totipotency_binary <- as.factor(ASRdata$totipotency_binary)
ASRdata <-dplyr::select(ASRdata, animal, totipotency_binary) 
 #98 species




#Read in sample of 400 phylogenetic trees
ant_trees_pruned <- read.tree(file ="/drives/4tb/Juliet/since_Nov_24/Sex_QN_ASR/ASR_QN_WRP_pruned.tree.tre")
#98 tips
#data is already pruned to tree and filtered for WRP data

### ESTIMATING ANCESTRAL STATES ###
setwd("/drives/4tb/Juliet/since_Nov_24/Sex_QN_ASR/rat2_ARD") #Assign file path for models to be saved
corHMM_multiphylo(multiphylo = ant_trees_pruned, trait_data = ASRdata, model_type = "ARD", rate_cat = 2)

