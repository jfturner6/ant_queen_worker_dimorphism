######## WRP and MF Regression Script #######
# Juliet January 2025
##  ALL SPECIES --- DESCRIBED AND INFERRED DATA ##
setwd("/drives/4tb/Juliet/since_Nov_24/regressions/MF_WRP")
.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))

library(tidyverse)
library(ape)
library(phytools)
library(ggplot2)
library(geiger)
library(MuMIn)
library(MCMCglmm)
#library(arm)
#library(lme4)
library(purrr)
library(mulTree)
#install.packages("remotes") #Packages needed to install mulTree
#remotes::install_github("TGuillerme/mulTree") #Packages needed to install mulTree
#Further advice for installing mulTree packages can be found on https://github.com/TGuillerme/mulTree


data <- read.csv("Nov24_data.csv", header=T) 
#trees from Economo https://datadryad.org/stash/dataset/doi:10.5061/dryad.g579t7k
NCuniform_stem_posterior <- read.tree("15k_NCuniform_stem_posterior.tre")
NCuniform_crown_posterior <- read.tree("15k_NCuniform_crown_posterior.tre")
FBD_stem_posterior <- read.tree("15K_FBD_stem_posterior.tre")
FBD_crown_posterior <- read.tree("15K_FBD_crown_posterior.tre")
all_400_trees <-c(NCuniform_stem_posterior, NCuniform_crown_posterior, FBD_stem_posterior, FBD_crown_posterior)

### QUEEN MATING FREQUENCY ###
## filter (I did not filter out gamergate species here)
## CONFIDENCE = >=0 meaning this analysis will run on all species 
data$queen.mating.frequency.original <- gsub(",", ".", data$queen.mating.frequency.original, fixed = T)
data$queen.mating.frequency.original<-as.numeric(data$queen.mating.frequency.original)
data <- dplyr::filter(data, type == 'ant', queen.mating.frequency.original >=0, Inherent_Constraints_Scale >=1, Constraints_Confidence >=0, Supercolonial.Helantera.review_no.stringent.nonstringent <1, social.parasite.Antwiki_no.yes < 1, clonal < 1, hybridisation < 1, Taxonomic_info.Antcat.registered_undescribed.species_Antcat.unregistered_incertae.sedis._subspecies == 0)
data <-dplyr::select(data, animal, Genus, type, queen.mating.frequency.original, Inherent_Constraints_Scale)
#data$log10.queen.mating.frequency.original<-log10(data$queen.mating.frequency.original)
data$animal <- gsub("_", ".", data$animal, fixed = T)
data<-filter(data, animal %in% all_400_trees[[1]]$tip.label)
data$Genus <- sapply(strsplit(data$Genus, "-"), "[[",1)
#View(data) # 60 species jan25

#Prune each of the 3 different trees
pruned_all_400_trees_WRP.MF<- drop.tip.multiPhylo(all_400_trees, setdiff(all_400_trees[[1]]$tip.label, data$animal)) 
#This is now a multiphylo object that has been pruned

#pruned_all_400_trees_WRP.MF[[1]] # 60 tips/species 


############################################################
#MCMCglmm with mulTree
############################################################
mulTree_data_WRP.MF <- as.mulTree(data = data, tree = pruned_all_400_trees_WRP.MF,
                                          taxa = "animal")
multree_formula_WRP.MF <- queen.mating.frequency.original ~ Inherent_Constraints_Scale  
nitt <-11000000 
burnin <- 1000000 
thin <- 5000 # parameters based on caste number methods
mul_parameters <- c(nitt, thin, burnin)
set.seed(123)

prior_inv_wishart <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002)
)

mulTree(mulTree.data = mulTree_data_WRP.MF, formula = multree_formula_WRP.MF, priors = prior_inv_wishart,
        parameters = mul_parameters, output = "400_MF_WRP_exp", ESS = 1000, convergence = 1.1,
        chains = 2, family = "exponential") 


