######## QUEEN NUMBER ~ Reproductive INHERENT Constraints Scale Regression Script #######
# Juliet January 2025

##  ALL SPECIES ##

#trees from Economo https://datadryad.org/stash/dataset/doi:10.5061/dryad.g579t7k


.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))
setwd("/drives/4tb/Juliet/since_Nov_24/regressions/QN_WRP")

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


ant_data <- read.csv("Nov24_data.csv", header=T) 

#Tree files:
NCuniform_stem_posterior <- read.tree("15k_NCuniform_stem_posterior.tre")
NCuniform_crown_posterior <- read.tree("15k_NCuniform_crown_posterior.tre")
FBD_stem_posterior <- read.tree("15K_FBD_stem_posterior.tre")
FBD_crown_posterior <- read.tree("15K_FBD_crown_posterior.tre")
#combine 
all_400_trees <-c(NCuniform_stem_posterior, NCuniform_crown_posterior, FBD_stem_posterior, FBD_crown_posterior)


#### QUEEN NUMBER ####.     
## filter (I did not filter out gamergate species)
## CONFIDENCE = 0 meaning this analysis will run on all species wiht data #removed polygyny.clean.new >=0,
ant_data$queen.number.continuous <- gsub(",", ".", ant_data$queen.number.continuous, fixed = T)
ant_data$queen.number.continuous<-as.numeric(ant_data$queen.number.continuous)
ant_data <- dplyr::filter(ant_data, type == 'ant', queen.number.continuous >=0, Inherent_Constraints_Scale >=1, Constraints_Confidence >=0, Supercolonial.Helantera.review_no.stringent.nonstringent <1, social.parasite.Antwiki_no.yes < 1, clonal < 1, hybridisation < 1, Taxonomic_info.Antcat.registered_undescribed.species_Antcat.unregistered_incertae.sedis._subspecies == 0)
ant_data <-dplyr::select(ant_data, animal, Genus, type, queen.number.continuous, Inherent_Constraints_Scale)
#data$log10.queen.number.continuous<-log10(data$queen.number.continuous)
ant_data$animal <- gsub("_", ".", ant_data$animal, fixed = T)
ant_data<-filter(ant_data, animal %in% all_400_trees[[1]]$tip.label)
ant_data$Genus <- sapply(strsplit(ant_data$Genus, "-"), "[[",1)
View(ant_data) #
pruned_all_400_trees_WRP_QN<- drop.tip.multiPhylo(all_400_trees, setdiff(all_400_trees[[1]]$tip.label, ant_data$animal)) 
#pruned_all_400_trees_WRP_QN[[1]] # 97 tips/species with both QN data and described constraints
ant_data$queen.number.continuous<-as.numeric(ant_data$queen.number.continuous)



############################################################
#MCMCglmm with mulTree
############################################################

#Prepare the data for the mulTree function
mulTree_data_WRP_QN <- as.mulTree(data = ant_data, tree = pruned_all_400_trees_WRP_QN,
                                          taxa = "animal")

#Set up model parameters 
# The formula 
multree_formula_WRP_QN <- queen.number.continuous ~ Inherent_Constraints_Scale
# The MCMC parameters (iterations, thining, burnin)
nitt <- 1100000 
burnin <- 100000 
thin <- 5000 # parameters based on caste number methods
mul_parameters <- c(nitt, thin, burnin)

prior_inv_wishart <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002)
)

mulTree(mulTree.data = mulTree_data_WRP_QN, formula = multree_formula_WRP_QN, priors = prior_inv_wishart,
        parameters = mul_parameters, output = "400_QN_WRP_exp", ESS = 1000, convergence = 1.1,
        chains = 2, family = "exponential") 

# Reading only one specific model
#one_model_WRP_QN <- read.mulTree("400_WRP_QN_poisson-tree1_chain1", model = TRUE)
#two_model_WRP_QN <- read.mulTree("400_WRP_QN_poisson-tree2_chain1", model = TRUE)
#summary(one_model_WRP_QN_poisson) 
#summary(two_model_WRP_QN_poisson)





