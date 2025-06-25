######## WRP and Colony Size Regression Script #######
# JANUARY 2025

##  This subset includes species with all confidence levels (inferred and described) ##
# Does colony size predict worker reproductive potential?
.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))
setwd("/drives/4tb/Juliet/since_Nov_24/regressions/CS_WRP")

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
data$colony.size <- gsub(",", ".", data$colony.size, fixed = T)
data$colony.size<- as.numeric(data$colony.size)

#### COLONY SIZE ####
## filter (I did not filter out gamergate species here)
## CONFIDENCE >= 0 meaning this analysis will run on all species 

data <- dplyr::filter(data, type == 'ant', colony.size >=0, Inherent_Constraints_Scale >=1, Constraints_Confidence >= 0, Supercolonial.Helantera.review_no.stringent.nonstringent <1, social.parasite.Antwiki_no.yes < 1, clonal < 1, hybridisation < 1, Taxonomic_info.Antcat.registered_undescribed.species_Antcat.unregistered_incertae.sedis._subspecies == 0)

data <-dplyr::select(data, animal, Genus, type, colony.size, Inherent_Constraints_Scale)
data$log10.colony.size<-log10(data$colony.size)
data$animal <- gsub("_", ".", data$animal, fixed = T)
data<-filter(data, animal %in% all_400_trees[[1]]$tip.label)
data$Genus <- sapply(strsplit(data$Genus, "-"), "[[",1)


pruned_all_400_trees<- drop.tip.multiPhylo(all_400_trees, setdiff(all_400_trees[[1]]$tip.label, data$animal)) 

# VISUALISING THE VARIANCE FOR PRIOR SELECITON (?)
#ggplot(data, aes(x = colony.size, y = Inherent_Constraints_Scale)) + geom_point() + geom_smooth(method = "lm")
#model <- lm(Inherent_Constraints_Scale ~ colony.size, data = data)
#data$residuals <- residuals(model)
#residuals_variance <- var(data$residuals)
#print(residuals_variance)

############################################################
#MCMCglmm with mulTree
############################################################

#Prepare the data for the mulTree function
mulTree_Reprod_data_CS_WRP <- as.mulTree(data = data, tree = pruned_all_400_trees,
                                          taxa = "animal")
multree_formula_CS_WRP <- log10.colony.size ~ Inherent_Constraints_Scale
nitt <- 11000000 
burnin <- 1000000 
thin <- 5000 # parameters based on caste number methods
mul_parameters <- c(nitt, thin, burnin)
set.seed(123)


prior_inv_wishart <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002)
)


mulTree(mulTree.data = mulTree_Reprod_data_CS_WRP, formula = multree_formula_CS_WRP, priors = prior_inv_wishart,
        parameters = mul_parameters, output = "400_CS_WRP", ESS = 1000, convergence = 1.1,
        chains = 2, family = "gaussian")

#one_model_WRP_CS <- read.mulTree("400_WRP_CS-tree1_chain1", model = TRUE)
#two_model_WRP_CS <- read.mulTree("400_WRP_CS-tree1_chain2", model = TRUE)

#summary(one_model_WRP_CS) # p = 0.172 pre-fisher
#summary(two_model_WRP_CS) # p = 0.008 ** pre-fisher





