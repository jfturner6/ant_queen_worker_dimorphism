#DIAGNOSTICS
.libPaths(c(.libPaths(), "/drives/4tb/modules/R")) 
library(tidyverse)
library(MCMCglmm)
library(ggplot2)
library(mulTree)

one_model_WRP_CS <- read.mulTree("400_CS_WRP-tree1_chain1", model = TRUE)
two_model_WRP_CS <- read.mulTree("400_CS_WRP-tree1_chain2", model = TRUE)

summary(one_model_WRP_CS) 
summary(two_model_WRP_CS) 



all_models_converged <-read.mulTree("400_CS_WRP", convergence = TRUE)

################################ DIAGNOSTICS ###################################

# https://www.statlect.com/fundamentals-of-statistics/Markov-Chain-Monte-Carlo-diagnostics

##################


#Custom functions


#Calculates SRF (Gelman-Rubin) estimate
##Input the convergence models from read.mulTree with convergence set to TRUE
# srf_function <- function(multree_conv) {
#   
#   # Initialize an empty vector to store the SRFs
#   list_psrf <- c()
#   
#   # Loop through the list of gelman.diag objects
#   for (i in 1:length(multree_conv)) {
#     # Get the first column of the ith object and calculate its range
#     ith_psrf <- multree_conv[[i]]$psrf[,2]
#     
#     list_psrf <- c(list_psrf, ith_psrf)
#   }
#   
#   #PSR statistics
#   return(quantile(list_psrf))
#   # View(as.data.frame(list_psrf_CS_Caste))
#   return(length(list_psrf))
# }

srf_function <- function(multree_conv) {
  
  # Initialize an empty vector to store the SRFs
  list_psrf <- c()
  
  # Loop through the list of gelman.diag objects
  for (i in 1:length(multree_conv)) {
    # Get the first column of the ith object and calculate its range
    ith_psrf <- multree_conv[[i]]$psrf
    
    list_psrf <- c(list_psrf, ith_psrf)
  }
  #PSR statistics
  return(list_psrf)
}



#Write a function that reads in all of the individual models
model_reader_function <- function() {
  
  ## get list of all rda files in the working directory
  model_files <- list.files(pattern = "chain")
  ## read in all models
  # Create an empty list to store the loaded data
  mcmc_list <- list()
  
  #Read in each MCMCglmm model
  mcmc_list <- lapply(model_files, function(file) {
    read.mulTree(file, model = TRUE)
  }) #mcmc_list is the list of all of the models read in
  return(mcmc_list)
}



#Write a function that tests levels of autocorrelation for the fixed effects for each model - provides the values in the second row of the autocorr.diag function. These should be <0.1
auto_cor_function_fix <- function(mcmc_list){
  
  # Apply the autocorr.diag function to each element of the list
  autocorr_list <- lapply(mcmc_list, function(x) {
    autocorr.diag(x$Sol)
  })
  
  # Use lapply to extract the second row of each matrix
  second_row_list <- lapply(autocorr_list, function(x) {
    x[2,]
  })
  
  # Convert the list to a vector
  second_row_vec <- unlist(second_row_list)
  return(second_row_vec)
}


#Write a function that tests levels of autocorrelation for the random effects for each model - provides the values in the second row of the autocorr.diag function. These should be <0.1
auto_cor_function_rand <- function(mcmc_list){
  
  # Apply the autocorr.diag function to each element of the list
  autocorr_list <- lapply(mcmc_list, function(x) {
    autocorr.diag(x$VCV)
  })
  
  # Use lapply to extract the second row of each matrix
  second_row_list <- lapply(autocorr_list, function(x) {
    x[2,]
  })
  
  # Convert the list to a vector
  second_row_vec <- unlist(second_row_list)
  return(second_row_vec)
}



#Two functions: one that calculates the ESS values for the fixed effects for each model across the full sample of models and one that calculates the ESS values for the random effects for each model across the full sample of models

##Fixed
ESS_function_fix <- function(mcmc_list) {
  
  ESS_list <- lapply(mcmc_list, function(x) {
    effectiveSize(x$Sol)
  })
  
  # Convert the list to a vector
  ESS_vec <- unlist(ESS_list)
  
  return(ESS_vec)
}
##Random
ESS_function_rand <- function(mcmc_list) {
  
  ESS_list <- lapply(mcmc_list, function(x) {
    effectiveSize(x$VCV)
  })
  
  # Convert the list to a vector
  ESS_vec <- unlist(ESS_list)
  
  return(ESS_vec)
}

##################


### CONVERGENCE ###
#Calculate whether scale-reduction factor (Gelman-Rubin) is <1.1 for all 400 models
View(as.data.frame(srf_function(all_models_converged)))     

### AUTOCORRELATION ###
#Test levels of autocorrelation for the fixed and random effects for each model
#Read in all of the individual models
#get list of all rda files in the working directory
model_files_CS_IRCS <- list.files(pattern = "chain")
#read in all models
#Create an empty list to store the loaded data
mcmc_list_CS_IRCS <- list()
#Read in each MCMCglmm model
mcmc_list_CS_IRCS <- lapply(model_files_CS_IRCS, function(file) {
  read.mulTree(file, model = TRUE)})


## Fixed effects ##
#Apply the autocorr.diag function to each element of the list
autocorr_list_CS_IRCS_fix <- lapply(mcmc_list_CS_IRCS, function(x) {
  autocorr.diag(x$Sol)
})
#Use lapply to extract the second row of each matrix
second_row_list_fix <- lapply(autocorr_list_CS_IRCS_fix, function(x) {
  x[2,]
})
#Convert the list to a vector
second_row_vec_fix <- unlist(second_row_list_fix)
quantile(second_row_vec_fix)
View(as.data.frame(second_row_vec_fix))
#  should not exceed 0.1 or be less than -0.1      

## Random effects ##
#Apply the autocorr.diag function to each element of the list
autocorr_list_CS_IRCS_rand <- lapply(mcmc_list_CS_IRCS, function(x) {
  autocorr.diag(x$VCV)
})

#lapply to extract the second row of each matrix
second_row_list_rand <- lapply(autocorr_list_CS_IRCS_rand, function(x) {
  x[2,]
})
#Convert the list to a vector
second_row_vec_rand <- unlist(second_row_list_rand)
quantile(second_row_vec_rand)
# 
View(as.data.frame(second_row_vec_rand)) 

# Ideally, levels of autocorrelation should not exceed 0.1 or be less than -0.1. 
# But it's OK if it's only slightly over, or only a small proportion of models 
# have a value that is only slighlty greater than 0.1 


### EFFECTIVE SAMPLE SIZE ###
#Calculate effective sample size for both the fixed and random effects for each model
#https://www.johndcook.com/blog/2017/06/27/effective-sample-size-for-mcmc/

## Fixed effects ##
#Test levels of autocorrelation for the fixed effects for each model 
#Apply the autocorr.diag function to each element of the list
ESS_list_CS_IRCS_fix <- lapply(mcmc_list_CS_IRCS, function(x) {
  effectiveSize(x$Sol)
})
#Convert the list to a vector
ESS_vec_CS_IRCS_fix <- unlist(ESS_list_CS_IRCS_fix)
quantile(ESS_vec_CS_IRCS_fix) # should be 1000.    --------- 2004

## Random effects ##
#Test levels of autocorrelation for the fixed effects for each model 
#Apply the autocorr.diag function to each element of the list
ESS_list_CS_IRCS_rand <- lapply(mcmc_list_CS_IRCS, function(x) {
  effectiveSize(x$VCV)
})
#Convert the list to a vector
ESS_vec_CS_IRCS_rand <- unlist(ESS_list_CS_IRCS_rand)
quantile(ESS_vec_CS_IRCS_rand)  
# should be 1000 ------ 2433

#Summary:
# Convergence = 
# Autocorrelation = 
# ESS = 



