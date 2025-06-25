########################Plot transitions of the corHMM ASR for WRP and calculate summary statistics for each transition type######################
#Juliet January 2025, adapted from Louis Bell-Roberts


.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))

#Load packages
library(tidyverse)
library(ape)
library(phytools)
library(ggplot2)
library(geiger)
library(corHMM)
library(phangorn)
library(gtools)
library(ggtree)
library(tidytree)

setwd("/drives/4tb/Juliet/since_Nov_24/SD_ASRs/Sizedim_ASR")

#Load custom functions

#Function that takes a matrix or data frame as input and returns a vector of the column indices corresponding to the maximum values in each row - replicates the behaviour of the max.col function
my_max_col <- function(x) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("Input must be a matrix or data frame")
  }
  apply(x, 1, which.max)
}

#The extract_frequency function is designed to extract frequencies from a list of data frames (table_list) based on a specified pattern found in the 'Var1' column of each data frame. Here's a summary of what the function does:
extract_frequency <- function(table_list, pattern) {
  lapply(table_list, function(df) {
    if (pattern %in% df$Var1) {
      subset_df <- df[grep(pattern, df$Var1), ]
      # Extract the 'Freq' column from the subsetted data frame
      freq <- subset_df$Freq
    } else {
      # If 'pattern' is not present, return 0
      freq <- 0
    }
    return(freq)
  })
}

###############################################################

#Read in data file
sxData <- read.csv("/drives/4tb/Juliet/since_Nov_24/SD_ASRs/Sizedim_ASR/ASR_sizedim_df.csv")
# 153 species with highest confidence data, only WRP info and species names

sxData <- sxData %>%
  rename(SD = dimorphism.binary)


ant_trees_pruned <- read.tree(file ="/drives/4tb/Juliet/since_Nov_24/SD_ASRs/Sizedim_ASR/ASR_sizedim_tree.tre")

###############
###Read in the 400 corHMM model results in ascending order
###############

# create an empty list to hold the results
ASR_models <- list()

# set folder path
folder_path <- ("/drives/4tb/Juliet/since_Nov_24/SD_ASRs/Sizedim_ASR/sizedim_output")

# get the file names
file_names <- list.files(folder_path)

# sort the file names based on the corresponding numbers
sorted_file_names <- mixedsort(file_names)

# read in the files and save them to a list
ASR_models <- lapply(sorted_file_names, function(x) readRDS(file.path(folder_path, x)))

###############
####################################################################################################
############ NEW CODE 2ND APRIL 25 -- ROOT STATE PROBABILTIES
# Check first few models
#lapply(ASR_models[1:5], function(model) str(model$root.p)) #would work if not using 'yang' method

head(ASR_models[[1]]$states) #two states, 1 rate cat = (1,R1)    (2,R1)
############# NO NEED TO COLLAPSE INTO 1 RATE CAT AS ONLY HAD ONE CAT FOR SD ASR

#states were coded in as 0 = LOW, 1 = HIGH in corHMM 

# OVERALL ROOT PROBABILITIES
# Create a matrix for the root state probabilities (first row = root) for each of 400 models
root_state_matrix <- do.call(rbind, lapply(ASR_models, function(model) model$states[1, ]))
# Print the matrix
View(root_state_matrix)  # probabilities for each state for every model. 400 rows
# Calculate the average probability for each state  across all models.
overall_root_probs <- colMeans(root_state_matrix)
overall_root_probs
#   (1,R1)     (2,R1) 
# 0.4640248  0.5359752 --- state 2 is more likely overall (HIGH)

# CONFIDENCE INTERVALS
# Create a list for the root state probabilities for each of 400 models
root_state_list <- lapply(ASR_models, function(model) model$states[1, ])
# Print the first few elements of the list for verification
head(root_state_list)
# Combine all root probabilities for all 400 models
combined_probs <- do.call(rbind, root_state_list)
# Calculate the 95% credible interval for each state 
credible_intervals <- apply(combined_probs, 2, function(x) {
  # Calculate the 2.5th and 97.5th percentiles for the root state probabilities
  quantile(x, probs = c(0.025, 0.975))
})
# Print the credible intervals for each root state
print(credible_intervals)
#.       (1,R1)    (2,R1)
#2.5%  0.3923244 0.5000000
#97.5% 0.5000000 0.6076756



# PLOTTING
# Load necessary library
library(ggplot2)

# Prepare data for plotting with modified state labels
root_states <- c("Low Dimorphism", "High Dimorphism")  
estimates <- overall_root_probs
lower_ci <- c(0.3923244, 0.5000000)
upper_ci <- c(0.5000000, 0.6076756)

# Ensure root_states are treated as a factor with specific order
root_states <- factor(root_states, levels = c("Low Dimorphism", "High Dimorphism"))

# Create a data frame for plotting
plot_data <- data.frame(
  RootState = root_states,
  Estimate = estimates,
  LowerCI = lower_ci,
  UpperCI = upper_ci
)

# Plot using ggplot2
ggplot(plot_data, aes(x = RootState, y = Estimate)) +
  geom_point(size = 4, color = "mediumpurple") +  # Plot points for the estimates
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2) +  # Add error bars for CI
  labs(
    x = "Root State",
    y = "Estimate"
  ) +
  theme_minimal() +  # Apply a minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for clarity
    panel.grid = element_blank(),  # Remove grid lines
    axis.title = element_text(size = 14)  # Increase axis title font size
  )






####################################################################################################










###############
###For each tree, estimate the most likely number of x for each node of the phylogeny
###############

#Convert ASR estimates to list of data frames
ASR_states_list <- lapply(ASR_models, function(df) as.data.frame(df$states))

# Create an empty list to store the results
TransDat_neg_rep_list <- vector(mode = "list", length = length(ASR_states_list))

##############
for (i in seq_along(ASR_states_list)) {
  ASR_states <- ASR_states_list[[i]]
  
  #Combine the probability scores for the two different rate classes for each caste number
  #need to check this works 08/01/25
  ASR_states_summed <- data.frame(low = ASR_states$`(1,R1)` , high = ASR_states$`(2,R1)`)
  
  #Identify which WRP is the highest probability for each node
  ASR_states_summed$max_col <- colnames(ASR_states_summed)[my_max_col(ASR_states_summed)] #originally was using the max.col() function however it was making rounding errors
  
  #Assign node number to each row of the data frame, given that node number corresponds to row number in the ASR_states_summed dataframe
  ASR_states_summed$node_number <- paste0("animal.Node", seq_len(nrow(ASR_states_summed)))
  
  #Create new list of data frames with just the node_state and node_number columns
  ASR_states_summed <- ASR_states_summed %>% dplyr::select(max_col, node_number)
  
  #######
  
  #Get data for each species tip and its WRP data
  ##Create 'nodecode' column" which pastes animal and then also the species tip
  sxData$nodecode <- paste("animal", sxData$animal, sep=".")
  # Create new data frame with just the WRP and nodecode columns
  sxData_minimal <- sxData %>% dplyr::select(SD, nodecode)
  #Make the row names equal to nodecode column
  rownames(sxData_minimal) <- sxData_minimal[, "nodecode"]
  
  #Join 'ASR_states_summed' on top of 'sxData_minimal' using rbind
  ##Use lapply to bind sxData_minimal to each element in ASR_states_summed
  ASR_states_summed <- ASR_states_summed %>% dplyr::rename(poly2 = max_col, animal = node_number)
  sxData_minimal <- sxData_minimal %>% dplyr::rename(poly2 = SD, animal = nodecode)
  
  sxData_minimal$poly2 <- ifelse(sxData_minimal$poly2 == 0, "low",
                                 ifelse(sxData_minimal$poly2 == 1, "high", NA))
  
  Poly <- rbind(ASR_states_summed, sxData_minimal)
  
  
  #Turn tree edges into a data frame - displays the parent node (V1) and the two offspring nodes (V2). 653 is the root node number as there are 652 tips in the tree
  TransDat <- as.data.frame(ASR_models[[i]]$phy$edge)
  
  #Name each of the nodes in the phylogeny with correct name
  TransDat$V1name <- paste("animal.Node", TransDat$V1-length(ASR_models[[i]]$phy$tip.label), sep = "")
  
  #Do the same but for the descendant nodes
  TransDat$V2name <- paste("animal.Node", TransDat$V2-length(ASR_models[[i]]$phy$tip.label), sep = "")
  
  #Number each of the species in the tree by their tip number
  treesp <- data.frame(tip.label=paste("animal.", ASR_models[[i]]$phy$tip.label,sep=""), no=1:length(ASR_models[[i]]$phy$tip.label))
  
  #Match species names with the descendant nodes
  TransDat_neg_rep <- data.frame(TransDat, ancestors=NA, descendents=treesp$tip.label[match(TransDat$V2, treesp$no)])
  
  # make sure these are of class 'character'
  TransDat_neg_rep$ancestors <- as.character(TransDat_neg_rep$ancestors)
  TransDat_neg_rep$descendents <- as.character(TransDat_neg_rep$descendents)
  
  # ancestors currently NAs, so give these V1name 
  TransDat_neg_rep$ancestors <- as.character(ifelse(!is.na(TransDat_neg_rep$ancestors),TransDat_neg_rep$ancestors,TransDat_neg_rep$V1name))
  
  # replace any descendents (i.e. those that aren't tips) with V2 name
  TransDat_neg_rep$descendents <- as.character(ifelse(!is.na(TransDat_neg_rep$descendents),TransDat_neg_rep$descendents,TransDat_neg_rep$V2name))
  
  
  ################
  #Step that assigns WRP to ancestor and descendant
  TransDat_neg_rep <- data.frame(TransDat_neg_rep, 
                                 ancPoly=Poly$poly2[match(TransDat_neg_rep$ancestors, Poly$animal)],
                                 desPoly=Poly$poly2[match(TransDat_neg_rep$descendents, Poly$animal)])
  
  #create new column with concatenated values from ancPoly and desPoly
  TransDat_neg_rep$trans <- paste0(TransDat_neg_rep$ancPoly, "_to_", TransDat_neg_rep$desPoly)
  TransDat_neg_rep$trans <- as.factor(TransDat_neg_rep$trans)
  
  #Append TransDat_neg_rep to a list
  TransDat_neg_rep_list[[i]] <- TransDat_neg_rep
}
#######


#############
#Calculate summary statistics for the number of transitions of each type over the 400 trees
#############

#For a single dataframe
table(TransDat_neg_rep_list[[1]]$trans)
# high_to_high  high_to_low  low_to_high   low_to_low 
# 40            6            7           43
#

#For all 400 data frames

# First, apply the table() function to the 'trans' column of each data frame
table_list <- lapply(TransDat_neg_rep_list, function(df) data.frame(table(df$trans)))

##WRP_1_to_WRP_1
#Calculate the mean frequency of occurrences of 'WRP_1_to_WRP_1' across the table_list list of dataframes. 
# If the column 'WRP_1_to_WRP_1' exists in a data frame, extract the frequencies associated with it and compute the mean. 
#If the column doesn't exist in a data frame, assign a frequency of 0 for that data frame.

freq_low_to_high <- extract_frequency(table_list, "low_to_high")
mean(unlist(freq_low_to_high)) # 9.375
median(unlist(freq_low_to_high)) # 7 
range(unlist(freq_low_to_high)) # 3 - 30

##WRP_1_to_WRP_2
freq_high_to_low <- extract_frequency(table_list, "high_to_low")
mean(unlist(freq_high_to_low)) # 7.52
median(unlist(freq_high_to_low)) # 6
range(unlist(freq_high_to_low)) # 3 - 29



##############
## NEW CODE ADDED 18/01/2024
# Combine transition frequencies into a summary dataframe
# Extract unique transition types across all tables
all_transition_types <- unique(unlist(lapply(table_list, function(df) as.character(df$Var1))))

# Create an empty dataframe to store summary statistics
summary_stats <- data.frame(
  Transition = all_transition_types,
  Mean = numeric(length(all_transition_types)),
  Median = numeric(length(all_transition_types)),
  Min = numeric(length(all_transition_types)),
  Max = numeric(length(all_transition_types)),
  stringsAsFactors = FALSE
)

# Calculate statistics for each transition type
for (i in seq_along(all_transition_types)) {
  transition <- all_transition_types[i]
  # Extract frequencies for the current transition type
  freqs <- extract_frequency(table_list, transition)
  freqs_unlist <- unlist(freqs)
  
  # Populate the summary dataframe
  summary_stats[i, "Mean"] <- mean(freqs_unlist)
  summary_stats[i, "Median"] <- median(freqs_unlist)
  summary_stats[i, "Min"] <- min(freqs_unlist)
  summary_stats[i, "Max"] <- max(freqs_unlist)
}

# View the summary dataframe
View(summary_stats)

# Optionally, save the summary dataframe to a CSV file
#write.csv(summary_stats, "transition_summary_statistics.csv", row.names = FALSE)


############# MODES

# Load necessary library
library(dplyr)

# Function to calculate the mode
calculate_mode <- function(x) {
  uniq_vals <- unique(x)
  uniq_vals[which.max(tabulate(match(x, uniq_vals)))]
}

# Extract unique transition types from the first data frame
unique_transitions <- unique(TransDat_neg_rep_list[[i]]$trans) #### modified 17JAN25

# Initialize an empty list to store the mode for each transition type
mode_frequencies <- list()

# Loop through each unique transition type, calculate the mode, and store it
for (transition in unique_transitions) {
  # Extract frequencies for the current transition type across all data frames
  freq_current_transition <- extract_frequency(table_list, transition)
  
  # Flatten the list into a single vector
  freq_vector <- unlist(freq_current_transition)
  
  # Calculate the mode for the current transition type
  mode_freq <- calculate_mode(freq_vector)
  
  # Store the result in the list
  mode_frequencies[[transition]] <- mode_freq
}

# Convert the list of mode frequencies to a dataframe
transition_summary_df <- data.frame(
  Transition = names(mode_frequencies),
  ModeFrequency = unlist(mode_frequencies)
)

# View the resulting dataframe
View(transition_summary_df)


###################################################################################




#############Plotting ASR with transition nodes###########
#Load ggtree
library(ggtree)
library(ggtreeExtra,
        lib.loc = "/drives/4tb/Juliet/miniconda3/envs/my_r/lib/R/library/")

p <- ggtree(ASR_models[[1]]$phy) # + geom_tiplab(size = 0.2)

#Change WRP number column name
sxData_rename <- sxData %>% dplyr::rename(`Queen-Worker Size Dimorphism` = SD)
sxData_rename$`Queen-Worker Size Dimorphism` <- as.factor(sxData_rename$`Queen-Worker Size Dimorphism`)

# OLD
# Add node labels
# 389 tips, 388 nodes
# tips 1 - 389, nodes 390 onwards
# code says start at 390 and create a sequence of numbers for an additional 387 nodes
#ASR_models[[1]]$phy$node.label <- c(seq(390,(390+387))) 

# NEW
# Add node labels
# 49 tips (species), so 48 internal nodes
# tips are numbered 1 - 49, so nodes should be 50 onwards (48)
# code says start at 50 and create a sequence of numbers for an additional 47 nodes
ASR_models[[1]]$phy$node.label <- c(seq(50,(50+47))) ## double check


#################################################### 
# SIMPLE TREE WITH LABELLED NODES FOR ADDING SUBFAMILY LABELS
simple_tree <- ggtree(ASR_models[[1]]$phy, 
                     layout = "fan", 
                     branch.length = "none", 
                     open.angle = 25, 
                     size = 0.2) + 
  geom_rootedge(TRUE, size = 0.1) + 
  geom_tiplab(size = 2, offset = 2) +  # Increase tip label size
  geom_nodelab(size = 2)             # Increase node label size

simple_tree


####################################################   


# FIGURE TREE

sxData_rename$`Queen-Worker Size Dimorphism` <- factor(sxData_rename$`Queen-Worker Size Dimorphism`,
                                                                    levels = c(0,1), 
                                                                    labels = c("low", "high"))

#basic tree
TEST <- ggtree(ASR_models[[1]]$phy, 
                    layout = "fan", 
                    branch.length = "none", 
                    open.angle = 25, 
                    size = 0.2) + 
  geom_rootedge(TRUE, size = 0.1) + 
  #geom_nodelab(size = 2)    
  geom_tiplab(size = 3, offset = 2) 

TEST

nicer_tree_rotated <- rotate_tree(TEST, angle=90)

tree_styled <- nicer_tree_rotated +
  geom_fruit(data=sxData_rename, # Data
             geom=geom_tile, # Plots 'tiles' (squares) onto each tip
             width=1.6, # Alters the length of the tiles
             position=position_identityx(hexpand=16), # Adjusts how far away the tiles are from the tree
             mapping=aes(y=animal, fill=`Queen-Worker Size Dimorphism`), # Analogous to ggplot aes()
             color = "white", # Colour of outline round the tiles (white makes it look like a gap)
             lwd = 0.2, # Width of line between tiles
             linetype = 1, # Default - other numbers make the line dashes, dotted etc.
             axis.params=list( # Add label to the geom_tile - by adding an x-axis
               axis="x",
               text = " ", # Label to plot
               text.size = 3.5, # Size of text
               hjust = 0, # Adjust position of text relative to the geom_tile
               vjust = 0.8,
               text.angle=360,
               line.colour="white"), # Set to white so axis line is not visible
             offset = 4, # Fine-scale adjustment of space between tree & layers
             pwidth=100) + # Width of whole plot
  scale_fill_manual(values=c("mediumpurple1", "darkorchid4")) + # Tile colours
  theme_tree(legend.position = "none") # This line removes the legend

tree_styled

############### ############### ############### ############### ############### ############### ############### 

############### NODE LABELS #####################

# Subset the first dataframe from the TransDat_neg_rep_list list of dataframes for rows with particular transition types

# GAINS OF HIGH DIMORPHISM
low_to_high <- subset(TransDat_neg_rep_list[[1]], grepl("low_to_high", trans))

# LOSSES OF HIGH DIMORPHISM
high_to_low <- subset(TransDat_neg_rep_list[[1]], grepl("high_to_low", trans))

low_to_high_trans <- low_to_high$V2
high_to_low_trans <- high_to_low$V2

library(ggplot2)
library(ggtree)

# Add transitions for losses and gains, and include a legend for nodes
tree_trans <- tree_styled + 
  # (Mapping color and shape for legend)
  ggtree::geom_point2(aes(subset = node %in% low_to_high_trans,
                          shape = 'Gains of high dimorphism', colour = 'Gains of high dimorphism'),
                      size = 3, fill = 'red', alpha = 0.65) +
  
  # (Mapping color and shape for legend)
  ggtree::geom_point2(aes(subset = node %in% high_to_low_trans,
                          shape = 'Losses of high dimorphism', colour = 'Losses of high dimorphism'),
                      size = 3, fill = 'black', alpha = 0.65) +
  
  # Fill scale for the tip tiles (reproductive potential)
  scale_fill_manual(values=c("mediumpurple1", "darkorchid4"),
                    name = "Queen-Worker Size Dimorphism") +
  
  # Color and shape scale for the node points (transitions)
  scale_colour_manual(name = "Transitions",
                      values = c('Losses of high dimorphism' = 'black', 'Gains of high dimorphism' = 'red')) +
  scale_shape_manual(name = "Transitions",
                    values = c('Losses of high dimorphism' = 20, 'Gains of high dimorphism' = 24)) +
  
  # Adding the legend for the plot
  theme(legend.position = "right") # Adjust legend position

tree_trans




################################################################################################
#### ADD SUBFAMILY LABELS

#might need to: remotes::install_github('YuLab-SMU/ggtree')


# specify node at base of subfamily clade
tree_clad_lab <- tree_trans + geom_cladelab(node=65, label="Pheidole\n(big-headed ants)", align=TRUE, fontsize = 3, angle="auto",
                                            offset = 2, offset.text = 0.5, textcolor='black', barcolor='black') 
1
# double check Attine ants
tree_clad_lab <- tree_clad_lab + geom_cladelab(node=67, label="Fungus-growing \nants (including \nAtta and\nAcromyrmex)", align=TRUE, fontsize = 3, angle="auto",
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=300.  , label="Carebara", align=TRUE, fontsize = 3, angle='auto',
 #                                              offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=84, label="Camponotus &\nCalomyrmex \n(including \ncarpenter ants)", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=87, label="Cataglyphis", align=TRUE, fontsize = 2, angle='auto',
 #                                              offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=88, label="Formica", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=90, label="Lasius", align=TRUE, fontsize = 2, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=70, label="Pogonomyrmex", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=61, label="Temnothorax", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=403, label="Crematogaster", align=TRUE, fontsize = 3, angle='auto',
 #                                              offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1, label="Ponerinae", align=TRUE, fontsize = 3, angle='auto',
 #                                              offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')
## need to check Ponerines

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=707, label="Dorylinae (army ants)", align=TRUE, fontsize = 3, angle='auto',
 #                                              offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')


#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=687, label="Dolichoderinae", align=TRUE, fontsize = 3, angle='auto',
 #                                              offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=582, label="Solenopsidini \n(including fire ants)", align=TRUE, fontsize = 3, angle='auto',
 #                                              offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=600, label="Stenammini", align=TRUE, fontsize = 3, angle='auto',
#                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=76, label="Myrmicini", align=TRUE, fontsize = 3, angle='auto',
                                              offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=622, label="Ectatommini", align=TRUE, fontsize = 3, angle='auto',
 #                                              offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1011, label="Polyrhachis", align=TRUE, fontsize = 3, angle='auto',
#                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1101, label="Brachymyrmex", align=TRUE, fontsize = 3, angle='auto',
 #                                              offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=697, label="Pseudomyrmecini", align=TRUE, fontsize = 3, angle='auto',
#                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=460, label="Metapone", align=TRUE, fontsize = 2, angle='auto',
 #                                              offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=681, label="Oecophylla", align=TRUE, fontsize = 2, angle='auto',
 #                                              offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=461, label="Tetramorium", align=TRUE, fontsize = 3, angle='auto',
 #                                              offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=89, label="Plagiolepis", align=TRUE, fontsize = 2, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=819, label="Strumigenys", align=TRUE, fontsize = 3, angle='auto',
#                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=698, label="Myrmeciini", align=TRUE, fontsize = 3, angle='auto',
#                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

#tree_clad_lab <- tree_clad_lab + geom_cladelab(node=879, label="Daceton & \nOrectognathus", align=TRUE, fontsize = 3, angle='auto',
 #                                              offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab

with_species <- tree_clad_lab + geom_tiplab(size = 1, offset = 2) 
with_species
