library(ape)
library(treeio)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read the NEXUS file with annotations
trees <- read.beast("mapper.trees")


fixed_tree <- ape::read.tree(text = "(t5:5.7,((t1:1,t2:2):1,(t3:3,t4:4):0.5):1.3):0.0;")
n_tips <- length(fixed_tree$tip.label)
n_nodes <- fixed_tree$Nnode
internal_ids <- (n_tips + 1):(n_tips + n_nodes)
internal_nodes <- tree@data$node[!tree@data$node %in% c(1,2,3,4,5)]
n_types <- 2

node_heights <- node.depth.edgelength(fixed_tree)[internal_ids]
pi_empirical <- matrix(0, nrow = length(internal_nodes), ncol = n_types)
rownames(pi_empirical) <- internal_nodes



type_for_node <- function(node_data, node_id) {
  match_row <- which(node_data$node == as.character(node_id))
  node_data$type[match_row]
}


# Loop over trees
for (tree in trees) {
  node_data <- tree@data
  
  for (node_id in internal_nodes) {
    type <- type_for_node(node_data, node_id)
    if (!is.na(type)) {
      pi_empirical[node_id, type + 1] <- pi_empirical[node_id, type + 1] + 1
    }
  }
}


# Normalize
pi_empirical <- pi_empirical / length(trees)
colnames(pi_empirical) <- c("π0", "π1")
print(pi_empirical)

# Last row = root probabilities



