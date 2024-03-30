library(readr)
library(readxl)
library(eulerr)
library(igraph)
Bocklandt_CpGs <- read_excel("Bocklandt CpGs.xls")
Bocklandt_CpGs <- Bocklandt_CpGs[, c(1)]
colnames(Bocklandt_CpGs)[1] <- "CPG"
Voisin_CpGs <- read_excel("Voisin CpGs.xlsx")
Voisin_CpGs <- Voisin_CpGs[-1, c(1)]
colnames(Voisin_CpGs)[1] <- "CPG"
Hannum_CpGs <- read_excel("Hannum CpGs.xlsx")
Hannum_CpGs <- Hannum_CpGs[, c(1)]
colnames(Hannum_CpGs) = "CPG"
Levine_CpGs <- read.csv("Levine CpGs.csv", stringsAsFactors=FALSE)
Levine_CpGs <- Levine_CpGs[-1, 1, drop=FALSE]
names(Levine_CpGs)[1] = "CPG"
Horvath1_CpGs = read.csv("Horvath 1 CpGs.csv", stringsAsFactors=FALSE)
Horvath1_CpGs = Horvath1_CpGs[-(1:3), 1, drop=FALSE]
colnames(Horvath1_CpGs)[1] = "CPG"
Weidner_CpGs <- read_excel("Weidner CpGs.xlsx")
Vidal_Bralo_CpGs <- read_excel("Vidal-Bralo CpGs.xlsx")
Horvath3_CpGs <- read.csv("Horvath3 CpGs.csv", stringsAsFactors=FALSE)
Horvath3_CpGs = Horvath3_CpGs[-1, 1, drop=FALSE]
colnames(Horvath3_CpGs)[1] = "CPG"
McEwen_CpGs <- read_excel("McEwen CpGs.xlsx")



list_of_dfs <- list(Bocklandt_CpGs, Hannum_CpGs, Horvath1_CpGs, Horvath3_CpGs,
                    Levine_CpGs, Vidal_Bralo_CpGs, Voisin_CpGs, McEwen_CpGs, Weidner_CpGs)

list_of_cpgs <- lapply(list_of_dfs, function(df) {
  unique(df$CPG)
})
edges <- data.frame(from = character(), to = character(), weight = numeric())
for (i in 1:(length(list_of_dfs) - 1)) {
  for (j in (i + 1):length(list_of_dfs)) {
    # Find the intersection of CpG sites between datasets i and j
    intersection <- intersect(list_of_cpgs[[i]], list_of_cpgs[[j]])
    
    # Add the edge if there is an intersection
    if (length(intersection) > 0) {
      edges <- rbind(edges, data.frame(from = paste0("Dataset", i), 
                                       to = paste0("Dataset", j), 
                                       weight = length(intersection)))
    }
  }
}

graph <- graph_from_data_frame(edges, directed = FALSE)

node_sizes <- sapply(list_of_cpgs, length)
max_size <- 100
min_size <- 10
normalized_node_sizes <- (node_sizes - min(node_sizes)) / (max(node_sizes) - min(node_sizes)) * (max_size - min_size) + min_size
V(graph)$name <- c("Bocklandt", "Hannum", "Horvath1", "Horvath3", "Levine", "Vidal-Bralo", "Voisin", "McEwen", "Weidner")

graph_layout <- layout_in_circle(graph)
scaling_factor <- 4  # Adjust this factor as needed
graph_layout <- graph_layout * scaling_factor
plot(graph, 
     layout = graph_layout,
     edge.width = E(graph)$weight / 5,  # Use the weight for edge width
     vertex.size = normalized_node_sizes / 1.6,
     vertex.label = V(graph)$name,
     vertex.label.cex = 0.7,
     vertex.color = "lightblue",
     edge.color = "grey",
     vertex.frame.color = NA)

png(filename = "connections.png", width = 4, height = 4, units = 'in', res = 300)
plot(graph, 
     layout = graph_layout,
     edge.width = E(graph)$weight / 5,  # Use the weight for edge width
     vertex.size = normalized_node_sizes / 1.6,
     vertex.label = V(graph)$name,
     vertex.label.cex = 0.7,
     vertex.color = "lightblue",
     edge.color = "grey",
     vertex.frame.color = NA)
dev.off()
