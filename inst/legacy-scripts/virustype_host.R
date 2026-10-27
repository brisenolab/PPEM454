# Types of viruses


nucleicacid <- c("ssRNA+", "ssRNA-", "dsRNA", "cRNA", 'RNA-RT', "ssDNA", "dsDNA", "DNA-RT") 
host <- c("archeobacteria", "bacteria", "algae", "protozoa", "oomycete", "fungi", "plant", "animal")
# matrix(nucleicacid, host)

M <- matrix(c(
  # ssRNA+
  0, 1, 1, 0, 0, 1, 1, 1,
  # ssRNA-
  0, 0, 0, 0, 1, 1, 1, 1,
  # dsRNA
  0, 1, 0, 1, 1, 1, 1, 1,
  # cRNA  (circular RNA viruses; e.g., HDV in animals; circular RNA viruses reported in fungi)
  0, 0, 0, 0, 0, 1, 1, 1,
  # retroRNA (RNA-RT; Baltimore VI: Retroviridae)
  0, 0, 0, 0, 0, 0, 0, 1,
  # ssDNA
  1, 1, 1, 0, 0, 1, 1, 1,
  # dsDNA
  1, 1, 1, 1, 0, 0, 1, 1,
  # DNA-RT
  0, 0, 0, 0, 0, 0, 1, 0  
), nrow = length(nucleicacid), byrow = TRUE,
dimnames = list(nucleicacid, host))

library(bipartite)

virusinfection <- sortweb(M, sort.order = "dec")

# Basic bipartite plot
plotweb(virusinfection, 
        method = "normal",   # layout
        col.interaction = "grey70",
        text.rot = 45,       # rotate host labels
        labsize = 2,
        # low.spacing = 1, 
        # high.spacing = 1,
        adj.high=.1, 
        adj.low=1,
        ybig = 1.5,          # spacing between the two levels
        col.high = "tomato", # virus nodes
        col.low = "steelblue") # host nodes

###--- Virus genome and virion

nucleicacid <- c("ssRNA+", "ssRNA-", "dsRNA", "cRNA", "RNA-RT", "ssDNA", "dsDNA", "DNA-RT")
virion_structure <- c("capsid", "pleomorphic", "naked")   # (typo fixed: "pleomorphic")

M <- matrix(c(
  # ssRNA+
  1, 1, 1,
  # ssRNA-
  0, 1, 0,
  # dsRNA
  1, 1, 1,
  # cRNA (circular RNA viruses)
  0, 1, 1,
  # retroRNA (Retroviridae)
  0, 1, 0,
  # ssDNA
  1, 1, 0,
  # dsDNA
  1, 1, 0,
  # DNA-RT
  1, 0, 0
), nrow = length(nucleicacid), byrow = TRUE,
dimnames = list(nucleicacid, virion_structure))


M
# Basic bipartite plot
plotweb(M, 
        method = "normal",   # layout
        col.interaction = "grey70",
        text.rot = 45,       # rotate host labels
        labsize = 2,
        # low.spacing = 1, 
        # high.spacing = 1,
        ybig = 1.5,          # spacing between the two levels
        col.high = "tomato", # virus nodes
        col.low = "steelblue") # host nodes


# install.packages(c("igraph"))
library(igraph)

# M is your incidence matrix (rows = genome types, cols = virion structures)
g <- graph_from_biadjacency_matrix(M)  # bipartite graph

# Colors/shapes by mode
V(g)$type <- bipartite_mapping(g)$type  # FALSE = rows, TRUE = cols
V(g)$color <- ifelse(V(g)$type, "#2b8cbe", "#e34a33")  # cols vs rows
V(g)$shape <- ifelse(V(g)$type, "square", "circle")

# Circle layout
lo <- layout_in_circle(g, order = V(g)$name)

plot(g, layout = lo,
     vertex.size = 18, vertex.label.cex = 0.9,
     edge.color = adjustcolor("grey40", 0.6), edge.width = 2,
     vertex.label.family = "sans")
