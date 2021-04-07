#ok wow the installation for ggtree is weird and intense
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
browseVignettes("ggtree")

library(ggplot2)
library(ggtree)
library(ape)

#read in data and phylogeny

data1 = read.csv("cumulative_data_final.csv", header=T)
str(data1)
colnames(data1)[3]="id"
data = data1[,-c(1:2)]

tree = read.nexus("consensusTree_10kTrees_Primates_allstreps.nex")
plot(tree)
str(tree$tip.label)

tree_plot = ggtree(tree) + geom_tiplab() 

tree_plot + coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 120, 6, 6))


# Make the original tree plot
p = ggtree(tree) + geom_tiplab(size=2)+ coord_cartesian(clip = 'off')

# Make a second plot with the original, naming the new plot "dot", 
# using the data you just created, with a point geom.
p2 = facet_plot(p, panel="Body mass (kg)", data=data, geom=geom_point, 
                 aes(x=age_sex_meanbodyweight, color = Mating_System)) 

p2 + theme_tree2(panel.spacing = unit(5, "lines"))





