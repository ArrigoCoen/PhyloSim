#-------------------------------------------------------------------------#
#                                                                         #
#                               Simulating_nets                           #
#                                                                         #
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
# OBSERVATIONS:

# 1. This code is for simulating phylogenetic networks in a Newick format
# 2. This code is a work in progress
# 3. The functions are in the file Fn_PhyloSim
# 4. The principal example is on "Example of a random tree generation" section
#
#-------------------------------------------------------------------------#

## Packages
# using Plots
using PhyloNetworks
using PhyloPlots

## Example of a random tree generation


# Including functions
include("Fn_PhyloSim.jl")
# Fixing seed
# Random.seed!(42)

# Different taxa examples
taxa = ["1ba","b","c","d"]
taxa = ["A","B"]
taxa = collect('a':'d')
taxa = collect('a':'z')

# Total branch lengths
total_bl = 10

# Return a network of a tree
is_network = true

# Generating a random tree topology and branch length in Newikc format
x = Sim_tree(taxa,total_bl,is_network)

print(join(x))

net = readTopology(join(x))
tipLabels(net)

# plot(net,:R)

plot(net,:R,showEdgeLength=true,showGamma=true)
# plot(net,:R,showEdgeLength=true)
net
## Some info about the network
printEdges(net)
