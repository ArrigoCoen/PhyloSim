#-------------------------------------------------------------------------#
#                                                                         #
#                               Simulating_nets                           #
#                                                                         #
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
# OBSERVATIONS:

# 1. This code is for simulating phylogenetic networks
# 2. This code is a work in progress
# 3. The functions are in the file Fn_PhyloSim
# 4. The principal example is on "Example of a random tree generation" section
#
#-------------------------------------------------------------------------#

using Plots
using PhyloNetworks

using PhyloPlots
## Needed

# using Gadfly, Cairo, Fontconfig
# using PhyloPlots
# using PhyloNetworks
#
# ## Extras packa
# using LinearAlgebra
# using DelimitedFiles
# using PhyloNetworks
# using Printf
# using Dates
# using Statistics
# using RCall

## Changing packages
using Pkg
# Pkg.add("PhyloPlots")

# Pkg.update() # get all latest versions
# Pkg.update("PhyloPlots")
# Pkg.update("PhyloNetworks")
# Pkg.rm("Plots")
# Pkg.rm("PhyloPlots")

## Example of a random tree generation

# Including functions
include("Fn_PhyloSim.jl")
# Fixing seed
Random.seed!(42)

# Different taxa examples
taxa = ["1ba","b","c","d"]
taxa = ["A","B"]
taxa = collect('a':'d')
taxa = collect('a':'z')

# Total branch lengths
total_bl = 10

# Generating a random tree topology and branch length in Newikc format
x = Sim_tree(taxa,total_bl)
print(join(x))

net = readTopology(join(x))
tipLabels(net)

plot(net,:R)
net
## Some info about the network
printEdges(net)
