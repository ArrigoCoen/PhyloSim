#-------------------------------------------------------------------------#
#                                                                         #
#                               Simulating_nets                           #
#                                                                         #
#-------------------------------------------------------------------------#

# V1 18 may 2020
#-------------------------------------------------------------------------#
# OBSERVATIONS:

# 1. Bla bla bla
# 2. Bla bla bla
# 3. Bla bla bla
#
#-------------------------------------------------------------------------#

## Pacakges
using LinearAlgebra
using PhyloPlots
using DelimitedFiles
using PhyloNetworks
using Printf
using Plots
using Dates
using Statistics
using RCall
import Plots
## Extra libraries
# using Pkg
# Pkg.add("RCall")

## Example of a random tree generation

# Including functions
include("Fn_PhyloSim.jl")
# Fixing seed
Random.seed!(42)

# Different taxa for testing
taxa = ["1ba","b","c","d"]
taxa = ["A","B"]
taxa = collect('a':'d')
taxa = collect('a':'z')

# Total branch lengths
total_bl = 10

# Generating a random tree topology and branch length in Newikc format
x = Sim_tree(taxa,total_bl)
print.(x); println("\n\n\n")

## To plot this tree

nettext = join(x)
net = readTopology(nettext)
p = plot(net, showEdgeLength=true)
draw(SVG("Plot_random_tree.svg"), p)
