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
using Printf
using Plots
using Dates
using Statistics
using RCall
using Plots

using DataFrames
using CSV


using PhyloPlots
using DelimitedFiles
using PhyloNetworks

## Extra libraries
# using Pkg
# Pkg.add("RecipesBase")


 # import RecipesBase: plot


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


## Runing

x = Sim_tree(taxa,total_bl)

nettext =  "(#H1:0.629867076426745::0.667790557153057,(((2:0.511301759033819,7:0.511301759033819):0.0982665963226665,((4:0.0106596293945877,8:0.0106596293945877)#H1:1.45281970574055::0.533109405718278,6:0.177021141802196):0.432547213554289):0.466181909365447,(3:0.732115608033462,5:0.732115608033462):0.343634656688471):2.60910392607468):6.31514580920339;"
nettext = "(10:9.6,(#H2:2.9::0.3,(1:7.2,(2:6.0,(((9:0.4)#H1:5.0::0.8,(3:4.4,(4:3.5,((5:0.2,6:0.2):2.1,(7:1.4,(8:0.4,#H1:0.0::0.2):1.0):0.9):1.2):0.9):1.0):0.1)#H2:0.5::0.7):1.2):1.2):1.2);"


net = readTopology(nettext)
p = plot(net, showEdgeLength=true)
# draw(SVG("resources/n10.svg"), p)


## aoheu

# Distance matrix
Sigma = sharedPathMatrix(net)
Vy = Sigma[:Tips]

# From Cecil
M = pairwiseTaxonDistanceMatrix(net, keepInternal=false)
tipLabels(net)

# Generating file to run njWillems2014
ind = parse.(Int, tipLabels(net))
# ind =  tipLabels(net)



Mext = hcat(ind,M)
n = length(tipLabels(net))


df = DataFrame(Mext)
df[!,:x1] = convert.(Int,df[!,:x1])

CSV.write("dist-n10.txt", df, writeheader=false, delim=' ')


# the first row with the number of tips
f=readlines("dist-n10.txt")
g = open("Input_Matrix_Distances.txt","w")
write(g,"$n \n")

for l in f
   write(g,l)
   write(g,"\n")
end
close(g)


rm("Network_Distances.txt")
run(`./njwillems2014 Input_Matrix_Distances.txt Parameters_Distances.txt`)




## atosenhua

using PhyloNetworks, PhyloPlots
using Gadfly, Cairo, Fontconfig
nettext = "(10:9.6,(#H2:2.9::0.3,(1:7.2,(2:6.0,(((9:0.4)#H1:5.0::0.8,(3:4.4,(4:3.5,((5:0.2,6:0.2):2.1,(7:1.4,(8:0.4,#H1:0.0::0.2):1.0):0.9):1.2):0.9):1.0):0.1)#H2:0.5::0.7):1.2):1.2):1.2);"
net = readTopology(nettext)
p = plot(net, showEdgeLength=true)
draw(SVG("resources/n10.svg"), p)

## atosenhua
