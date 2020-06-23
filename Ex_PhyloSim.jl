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

## Example from
# http://crsl4.github.io/PhyloNetworks.jl/latest/man/inputdata/#Input-for-SNaQ-1

# To show the info on the file raxmltrees.tre, run the next
# less("raxmltrees.tre") # to exit use "q"

raxmltrees = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","raxmltrees.tre")
genetrees = readMultiTopology(raxmltrees);
genetrees[3]

plot(genetrees[3], :R); # tree for 3rd gene
plot(genetrees[3], showEdgeLength=true); # tree for 3rd gene
plot(genetrees[3]); # tree for 3rd gene



## tanehu
using Plots
x = 1:10; y = rand(10); # These are the plotting data
plot(x, y)

## To plot this tree
x = Sim_tree(taxa,total_bl)
nettext = join(x)
net = readTopology(nettext)

plot(net)

p = plot(net, showEdgeLength=true)
draw(SVG("Plot_random_tree.svg"),   p)


## Runing

x = Sim_tree(taxa,total_bl)

nettext =  "(#H1:0.629867076426745::0.667790557153057,(((2:0.511301759033819,7:0.511301759033819):0.0982665963226665,((4:0.0106596293945877,8:0.0106596293945877)#H1:1.45281970574055::0.533109405718278,6:0.177021141802196):0.432547213554289):0.466181909365447,(3:0.732115608033462,5:0.732115608033462):0.343634656688471):2.60910392607468):6.31514580920339;"
nettext = "(10:9.6,(#H2:2.9::0.3,(1:7.2,(2:6.0,(((9:0.4)#H1:5.0::0.8,(3:4.4,(4:3.5,((5:0.2,6:0.2):2.1,(7:1.4,(8:0.4,#H1:0.0::0.2):1.0):0.9):1.2):0.9):1.0):0.1)#H2:0.5::0.7):1.2):1.2):1.2);"


net = readTopology(nettext)
p = plot(net, showEdgeLength=true)
draw(SVG("n10.svg"), p)


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
