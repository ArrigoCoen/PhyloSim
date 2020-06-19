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
## Extras packa
using DelimitedFiles
using PhyloNetworks
using Printf
using Plots


using Dates
# import Dates
using Statistics
using RCall

## Extra libraries
using Pkg
# Pkg.add("RCall")


## FUNCTIONS

using Random

"""
    testing if a one_taxa string is valid
        This is a function in progress
"""
function is_valid_1taxa(one_taxa::String)
    test = true
    if(tryparse(Float64, SubString(one_taxa, 1,1)) !== nothing)
        error("taxa should not start with a number")
    end
    return test
end

"""
    Simulation of a phylogenetic tree topology
        We start with a two taxa tree and start to paste subtrees by choosing
        randomly taxa to paste.
        Currently the branch lengths are uniformly distributed (0,1)
        Returs a the topology in Newick format
"""
function random_tree_topology(taxa)
    shuffle!(taxa)
    # Starting the tree
    x = Any["(",taxa[1],",",taxa[2],")",";"]

    for i in 3:length(taxa)
        # taxa already in the tree and taxa not in the tree
        used_taxa = taxa[1:(i-1)]
        left_taxa = taxa[i:end]
        # choosing a taxa in the tree to add new branches
        taxa_to_change = rand(used_taxa)
        # possition to add branches
        position_of_change = findfirst(x.==taxa_to_change)
        # updating tree
        insert!(x,position_of_change+1,")")
        insert!(x,position_of_change,",")
        insert!(x,position_of_change,taxa[i])
        insert!(x,position_of_change,"(")
    end
    return x
end

"""
    Adds random branch lengths to a phylogenetic tree topology
        Currently the branch lengths are uniformly distributed (0,1)
        Returs a the topology in Newick format
"""
function add_random_bl(x,taxa)
    # We get the index of the taxa and right parenthesis (to add the branch lengths)
    idx_right_parenthesis = findall(x.==")")
    idx_taxa = zeros(Int64,length(taxa))
    for i in 1:length(taxa)
        idx_taxa[i] = findfirst(x.==taxa[i])
    end
    idx_taxa_and_right_p = vcat(idx_taxa, idx_right_parenthesis)
    # To paste the new branch length we need to sort in decreasing order
    sort!(idx_taxa_and_right_p, rev=true)
    # Random uniform branch length
    bl = rand(length(idx_taxa_and_right_p))
    # Adding the branch length
    for i in 1:length(idx_taxa_and_right_p)
        insert!(x,idx_taxa_and_right_p[i]+1,bl[i])
        insert!(x,idx_taxa_and_right_p[i]+1,":")
    end
    return x
end


"""
    Given a topology and branch lengths, this function change the branch lengths
    in a way that the tree is now (insert word to say that all the paths from a
    leaf to the root are of the same length)
        Only keeps the topology and with a stick break methodology gives the
        new branch lengths.
        Returs a the topology in Newick format
"""
function balance_a_tree(x,total_bl)
    # Case of taxa-colon-branch length-semicolon; eg A:0.32;
    if length(x)==4
        x[3] = total_bl[1]
        return x
    end

    # Case of two taxa
    if length(x)==12
      partition = rand(1)
      # bl = total_bl.*partition
      x[[4,8]] .= total_bl.*partition
      aux = total_bl.*(Float64[1]-partition)
      x[11] = aux[1]
      return x
    end

    # We need to find the place of the comma that divides the tree
    aux = (x.==")") - (x.=="(")
    pop!(aux)
    popfirst!(aux)
    central_comma_place = findfirst(cumsum(aux).==0) +4

    # Stuff of the left tree
    left_tree = vcat(x[2:(central_comma_place-1)],";")
    partition = rand(1)
    left_tree = balance_a_tree(left_tree,total_bl.*partition)

    # Stuff of the right tree
    right_tree = vcat(x[(central_comma_place+1):(length(x)-4)],";")
    right_tree[length(right_tree)-1] = left_tree[length(left_tree)-1]
    right_tree = balance_a_tree(right_tree,total_bl.*partition)

    # Pasting the left and right trees
    pop!(left_tree)
    pop!(right_tree)
    x = vcat("(",left_tree,",",right_tree,")",":",total_bl.*(Float64[1]-partition),";")

    return x
end
"""
    Simulating a phylogenetic tree
        Returs a the topology in Newick format
"""
function Sim_tree(taxa,total_bl)
    # Getting the topology
    x = random_tree_topology(taxa)
    # Adding some branch lengths
    x = add_random_bl(x,taxa)
    # Fixing the branch lengths to be "equal"
    x = balance_a_tree(x,total_bl)
    return x
end

## aoestuh
