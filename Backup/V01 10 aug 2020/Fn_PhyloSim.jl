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
# using PhyloNetworks

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
    Addition of a random hybrid node to a phylogenetic tree
        We start a tree "x" with taxa "taxa" and select two taxa to add the hybrid
        node.
        Currently the branch length of the new hybrid and its gamma parameter are
        uniformly distributed (0,1). The first node that is selected we currently
        TRANSFORM it to a hybrid node. For the second node we only paste the
        other hybrid node.
        For instance, if we have the tree "((a,b),(c,d));" and select "b" and "c",
        we get "((a,#H1),((c)#H1,d));", without considering the branch lengths and
        the gamma parameter.
        Currently, sometimes the hybrid is in the same branch.
        Returs the tree with the added hybrid, in Newick format.

        Future improvements:
            - Extra variable with the name of the hybrid, now is H1
            - Correcting that in some simulations the hybrid nodes are in the same
            path an this probuces a weird hybrid that follows the same path as
            the original tree
"""
function add_random_hybrid(x,taxa)
   # sampling two taxa to insert the hybrid
   taxa_to_change1,taxa_to_change2 = taxa[randperm(length(taxa))[1:2]]
   # We get the branch length and the gamma parameter for the new hybrid
   bl_hybrid = rand(1)
   gamma = rand(1)
   # Finding the index of the taxa to change
   idx_node1 = findfirst(x.==taxa_to_change1)
   idx_node2 = findfirst(x.==taxa_to_change2)
   # We allways prefere the index of node 1 to be smaller
   if idx_node1 > idx_node2
      idx_node1,idx_node2 = idx_node2,idx_node1
   end
   # Now we can make the changes on the array x
   x[idx_node2] = "#H1"

   insert!(x,idx_node1+3,gamma[1])
   insert!(x,idx_node1+3,":")
   insert!(x,idx_node1+3,":")
   insert!(x,idx_node1+3,bl_hybrid[1])
   insert!(x,idx_node1+3,":")
   insert!(x,idx_node1+3,"#H1")
   insert!(x,idx_node1+3,")")
   insert!(x,idx_node1,"(")
   return x
end

"""
    Simulating a phylogenetic tree/network
        Returs a the topology in Newick format

        Future improvements:
            - To add the hybrid a node is eliminated, this could be solved by
            since the begining ask the name of the hybrid and past it as part of
            the tree, but remember that this change also affects the function
            add_random_hybrid; in particular, since this function chose randomly
            which node of the tree is transformed to hybrid.
"""
function Sim_tree(taxa,total_bl,is_network)
    # Getting the topology
    x = random_tree_topology(taxa)
    # Adding some branch lengths
    x = add_random_bl(x,taxa)
    # Fixing the branch lengths to be "equal"
    x = balance_a_tree(x,total_bl)
    # Adding a random hybrid
    if is_network
        x = add_random_hybrid(x,taxa)
    end
    return x
end

## Funcitons of njWillems2014

function Julia_willems2014(x,name_file_dist_txt,
                           name_file_Input_Matrix_Distances,
                           name_file_Parameters_Distances,
                           name_file_new_Network_Distances,
                           folder_output_distances_willems)
   # Folder to save the output of willems
   output_willems = folder_output_distances_willems*"/"*name_file_new_Network_Distances
   # Getting the netwokr from the Newick format
   net = readTopology(join(x))
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

   CSV.write(name_file_dist_txt, df, writeheader=false, delim=' ')


   # the first row with the number of tips
   f=readlines(name_file_dist_txt)
   g = open(name_file_Input_Matrix_Distances,"w")
   write(g,"$n \n")

   for l in f
      write(g,l)
      write(g,"\n")
   end
   close(g)



   # This "run" generates the file "Network_Distances.txt"
   text_to_run = `./njwillems2014 $name_file_Input_Matrix_Distances $name_file_Parameters_Distances`
   run(text_to_run)

   # We copy the file "Network_Distances.txt", that generates njwillems2014
   cp("Network_Distances.txt", output_willems, force=true)

   net_dist = readdlm("Network_Distances.txt")
   # rm(name_file_new_Network_Distances)

   # txt_file = open("Network_Distances.txt","r")
   # readline(txt_file)
   # close(txt_file)

   # print(join(net_dist[1:10]))
   # println(net_dist)

   idx_Main_text = findfirst(net_dist .== "Main")
   println(net_dist[CartesianIndex(1, 1):(idx_Main_text)])

   # println((net_dist[1:20]))
   return net_dist
end
## End of functions
