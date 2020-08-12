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
    if (tryparse(Float64, SubString(one_taxa, 1, 1)) !== nothing)
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
    x = Any["(", taxa[1], ",", taxa[2], ")", ";"]

    for i = 3:length(taxa)
        # taxa already in the tree and taxa not in the tree
        used_taxa = taxa[1:(i-1)]
        left_taxa = taxa[i:end]
        # choosing a taxa in the tree to add new branches
        taxa_to_change = rand(used_taxa)
        # possition to add branches
        position_of_change = findfirst(x .== taxa_to_change)
        # updating tree
        insert!(x, position_of_change + 1, ")")
        insert!(x, position_of_change, ",")
        insert!(x, position_of_change, taxa[i])
        insert!(x, position_of_change, "(")
    end
    return x
end

"""
    Adds random branch lengths to a phylogenetic tree topology

        x = a topology in Newick format; eg. Any["(", 'b', ",", "(", 'c', ",", 'a', ")", ")", ";"]
        taxa = names of taxa; eg. Any['b', 'c', 'a']
        digits_bl = an integer that restricts the number of digits

        Currently the branch lengths are uniformly distributed (0,1).
        Returs a the topology in Newick format
"""

function add_random_bl(x, taxa, digits_bl)
    # We get the index of the taxa and right parenthesis (to add the branch lengths)
    idx_right_parenthesis = findall(x .== ")")
    idx_taxa = zeros(Int64, length(taxa))
    for i = 1:length(taxa)
        idx_taxa[i] = findfirst(x .== taxa[i])
    end
    idx_taxa_and_right_p = vcat(idx_taxa, idx_right_parenthesis)
    # To paste the new branch length we need to sort in decreasing order
    sort!(idx_taxa_and_right_p, rev = true)
    # Random uniform branch length
    # bl = rand(length(idx_taxa_and_right_p)) # before
    bl = round.(rand(length(idx_taxa_and_right_p)), digits=digits_bl)
    # Adding the branch length
    for i = 1:length(idx_taxa_and_right_p)
        insert!(
            x,
            idx_taxa_and_right_p[i] + 1,
            round(bl[i], digits = digits_bl),
        )
        insert!(x, idx_taxa_and_right_p[i] + 1, ":")
    end
    return x
end

# OLD VERSION WITHOUT DIGITS RESTRICTION
# function add_random_bl(x,taxa)
#     # We get the index of the taxa and right parenthesis (to add the branch lengths)
#     idx_right_parenthesis = findall(x.==")")
#     idx_taxa = zeros(Int64,length(taxa))
#     for i in 1:length(taxa)
#         idx_taxa[i] = findfirst(x.==taxa[i])
#     end
#     idx_taxa_and_right_p = vcat(idx_taxa, idx_right_parenthesis)
#     # To paste the new branch length we need to sort in decreasing order
#     sort!(idx_taxa_and_right_p, rev=true)
#     # Random uniform branch length
#     bl = rand(length(idx_taxa_and_right_p))
#     # Adding the branch length
#     for i in 1:length(idx_taxa_and_right_p)
#         insert!(x,idx_taxa_and_right_p[i]+1,bl[i])
#         insert!(x,idx_taxa_and_right_p[i]+1,":")
#     end
#     return x
# end




"""
    Given a topology and branch lengths, this function change the branch lengths
    in a way that the tree is now ultrametric meaning that all the paths from a
    leaf to the root are of the same length)

        x = a topology in Newick format; eg. Any["(", 'b', ",", "(", 'c', ",", 'a', ")", ")", ";"]
        total_bl = total length of the tree from root to any leaf
        digits_bl = an integer that restricts the number of digits
        is_unrooted = true or false

        Only keeps the topology and with a stick break methodology gives the
        new branch lengths. The restriction on the digits is useful only when
        simulating small trees (eg. number of taxa<5), since for balancing it
        starts to need more and more digits.
        Returs a the a ultrametric rooted tree in Newick format
"""
function balance_a_tree(x, total_bl, digits_bl,is_unrooted)
    # Case of taxa-colon-branch length-semicolon; eg A:0.32;
    if length(x) == 4
        x[3] = total_bl[1]
        return x
    end

    # Case of two taxa
    if length(x) == 12
        # partition = rand(1) # Before
        partition = round(rand(1)[1], digits = digits_bl)
        # bl = total_bl.*partition
        x[[4, 8]] .= total_bl .* partition
        # aux = total_bl.*(Float64[1]-partition) # Before
        aux = total_bl .* (1 - partition)
        x[11] = aux[1]
        return x
    end

    # We need to find the place of the comma that divides the tree
    aux = (x .== ")") - (x .== "(")
    pop!(aux)
    popfirst!(aux)
    central_comma_place = findfirst(cumsum(aux) .== 0) + 4

    # partition = rand(1) # before)
    partition = round(rand(1)[1], digits = digits_bl)
    if is_unrooted
        partition = 1
        is_unrooted = false
    end

    # Stuff of the left tree
    left_tree = vcat(x[2:(central_comma_place-1)], ";")
    left_tree = balance_a_tree(left_tree, total_bl .* partition, digits_bl,is_unrooted)

    # Stuff of the right tree
    right_tree = vcat(x[(central_comma_place+1):(length(x)-4)], ";")
    right_tree[length(right_tree)-1] = left_tree[length(left_tree)-1]
    right_tree = balance_a_tree(right_tree, total_bl .* partition, digits_bl,is_unrooted)

    # Pasting the left and right trees
    pop!(left_tree)
    pop!(right_tree)
    # x = vcat("(",left_tree,",",right_tree,")",":",total_bl.*(Float64[1]-partition),";") # Before
    x = vcat(
        "(",
        left_tree,
        ",",
        right_tree,
        ")",
        ":",
        total_bl .* (1 - partition),
        ";",
    )

    return x
end


# OLD VERSION WITHOUT DIGITS RESTRICTION
# function balance_a_tree(x,total_bl)
#     # Case of taxa-colon-branch length-semicolon; eg A:0.32;
#     if length(x)==4
#         x[3] = total_bl[1]
#         return x
#     end
#
#     # Case of two taxa
#     if length(x)==12
#       partition = rand(1)
#       # bl = total_bl.*partition
#       x[[4,8]] .= total_bl.*partition
#       aux = total_bl.*(Float64[1]-partition)
#       x[11] = aux[1]
#       return x
#     end
#
#     # We need to find the place of the comma that divides the tree
#     aux = (x.==")") - (x.=="(")
#     pop!(aux)
#     popfirst!(aux)
#     central_comma_place = findfirst(cumsum(aux).==0) +4
#
#     # Stuff of the left tree
#     left_tree = vcat(x[2:(central_comma_place-1)],";")
#     partition = rand(1)
#     left_tree = balance_a_tree(left_tree,total_bl.*partition)
#
#     # Stuff of the right tree
#     right_tree = vcat(x[(central_comma_place+1):(length(x)-4)],";")
#     right_tree[length(right_tree)-1] = left_tree[length(left_tree)-1]
#     right_tree = balance_a_tree(right_tree,total_bl.*partition)
#
#     # Pasting the left and right trees
#     pop!(left_tree)
#     pop!(right_tree)
#     x = vcat("(",left_tree,",",right_tree,")",":",total_bl.*(Float64[1]-partition),";")
#
#     return x
# end


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


function add_random_hybrid(x, taxa,digits_bl)
    # sampling two taxa to insert the hybrid

    taxa_to_change1, taxa_to_change2 = taxa[randperm(length(taxa))[1:2]]
    # We get the branch length and the gamma parameter for the new hybrid

    println([taxa_to_change1,taxa_to_change2])

    # Finding the index of the taxa to change
    idx_node1 = findfirst(x .== taxa_to_change1)
    idx_node2 = findfirst(x .== taxa_to_change2)
    # We allways prefere the index of node 1 to be smaller
    if idx_node1 > idx_node2
        idx_node1, idx_node2 = idx_node2, idx_node1
        taxa_to_change1, taxa_to_change2 = taxa_to_change2, taxa_to_change1
    end

    #  Obs: The node1 will be added parenthesis and the node2 will be errased
    bl_leaf_parenthesis = x[idx_node1+2]
    bl_leaf_errased = x[idx_node2+2]

    gamma = round(rand(1)[1], digits = digits_bl)
    # gamma_bl_partition will divide the bl_leaf_parenthesis for the hybrid and the
    # new branch length for the node-parenthesis (this is for been ultrametric)
    gamma_bl_partition = round(rand(1)[1], digits = digits_bl)

    new_bl_leaf_parenthesis = bl_leaf_parenthesis*gamma_bl_partition
    bl_hybrid_leaf_parenthesis = bl_leaf_parenthesis*(1-gamma_bl_partition)
    bl_hybrid_leaf_errased = bl_leaf_errased - new_bl_leaf_parenthesis

    # Now we can make the changes on the array x
    # Errasing previous node
    x[idx_node2] = "#H1"
    x[idx_node2+2] = bl_hybrid_leaf_errased[1]
    insert!(x, idx_node2 + 3, 1-gamma[1])
    insert!(x, idx_node2 + 3, ":")
    insert!(x, idx_node2 + 3, ":")

    # Changes to node parenthesis
    x[idx_node1+2] = new_bl_leaf_parenthesis[1]
    insert!(x, idx_node1 + 3, gamma[1])
    insert!(x, idx_node1 + 3, ":")
    insert!(x, idx_node1 + 3, ":")
    insert!(x, idx_node1 + 3, bl_hybrid_leaf_parenthesis[1])
    insert!(x, idx_node1 + 3, ":")
    insert!(x, idx_node1 + 3, "#H1")
    insert!(x, idx_node1 + 3, ")")
    insert!(x, idx_node1, "(")
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
function Sim_tree(taxa, total_bl, is_network,digits_bl,is_unrooted)
    # Getting the topology
    x = random_tree_topology(taxa)
    # Adding some branch lengths
    x = add_random_bl(x, taxa,digits_bl)
    # Fixing the branch lengths to be "equal"
    x = balance_a_tree(x, total_bl,digits_bl,is_unrooted)
    # Adding a random hybrid
    if is_network
        x = add_random_hybrid(x, taxa,digits_bl)
    end
    return x
end

## Funcitons of njWillems2014

function Julia_willems2014(
    x,
    name_file_dist_txt,
    name_file_Input_Matrix_Distances,
    name_file_Parameters_Distances,
    name_file_new_Network_Distances,
    folder_output_distances_willems,
)
    # Folder to save the output of willems
    output_willems =
        folder_output_distances_willems * "/" * name_file_new_Network_Distances
    # Getting the netwokr from the Newick format
    net = readTopology(join(x))
    # Distance matrix
    Sigma = sharedPathMatrix(net)
    Vy = Sigma[:Tips]

    # From Cecil
    M = pairwiseTaxonDistanceMatrix(net, keepInternal = false)
    tipLabels(net)

    # Generating file to run njWillems2014
    ind = parse.(Int, tipLabels(net))
    # ind =  tipLabels(net)

    Mext = hcat(ind, M)
    n = length(tipLabels(net))


    df = DataFrame(Mext)
    df[!, :x1] = convert.(Int, df[!, :x1])

    CSV.write(name_file_dist_txt, df, writeheader = false, delim = ' ')


    # the first row with the number of tips
    f = readlines(name_file_dist_txt)
    g = open(name_file_Input_Matrix_Distances, "w")
    write(g, "$n \n")

    for l in f
        write(g, l)
        write(g, "\n")
    end
    close(g)



    # This "run" generates the file "Network_Distances.txt"
    text_to_run = `./njwillems2014 $name_file_Input_Matrix_Distances $name_file_Parameters_Distances`
    run(text_to_run)

    # We copy the file "Network_Distances.txt", that generates njwillems2014
    cp("Network_Distances.txt", output_willems, force = true)

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
