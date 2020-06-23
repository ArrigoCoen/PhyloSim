
## Test

using PhyloNetworks
using PhyloPlots


net = readTopology("((a,d),(b,c));")

plot(net,:R)



## End of test

"(a:2,(b:4,c:3):23);"
Pkg.status()

## Error

Pkg.add("Distributions")




## Not error
Pkg.rm("RecipesBase")



##
