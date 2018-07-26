using Temporal
using JuliaDB
include("strategies.jl")
include("simulation.jl")
include("learning_algorithms.jl")

function reduce_outliers(entry)
    if (entry >= 1.3) | (entry <= 0.7)
        return 1
    else
        return entry
    end
end

dir = "/Documents/Paper/new_sasa_paper/data"
files = glob("*.csv",dir)
files[2]
JSE = tsread(files[2])
JSE[1:5]
JSE.values = reduce_outliers.(JSE.values)
a = online_zbcrp(JSE[1:5].values, 2, false,  false)
d = online_anticor(JSE[1:5].values, 2, 3,  ,  false)
simulation(false, 1, 10, JSE, [5,5,5], [1,3], 1, "Uni")

JSE[1:5].values * ones(36)
