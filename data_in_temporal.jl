using Temporal
using Glob
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
files = readdir(dir)
files_name = string(dir, "/",files[3])
JSE = tsread(files_name)
JSE[1:5]
JSE.values = reduce_outliers.(JSE.values)

b_JSE1, q_JSE1, PLt_JSE1, PLnt_JSE1, H_JSE1 = simulation(false, 1, 5823, JSE, [5,5,5], [1,3], 0.01, "Uni")
b_JSE2, q_JSE2, PLt_JSE2, PLnt_JSE2, H_JSE2 = simulation(false, 1, 5823, JSE, [5,5,5], [1,3], 0.01, "Uni EWMA")
b_JSE, q_JSE, PLt_JSE, PLnt_JSE, H_JSE = simulation(false, 500, 1500, JSE, [5,5,5], [1,3], 0.01, "Uni")

b_JSE1, q_JSE1, PLt_JSE1, PLnt_JSE1, H_JSE1 = simulation(false, 1, 4000, JSE, [5,5,5], [1,3], 0.01, "Uni")

b_JSE2, q_JSE2, PLt_JSE2, PLnt_JSE2, H_JSE2 = simulation(false, 1, 1000, JSE, [5,5,5], [1,3], 0.01, "Uni")
b_JSE3, q_JSE3, PLt_JSE3, PLnt_JSE3, H_JSE3 = simulation(false, 1000, 2000, JSE, [5,5,5], [1,3], 0.01, "Uni")
b_JSE4, q_JSE4, PLt_JSE4, PLnt_JSE4, H_JSE4 = simulation(false, 2000, 3000, JSE, [5,5,5], [1,3], 0.01, "Uni")
b_JSE5, q_JSE5, PLt_JSE5, PLnt_JSE5, H_JSE5 = simulation(false, 3000, 4000, JSE, [5,5,5], [1,3], 0.01, "Uni")


using Plots
plotly()

plot(PLt_JSE1)
plot(PLt_JSE2)
