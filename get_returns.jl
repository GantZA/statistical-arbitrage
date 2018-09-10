using CSV
using Temporal
using Glob
using DataFrames

dir = "/home/gant/Documents/UCT 2018/Paper/sasa statarb/data"
files = glob("*prices.csv", dir)

price_data = CSV.read(files[1], missingstring = "NA")

function price_relatives(data)
    return data[2:end] ./ data[1:(end-1)]
end

function get_price_relatives(data)
    return_data = Matrix(size(data)[1]-1,size(data)[2]-1)
    for i in 2:size(data)[2]
        return_data[:,i-1] = price_relatives(data[:,i])
    end

    return_data = DataFrame(hcat(data[2:end,1], return_data))
    names!(return_data, names(data))
    return return_data
end

return_data = get_price_relatives(price_data)
CSV.write("JSE_top40_2003_2018.csv",return_data)
