using Temporal

dir = "/Documents/Paper/new_sasa_paper/data"
files = glob("*.csv",dir)
files[2]
JSE = tsread(files[2])


JSE[2]
JSE[end-100:end,1:4]

JSE[[:AGL, :ANG]]

JSE["2015/2016", :AGL]

function reduce_outliers(entry)
    if (entry >= 1.3) | (entry <= 0.7)
        entry = 1
    end
end

JSE[1] = 1
JSE[1]
