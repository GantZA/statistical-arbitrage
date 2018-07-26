##### Testing environment
using TimeSeries

using Temporal

dir = "/Documents/Paper/new_sasa_paper/data"
files = glob("*.csv",dir)
JSE_price_hist = loadtable(files, delim = ';')
files

using CSV

temp_file = CSV.read(files[1]; delim = ';')

temp_file[:,1][1]

temp_file[:,1] = Date.(temp_file[:,1], "yyyy/mm/dd")

temp_file
CSV.write("/Documents/Paper/new_sasa_paper/data/JSE_TOP_19942017_002.csv", temp_file)
