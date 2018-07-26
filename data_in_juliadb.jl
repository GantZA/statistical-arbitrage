using JuliaDB

dir = "/Documents/Paper/new_sasa_paper/data"
files = glob("*.csv",dir)
JSE_price_hist_db = loadtable(files, delim = ';')

JSE_price_hist_db
