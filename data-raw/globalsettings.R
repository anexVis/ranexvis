#dbpath = list(gtex="/opt/DB/GTEx/GTEx_V6.h5")
dbpath = list(gtex="/var/www/html/data/GTEx_V6.h5")

devtools::use_data(dbpath, internal=TRUE, overwrite=TRUE)
