#dbpath = list(gtex="/opt/DB/GTEx/GTEx_V6.h5")
dbpath = list(gtex="/var/www/html/data/GTEx_V6.h5")
container = new.env(parent=emptyenv())

devtools::use_data(dbpath, container, internal=TRUE, overwrite=TRUE)
