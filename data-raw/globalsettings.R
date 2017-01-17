hgnc2ensembl <- function(hgnc_symbol) {
    # the var name will be made col name in data.table. keep it the same to merge
    ensembl = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                               dataset = 'hsapiens_gene_ensembl',
                               host='www.ensembl.org')
    rep_ensembl = biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                        filter='hgnc_symbol', values=hgnc_symbol, mart=ensembl)
    ids = data.table::data.table(hgnc_symbol)
    names(ids) = c('hgnc_symbol')
    tmp = merge(ids, rep_ensembl, by='hgnc_symbol', sort=FALSE)
    return(tmp[['ensembl_gene_id']])
}


loadGeneSets = function() {
	output = list()
    mapping <- read.table("data-raw/geneSets/labels", sep='\t', header=T)
    for (i in 1:nrow(mapping)) {
        name = sub(".txt", "", mapping[i,1])
        geneNames  = readLines(paste("data-raw/geneSets/", mapping[i,2], sep=""))
        ensemblIds  = hgnc2ensembl(geneNames)
        output[[i]] = list('name'= name, 'value'= ensemblIds)
    }
    return(output)
}

geneSets = loadGeneSets()
# dbpath = list(gtex="/opt/DB/GTEx/GTEx_V6.h5")
dbpath = list(gtex="/var/www/html/data/GTEx_V6.h5")
container = new.env(parent=emptyenv())

devtools::use_data(dbpath, container, geneSets, internal=TRUE, overwrite=TRUE)



