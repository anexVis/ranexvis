# GTEx data survey
options(stringsAsFactors = FALSE)
sattrs = read.csv('/opt/DB/GTEx/GTEx_Analysis_V6_RNA-seq/annotations/varlist.csv', sep="\t", stringsAsFactors=FALSE)
sampleInfo = readRDS("/opt/DB/GTEx/GTEx_Analysis_V6_RNA-seq/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.RDS")
print(names(sampleInfo))
print(sattrs[['VARNAME']])
setdiff(names(sampleInfo), sattrs[['VARNAME']])
setdiff(sattrs[['VARNAME']], names(sampleInfo))
