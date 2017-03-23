# GTEx data survey
options(stringsAsFactors = FALSE)

## sample attributes
sattrs = read.csv('/opt/DB/GTEx/GTEx_Analysis_V6_RNA-seq/annotations/varlist.csv', sep="\t", stringsAsFactors=FALSE)
sampleInfo = readRDS("/opt/DB/GTEx/GTEx_Analysis_V6_RNA-seq/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.RDS")
print(names(sampleInfo))
print(sattrs[['VARNAME']])
setdiff(names(sampleInfo), sattrs[['VARNAME']])
setdiff(sattrs[['VARNAME']], names(sampleInfo))

## subject phenotypes
datadir = "/home/trang/ncbi/dbGaP-10719/53961/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/PhenotypeFiles/decrypted/"
subjectInfo = data.table::fread(paste0(datadir,"phs000424.v6.pht002742.v6.p1.c1.GTEx_Subject_Phenotypes.GRU.txt"), header=T, sep="\t",skip = 10)
mapSubjectSample = data.table::fread(paste0(datadir,"phs000424.v6.pht002741.v6.p1.GTEx_Sample.MULTI.txt"),header = T, sep="\t", skip = 10)

### Replace encoded values by friendly values for some fields of interest
replaceEncoding = function(x, column.name, labels) {
    tmp = factor(x[[column.name]], labels=labels)
    x[,column.name] = as.character(tmp)
    return(x)
}

subjectInfo = replaceEncoding(subjectInfo, 'GENDER', labels=c("Male", "Female"))
subjectInfo = replaceEncoding(subjectInfo, 'RACE',  labels=c("Asian", "Black or African American", "White", "American Indian or Alaska Native", "Not Reported", "Unknown"))
subjectInfo = replaceEncoding(subjectInfo, 'ETHNCTY',  c("Not Hispanic or Latino", "Hispanic or Latino", "CEPH","Not Reported", "Unknown"))

library(rhdf5)
hf = rhdf5::H5Fopen('/var/www/html/data/GTEx_V6.h5')
rhdf5::h5writeDataset(subjectInfo,hf,'/metadata/subject',DataFrameAsCompound=TRUE)
H5Fflush(hf)

rhdf5::h5writeDataset(mapSubjectSample, hf, "/metadata/mapSubjectSample", DataFrameAsCompound=TRUE)
H5Fflush(hf)

H5Fclose(hf)


## Reversing logarithm on the data downloaded from Xena HUB https://xenabrowser.net/datapages/?cohort=GTEX
### FPKM data set
hf = rhdf5::H5Fopen('/var/www/html/data/GTEx_V6.h5')
afterLog = rhdf5::h5read(hf, '/toil-rsem/gene/fpkm-log2')
beforeLog = (2^afterLog) - 0.001
# #### Check for negative values
# negs = beforeLog[beforeLog<0]
# for (tol in c(1e-7, 1.5e-7, 1e-8) ) {
#     print(paste("The number   of values <", -tol, ":", sum(negs< -tol)))
#     print(paste("The fraction of values <", -tol, ":", sum(negs< -tol) / length(negs)*1.0))
# }
### Since all of the negative values are -1.0907E-8, it is safe to zero out all of them
beforeLog[beforeLog<0] = 0
rhdf5::h5writeDataset(beforeLog, hf, '/toil-rsem/gene/fpkm')
H5Fflush(hf)

### TPM data set
afterLog = rhdf5::h5read(hf, '/toil-rsem/gene/tpm-log2')
beforeLog = (2^afterLog) - 0.001
rhdf5::h5writeDataset(beforeLog, hf, '/toil-rsem/gene/tpm')
beforeLog[beforeLog<0] = 0
rhdf5::h5writeDataset(beforeLog, hf, '/toil-rsem/gene/tpm')
H5Fflush(hf)

H5Fclose(hf)


## Why the expression levels in TOIL-processed data seems to be assigned into discrete bins?
## Such visually discrete data points should be expected for low expression levels.
library(rglyvis)
genes = c('ENSG00000002587', 'ENSG00000182601')   # HS3ST1, HS3ST4
# mj = rglyvis::getGeneExpressionMatrix(genes=genes, sampleGroups=c("Nerve"), expect='json')
# write(prettify(mj),file='test.json')
m1 = getGeneExpressionMatrix(genes=genes, sampleGroups=c("Nerve"), expect='datatable')
m2 = getGeneExpressionMatrix(genes=genes, sampleGroups=c("Blood"), expect='datatable')
m3 = getGeneExpressionMatrix(genes=genes, sampleGroups=c("Liver"), expect='datatable')
paste("Nerve tissue has", nrow(m1), "samples")
count1 = table(m1$ENSG00000002587)
count2 = table(m2$ENSG00000002587)
same12 = intersect(m1$ENSG00000002587, m2$ENSG00000002587)
same13 = intersect(m1$ENSG00000002587, m3$ENSG00000002587)
same23 = intersect(m2$ENSG00000002587, m3$ENSG00000002587)
same123 = intersect(same12, m3$ENSG00000002587)


sameAB1= intersect(m1$ENSG00000002587, m1$ENSG00000182601)
sameAB3= intersect(m3$ENSG00000002587, m3$ENSG00000182601)
