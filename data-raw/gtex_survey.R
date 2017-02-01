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
