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
