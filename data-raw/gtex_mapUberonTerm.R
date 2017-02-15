library(ontoCAT)
options(java.parameter='-Xms2g -Xms2g')

oldFile = "/opt/DB/GTEx/GTEx_Analysis_V6_RNA-seq/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
newFile = "/opt/DB/GTEx/GTEx_Analysis_V6_RNA-seq/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS-fullMapping.txt"
sinfo = data.table::fread(oldFile,stringsAsFactors = FALSE,colClasses = "character")

uberonAnn = sinfo[['SMUBRID']]
uberon = getOntologyNoReasoning('/opt/DB/Uberon/uberon.owl')
efo = getOntologyNoReasoning('/opt/DB/EFO/efo.owl')

# The EFO_0000496 annotates fibroblast, which should be CL:0000057 in EFO
ontologyTerms = sapply(uberonAnn, function(x) {
    if (x=="") return("")
    else if (grepl("^\\d+", x)) return(getTermNameById(uberon,paste0("UBERON:",x)))
    else if (grepl("^EFO_0000496",x)) return(getTermNameById(efo,"CL:0000057"))
    else if (grepl("^EFO_\\d+",x)) return(getTermNameById(efo,x))
    else return(x)
})

sinfo$ONTOTERM <- ontologyTerms
write.table(sinfo,newFile,sep = "\t",quote = FALSE,row.names = FALSE)
# sinfo_new  = data.table::fread(newFile,stringsAsFactors = FALSE, colClasses = "character")
