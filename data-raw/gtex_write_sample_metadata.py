import pandas as pd
import numpy as np
import h5py

sinfo = pd.read_csv("/opt/DB//GTEx/GTEx_Analysis_V6_RNA-seq/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt",
                    sep='\t',header=0,
                   encoding='ascii')

svarlist = pd.read_csv("/opt/DB//GTEx/GTEx_Analysis_V6_RNA-seq/annotations/varlist_relevant.csv",
                      sep="\t", header=0)

selectedFields = sinfo.loc[:,svarlist.loc[:,'VARNAME']]
pd.to_datetime(selectedFields.iloc[:,4])

types = []
for i in range(len(svarlist)):
    fieldName = svarlist.loc[i,'VARNAME']
    fieldType = svarlist.loc[i,'NPTYPELEN']
    if (fieldType[:8] == 'datetime'):
        types.append((fieldName, 'S10'))
    else:
        types.append((fieldName,fieldType))

sarray = np.zeros((len(selectedFields), ), dtype=types)
for field,ftype in types:
    sarray[field] = selectedFields.loc[:,field].values

dfile = h5py.File('/opt/DB/GTEx/GTEx_V6.h5', 'r+')
metagrp = dfile['metadata']


dset = metagrp.create_dataset("sample-ontoterm", (len(sarray),), np.dtype(types))
dset[...] = sarray

for i in range(len(svarlist)):
    attr_name = svarlist.loc[i,'VARNAME']
    attr_value = svarlist.loc[i, 'VARDESC']
    dset.attrs[attr_name] = attr_value


dfile.flush()
dfile.close()
