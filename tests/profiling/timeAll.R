library(rglyvis)
library(microbenchmark)

selectedIds = getGeneSets(expect='dt')[[2]]$value
symbols = c('XYLT1', 'XYLT2','B4GALT7', 'B3GALT6',
                'B3GAT3', 'EXTL2', 'EXTL3', 'EXTL1',
                'EXT1', 'EXT2', 'NDST1', 'NDST2',
                'NDST3', 'NDST4', 'GLCE', 'HS2ST1',
                'HS6ST1', 'HS6ST2', 'HS6ST3', 'HS3ST1',
                'HS3ST2', 'HS3ST3A1', 'HS3ST3B1', 'HS3ST4',
                'HS3ST5', 'HS3ST6')
print("Timing query functions")
microbenchmark(geneList.1= getGeneList("gtex", c("EnsemblID", "HGNC"), expect="df", read.from.redis=TRUE),
               geneList.2 =getGeneList("gtex", c("EnsemblID", "HGNC"), withEnsemblVersion=FALSE,expect="df", read.from.redis = TRUE),
               sampleMetadata = getSampleMetadata("gtex", cols=c('SAMPID', "SMTS", "SMTSD", "AGE", "GENDER", "RACE"), expect="dt"),
               sampleGroupingList = getSampleGroupingList(expect='df',grouping='SMTS'),
               scatter.json = getScatterData(x = "ENSG00000176022.3", y = "ENSG00000027847.9", sampleGroups="Bladder", sampleGrouping="SMTSD", sampleMetaFields = c("SMTS", "SMTSD", "GENDER", "RACE"), db = "gtex", processing = "toil-rsem", unit="tpm", read.from.redis=TRUE),
               expr.dtable.10cols =  getGeneExpressionMatrix(genes = "ENSG00000176022.3", sampleGroups=c("Bladder"), sampleGrouping="SMTS", sampleMetaFields=c("SMTS", "SMTSD"), db = "gtex", processing = "toil-rsem", unit="tpm", expect='datatable', read.from.redis = TRUE),
               expr.dtable.1000cols =  getGeneExpressionMatrix(genes = "ENSG00000176022.3", sampleGroups=c("Brain"), sampleGrouping="SMTS", sampleMetaFields=c("SMTS", "SMTSD"), db = "gtex", processing = "toil-rsem", unit="tpm", expect='datatable', read.from.redis = TRUE),
               expr.dtable.30rows =  getGeneExpressionMatrix(genes = selectedIds, sampleGroups=c("Brain"), sampleGrouping="SMTS", sampleMetaFields=c("SMTS", "SMTSD"), db = "gtex", processing = "toil-rsem", unit="tpm", expect='datatable', read.from.redis = TRUE),
               expr.json =  getGeneExpressionMatrix(genes = "ENSG00000176022.3", sampleGroups=c("Bladder"), sampleGrouping="SMTS", sampleMetaFields=c("SMTS", "SMTSD"), db = "gtex", processing = "toil-rsem", unit="tpm", expect='json', read.from.redis = TRUE),
               times=10)

print("Timing Utility Functions")
microbenchmark(ensembl2hgnc = rglyvis:::ensembl2hgnc(selectedIds),
               hgnc2ensembl = rglyvis:::hgnc2ensembl(symbols),
               times = 10)
