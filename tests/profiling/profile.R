library(rglyvis)

# Rprof('prof-getScatterData.out', interval=0.0002)
# getScatterData(x = "ENSG00000176022",y="ENSG00000027847", sampleGroups="Brain", sampleGrouping="SMTS", sampleMetaFields=c("SMTS"), db="gtex", processing="toil-rsem", unit="tpm", read.from.redis=TRUE)
# Rprof(NULL)
# sink('sum-getScatterData.out')
# summaryRprof('prof-getScatterData.out')
# sink()

prof_and_sum <- function(output, func, ...)  {
    prof_out = paste0('prof-', output, '.out')
    sum_out  = paste0('sum-', output, '.out')
    Rprof(prof_out, interval=0.0002)
    func(...)
    Rprof(NULL)
    sink(sum_out)
    print(summaryRprof(prof_out))
    sink()
}

start=proc.time()
# prof_and_sum('getScatterData', getScatterData,x = "ENSG00000176022",y="ENSG00000027847", sampleGroups="Brain", sampleGrouping="SMTS", sampleMetaFields=c("SMTS", "AGE", "SMTSD", "GENDER", "RACE", "ONTOTERM"), db="gtex", processing="toil-rsem", unit="tpm", read.from.redis=TRUE)
prof_and_sum('getGeneExpressionMatrix',
             getGeneExpressionMatrix,
             genes = c("ENSG00000176022", "ENSG00000103489","ENSG00000015532","ENSG00000027847"),
             sampleGroups=c("Brain"),
             sampleGrouping="SMTS",
             sampleMetaFields=c("SMTS", "SMTSD", "AGE","GENDER"),
             db = "gtex", processing = "toil-rsem", unit="tpm", expect='json', read.from.redis = TRUE)
duration=proc.time() - start
