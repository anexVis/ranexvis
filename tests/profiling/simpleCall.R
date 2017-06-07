start=proc.time()

expr_json = getScatterData(x = "ENSG00000176022",
                               y = "ENSG00000176022",
                               sampleGroups="Bladder",
                               sampleGrouping="SMTSD",
                               sampleMetaFields = c("SMTSD"),
                               db = "gtex",
                               processing = "toil-rsem",
                               unit="tpm",
                               read.from.redis=TRUE
    )

duration = proc.time() - start

