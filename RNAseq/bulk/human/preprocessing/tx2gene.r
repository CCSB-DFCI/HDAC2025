library(GenomicFeatures)
txdb <- makeTxDbFromGFF("gencode.v38.annotation.gtf", format = "gtf")
tx2gene <- select(txdb, keys(txdb, "TXNAME"), "GENEID", "TXNAME")
write.csv(tx2gene, "tx2gene_gencode_v38.csv", row.names = FALSE)
