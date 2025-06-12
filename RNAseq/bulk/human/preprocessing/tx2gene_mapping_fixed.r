#Set file paths
tx2gene_path <- "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_reference_genome/GTF_file/tx2gene_gencode_v38.csv"

#quant_file <- "/path/to/any/sample/quant.sf"  #Use any one quant.sf file for checking
quant_file <- "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_16hrs_batch3_27022025/salmon_output_batch3/E6R2_16hrs/quant.sf"

#Load files
tx2gene <- read.csv(tx2gene_path, stringsAsFactors = FALSE)
quant <- read.table(quant_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Strip version numbers from both if needed
tx2gene$TXNAME_noversion <- sub("\\..*", "", tx2gene$TXNAME)
quant$Name_noversion <- sub("\\..*", "", quant$Name)

#Check overlap
cat("Sample transcript IDs from quant.sf:\n")
print(head(quant$Name))
cat("\nSample transcript IDs from tx2gene:\n")
print(head(tx2gene$TXNAME))

cat("\nChecking overlap WITH version numbers:\n")
cat("Number of matching IDs: ", sum(quant$Name %in% tx2gene$TXNAME), "out of", nrow(quant), "\n")

cat("\nChecking overlap WITHOUT version numbers:\n")
cat("Number of matching IDs: ", sum(quant$Name_noversion %in% tx2gene$TXNAME_noversion), "out of", nrow(quant), "\n")

#If more matches WITHOUT version, save a version-stripped tx2gene
if (sum(quant$Name_noversion %in% tx2gene$TXNAME_noversion) > sum(quant$Name %in% tx2gene$TXNAME)) {
  cat("\n*** More matches without version numbers! Creating version-stripped tx2gene file... ***\n")
  tx2gene_noversion <- tx2gene
  tx2gene_noversion$TXNAME <- tx2gene_noversion$TXNAME_noversion
  tx2gene_noversion$GENEID <- sub("\\..*", "", tx2gene_noversion$GENEID)
  tx2gene_noversion <- tx2gene_noversion[, c("TXNAME", "GENEID")]
  write.csv(tx2gene_noversion, "tx2gene_noversion.csv", row.names = FALSE)
  cat("Saved: tx2gene_noversion.csv\n")
} else {
  cat("\n*** Your versioned tx2gene mapping is fine. ***\n")
}
