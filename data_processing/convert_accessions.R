setwd("/home/xabier/Dropbox/induction")

chipseq_datasets <- read.table("chipseq_datasets.csv", header = TRUE, sep = ",", fill = TRUE)

gsm_ids <- unlist(chipseq_datasets)
gsm_ids <- unlist(strsplit(gsm_ids, "\t"))


gsm_ids <- gsm_ids[grepl("GSM", gsm_ids)]



gsm_ids <- gsm_ids[!is.na(gsm_ids)]

gsm_ids<-as.vector(gsm_ids)

library(GEOquery)


srx_accessions <- sapply ( gsm_ids, function(gsm_id) {
    geo_data <- getGEO(gsm_id, destdir = tempdir())
    srx_id <- Meta(geo_data)$relation[grepl("SRX", Meta(geo_data)$relation)]
    srx_id <- sub(".*(SRX[0-9]+).*", "\\1", srx_id)
    return(srx_id)
} )




write.table(srx_accessions, file = "~/srx_accessions.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


srx_accessions_df <- data.frame(GSM = gsm_ids, SRX = srx_accessions)

write.table(srx_accessions_df, file = "~/srx_accessions_df.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)


system("
    export_shiva srx\*txt proneural/chip-seq
")