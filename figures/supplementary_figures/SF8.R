################## gene ontologies ################
###################################################
####### Supplementary figure 8




# Subdivide fc into mouse studies: select studies that are not human
fc_mouse <- fc[!grepl("hESC|IMEC|LNCaP|U2OS|SH-SY5Y|G523NS|GI-MEN|Human fibroblast|MRC5", fc$study), ]

fc_mouse <- fc_mouse %>% dplyr::select(gene_name, log2FoldChange, study)

fc_human <- fc[grep("hESC|IMEC|LNCaP|U2OS|SH-SY5Y|G523NS|GI-MEN|Human fibroblast|MRC5", fc$study), ]
fc_human <- fc_human %>% dplyr::select(gene_name, log2FoldChange, study)

library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)



# Helper function for GO enrichment per study
run_go_enrichment <- function(df, species_db, species_name) {
    studies <- unique(df$study)
    go_results <- list()
    for (study in studies) {
        cat("Running GO enrichment for study:", study, "\n")
        subdf <- df[df$study == study & !is.na(df$log2FoldChange), ]
        genes <- unique(subdf$gene_name)
        gene_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = species_db)
        if (is.null(gene_ids) || nrow(gene_ids) == 0) next
        gene_list <- subdf$log2FoldChange[match(gene_ids$SYMBOL, subdf$gene_name)]
        names(gene_list) <- gene_ids$ENTREZID
        gene_list <- sort(gene_list, decreasing = TRUE)
        ego <- tryCatch({
            gseGO(
                geneList = gene_list,
                OrgDb = species_db,
                ont = "BP",
                keyType = "ENTREZID",
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 1,
                verbose = FALSE
            )
        }, error = function(e) NULL)
        go_results[[study]] <- ego
    }
    return(go_results)
}

# Human studies GO enrichment
go_human <- run_go_enrichment(fc_human, org.Hs.eg.db, "human")

saveRDS(go_human, file = "~/proneural/go_human_gsea.rds")

# Mouse studies GO enrichment
fc_mouse<- fc_mouse %>% filter(!is.na(study))
go_mouse <- run_go_enrichment(fc_mouse, org.Mm.eg.db, "mouse")
saveRDS(go_mouse, file = "~/proneural/go_mouse_gsea.rds")




#### BALLOON PLOT ########################################################
############################################# 



############## mouse + human

go_mouse <- readRDS("~/proneural/go_mouse_gsea.rds")
go_human <- readRDS("~/proneural/go_human_gsea.rds")

## join mouse and human gene ontology objects. go_mouse and go_human objects
go_all <- c(go_mouse, go_human)


## get top ontologies
ontologies_all <- unique(unlist(
    lapply(go_all, function(res) {
        if (is.null(res) || nrow(res@result) == 0) return(NULL)
        head(res@result$Description, 10)
    })
))

### dataframe. one column description, the other GO:ID

df_ontologies<-data.frame(
    Description = ontologies_all,
    GO_ID = unlist(lapply(ontologies_all, function(desc) {
        # Find the first GO ID for this description in any result
        for (res in go_all) {
            if (!is.null(res) && nrow(res@result) > 0) {
                idx <- which(res@result$Description == desc)
                if (length(idx) > 0) return(res@result$ID[idx[1]])
            }
        }
        return(NA)
    })),
    stringsAsFactors = FALSE
)

write.table(df_ontologies, file="~/ontology_balloon_all_gsea_descriptions_ids.txt", sep="\t", row.names=FALSE, quote=FALSE)

## i ordered outside the GO:IDS. now read them
ordered_ids<-read.table("~/tips_order.txt")$V1

ordered_descriptions<-df_ontologies$Description[match(ordered_ids, df_ontologies$GO_ID)]


balloon_data_all <- do.call(rbind, lapply(names(go_all), function(study) {
    res <- go_all[[study]]
    if (is.null(res) || nrow(res@result) == 0) return(NULL)
    df <- res@result
    data.frame(
        Study = study,
        Ontology = df$Description,
        NES = df$NES,
        log10padj = -log10(df$p.adjust)
    )
}))
balloon_data_all <- balloon_data_all %>% filter(Ontology %in% ontologies_all)

### cluster ontologies by NES
# Prepare matrix: rows = ontologies, columns = studies, values = NES
ontology_matrix_balloon_all <- balloon_data_all %>%
    dplyr::select(Study, Ontology, NES) %>%
    pivot_wider(names_from = Study, values_from = NES)

ontology_matrix_balloon_all <- as.data.frame(ontology_matrix_balloon_all)
rownames(ontology_matrix_balloon_all) <- ontology_matrix_balloon_all$Ontology
ontology_matrix_balloon_all$Ontology <- NULL

# Cluster rows and columns using hierarchical clustering
row_dist_all <- dist(ontology_matrix_balloon_all, method = "euclidean")
row_clust_all <- hclust(row_dist_all)
row_order_all <- row_clust_all$order

# Get ordered ontology and study names
ordered_ontologies_all <- rownames(ontology_matrix_balloon_all)[row_order_all]

balloon_data_all$Ontology <- factor(balloon_data_all$Ontology, levels = ordered_ontologies_all)
# Plot dendrogram of the ordered ontologies

dend_all <- as.dendrogram(row_clust_all)


pdf("~/Dropbox/induction/7_expression/ontology_dendrogram_all_gsea.pdf", width = 12, height = 30)
par(mar = c(2, 3, 2, 30)) # Menos margen a la izquierda (segundo valor)
plot(dend_all, main = "Dendrogram of GO BP Ontologies (all studies)", horiz = TRUE, cex = 1.2)
dev.off()



orden_all <- readRDS("~/proneural/chip-seq/orden.RDS")
orden_all <- orden_all[orden_all %in% unique(balloon_data_all$Study)]
balloon_data_all$Study <- factor(balloon_data_all$Study, levels = orden_all)




# Plot balloon plot
pdf("~/Dropbox/induction/7_expression/ontology_balloon_all_gsea.pdf", width = 14, height = 40)
ggplot(balloon_data_all, aes(x = Study, y = Ontology)) +
    geom_point(aes(size = log10padj, color = NES)) +
    scale_size_continuous(name = "-log10 adj p") +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "NES") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "GO BP enrichment in all studies (balloon plot)",
         x = "Study", y = "Ontology")
dev.off()



