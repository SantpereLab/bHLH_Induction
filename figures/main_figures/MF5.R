
##################################################################
##### fold change accessibility ~ motifs  #######
##################################################################
####### Main figure 5a and supplementary figure 7a


######### initial processing of the data
######################################

### ebox counts with flanking regions
peak_and_flanking_counts<- readRDS("~/peak_and_flanking_counts.RDS")

##  to get peakid new --->> peakid
canonical_eboxes <- readRDS("~/proneural/chip-seq/combined_canonical_eboxes.RDS")
canonical_eboxes$peakid_new <- paste0(canonical_eboxes$chrom,"-",canonical_eboxes$summit)

peakidnew_peakid<- canonical_eboxes %>% select(peak_id, peakid_new,study) %>% distinct()

## rename study
study_mapping<-readRDS("~/proneural/chip-seq/study_mapping.RDS")
peakidnew_peakid$study<- study_mapping[match(peakidnew_peakid$study, study_mapping$study), "renamed_study"]

######### chromatin fold change data
combined_foldchange_data_join<-readRDS("~/proneural/foldchange_accessibility_with_tss_distance.RDS")

accessibility_info_by_peak <- combined_foldchange_data_join %>% select(peakid_new, study, fold_change, distanceToTSS) %>% distinct()

## rename study
accessibility_info_by_peak$study <- study_mapping[match(accessibility_info_by_peak$study, study_mapping$study), "renamed_study"]

## add peak_id column
accessibility_info_by_peak <- accessibility_info_by_peak %>%
  left_join(peakidnew_peakid, by = c("peakid_new", "study")) %>% dplyr::filter(!is.na(peak_id)) %>% distinct()


## add accessibility_info_by_peak info to peak_and_flanking_counts
peak_and_flanking_counts_with_accessibility <- peak_and_flanking_counts %>%
  dplyr::left_join(accessibility_info_by_peak, by = c("peak_id", "study")) %>%
  dplyr::filter(!is.na(fold_change)) %>% distinct()


saveRDS(peak_and_flanking_counts_with_accessibility, "~/peak_and_flanking_counts_with_accessibility.RDS")

peak_and_flanking_counts_with_accessibility<- readRDS("~/peak_and_flanking_counts_with_accessibility.RDS")

peak_and_flanking_counts_with_accessibility$log2FoldChange<- log2(peak_and_flanking_counts_with_accessibility$fold_change + 1e-6)

##################################################################
##### fold change accessibility ~ motifs in promoters  ###########
##################################################################

peak_and_flanking_counts_with_accessibility_promoters<- peak_and_flanking_counts_with_accessibility %>%
  dplyr::filter( abs (distanceToTSS) <= 2000 )


### filter the studies that have at least 1000 peaks
peaks_per_study_acc <- peak_and_flanking_counts_with_accessibility_promoters %>% select(study, peak_id) %>% distinct() %>%
  group_by(study) %>%
  summarise(n_peaks = n(), .groups = "drop")

peaks_per_study_acc %>% arrange(n_peaks) %>% print(n=Inf, width=Inf)

peak_and_flanking_counts_with_accessibility_promoters_filtered<- peak_and_flanking_counts_with_accessibility_promoters %>%
  dplyr::filter(study %in% (peaks_per_study_acc %>% filter( n_peaks >= 1000) %>% pull(study)) )







# Linear models with accessibility log2FoldChange in promoters
  # (identical in construction to the expression models)
  library(dplyr)
  library(tidyr)
  library(broom)

  motif_levels <- c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAC-CAC",
                    "CAT-CAC","CAA-CAT","CAA-CAG","CAA-CAC","CAA-CAA")

  df_prom <- peak_and_flanking_counts_with_accessibility_promoters_filtered %>%
    select(study, peak_id, central_dinucleotide, motif_enrichment, log2FoldChange,
           dplyr::any_of("log_accessibility")) %>%
    distinct() %>%
    mutate(central_dinucleotide = factor(central_dinucleotide, levels = motif_levels)) %>%
    pivot_wider(
      names_from = central_dinucleotide,
      values_from = motif_enrichment,
      names_prefix = "motif_enrichment_"
    )

lm_fits <- lapply(split(df_prom, df_prom$study), function(df_study) {
    curr_study <- unique(df_study$study)
    
    # 1. Identify predictors
    preds <- grep("^motif_enrichment_", colnames(df_study), value = TRUE)
    
    # 2. CRITICAL CLEANING: Keep only finite and real numerical values
    df_study_clean <- df_study %>% 
        dplyr::filter(is.finite(log2FoldChange))
    
    # 3. Variance validation in predictors
    preds_validos <- preds[vapply(df_study_clean[preds], function(x) {
        length(unique(x[!is.na(x)])) > 1
    }, logical(1))]
    
    # Safety filters
    if (length(preds_validos) == 0) return(NULL)
    if (nrow(df_study_clean) <= (length(preds_validos) + 1)) return(NULL)

    # 4. Imputation of NAs in predictors
    df_study_clean[preds_validos] <- lapply(df_study_clean[preds_validos], 
                                            function(x) replace(x, is.na(x), 0))

    # 5. Formula construction
    preds_escaped <- paste0("`", preds_validos, "`")
    fm <- as.formula(paste("log2FoldChange ~", paste(preds_escaped, collapse = " + ")))
    
    # 6. Model execution
    fit <- tryCatch(lm(fm, data = df_study_clean), error = function(e) NULL)
    
    if (is.null(fit)) return(NULL)

    # 7. Result formatting
    broom::tidy(fit) %>%
        dplyr::filter(term != "(Intercept)") %>%
        dplyr::mutate(
            term = gsub("`", "", term),
            study = curr_study,
            central_dinucleotide = sub("^motif_enrichment_", "", term)
        )
}) %>% dplyr::bind_rows()

  print(lm_fits)
  lm_prom_fits <- lm_fits %>%
    dplyr::mutate(q_value = p.adjust(p.value, method = "BH"))
  saveRDS(lm_prom_fits, "~/accessibility_fc~motif_enrichment_lm_promoters_fits.RDS")

  # LM in promoters controlling for initial accessibility
  fits_list <- lapply(split(df_prom, df_prom$study), function(df_study) {
    if (!"log_accessibility" %in% colnames(df_study)) return(NULL)
    curr_study <- unique(df_study$study)

    preds <- grep("^motif_enrichment_", colnames(df_study), value = TRUE)

    df_clean <- df_study %>%
      dplyr::filter(is.finite(log2FoldChange), is.finite(log_accessibility))

    if (nrow(df_clean) == 0 || length(unique(df_clean$log_accessibility)) <= 1) return(NULL)

    preds_valid <- preds[vapply(df_clean[preds], function(x) {
      length(unique(x[!is.na(x)])) > 1
    }, logical(1))]

    if (length(preds_valid) == 0) return(NULL)
    if (nrow(df_clean) <= (length(preds_valid) + 1)) return(NULL)

    df_clean[preds_valid] <- lapply(df_clean[preds_valid], function(x) replace(x, is.na(x), 0))

    fm <- as.formula(
      paste(
        "log2FoldChange ~",
        paste(c(paste0("`", preds_valid, "`"), "`log_accessibility`"), collapse = " + ")
      )
    )

    fit <- tryCatch(lm(fm, data = df_clean), error = function(e) NULL)
    if (is.null(fit)) return(NULL)

    broom::tidy(fit) %>%
      dplyr::filter(term != "(Intercept)") %>%
      dplyr::mutate(
        term = gsub("`", "", term),
        study = curr_study,
        central_dinucleotide = sub("^motif_enrichment_", "", term)
      )
  })

  fits_list <- Filter(Negate(is.null), fits_list)
  lm_prom_acc_fits <- if (length(fits_list)) dplyr::bind_rows(fits_list) else NULL

  if (!is.null(lm_prom_acc_fits)) {
    lm_prom_acc_fits <- lm_prom_acc_fits %>%
      dplyr::mutate(q_value = p.adjust(p.value, method = "BH"))
    saveRDS(lm_prom_acc_fits, "~/accessibility_fc~motif_enrichment_plus_initial_accessibility_lm_promoters_fits.RDS")
  }

  # Ensure 'slope' column for heatmaps
  if ("estimate" %in% colnames(lm_prom_fits) && !("slope" %in% colnames(lm_prom_fits))) {
    lm_prom_fits <- dplyr::mutate(lm_prom_fits, slope = estimate)
  }
  if (!is.null(lm_prom_acc_fits) &&
      ("estimate" %in% colnames(lm_prom_acc_fits)) &&
      !("slope" %in% colnames(lm_prom_acc_fits))) {
    lm_prom_acc_fits <- dplyr::mutate(lm_prom_acc_fits, slope = estimate)
  }


  # Compact version of the heatmap for promoters (smaller PDF and dynamic according to the number of studies)
  plot_model_heatmap_small <- function(df, title, outfile) {
    df <- df %>%
      dplyr::mutate(
        study = factor(study, levels = orden),
        central_dinucleotide = factor(central_dinucleotide, levels = motif_levels)
      )

    mat <- df %>%
      dplyr::select(study, central_dinucleotide, slope) %>%
      tidyr::pivot_wider(names_from = central_dinucleotide, values_from = slope) %>%
      as.data.frame()
    rownames(mat) <- mat$study
    mat <- as.matrix(mat[, -1])

    rows <- intersect(orden, rownames(mat))
    cols <- intersect(motif_levels, colnames(mat))
    if (length(rows) == 0 || length(cols) == 0) {
      warning("No rows or columns to plot for: ", title)
      return(invisible(NULL))
    }
    mat <- mat[rows, cols, drop = FALSE]

    sig <- df %>%
      dplyr::mutate(sig = ifelse(!is.na(q_value) & q_value < 0.01, "*", "")) %>%
      dplyr::select(study, central_dinucleotide, sig) %>%
      tidyr::pivot_wider(names_from = central_dinucleotide, values_from = sig) %>%
      as.data.frame()
    rownames(sig) <- sig$study
    sig <- as.matrix(sig[, -1])
    sig <- sig[rows, cols, drop = FALSE]

    spec <- make_diverging(mat, n = 200)

    # Dynamic size: less height if there are fewer studies
    n_rows <- length(rows)
    n_cols <- length(cols)
    pdf_height <- max(4, min(10, 0.28 * n_rows + 3))  # between 4 and 10
    pdf_width  <- max(5, min(8,  0.35 * n_cols + 2))  # between 5 and 8

    title_fontsize <- 6

    pdf(outfile, height = pdf_height, width = pdf_width, useDingbats = FALSE)
    ph <- pheatmap::pheatmap(
      mat,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = spec$colors,
      breaks = spec$breaks,
      main = title,
      display_numbers = sig,
      number_color = "black",
      fontsize_number = 12,
      fontsize_row = 7,
      fontsize_col = 7,
      scale = "none",
      silent = TRUE
    )

    # Search for and reduce the title font size robustly
    title_idx <- which(grepl("title", ph$gtable$layout$name, ignore.case = TRUE) |
               grepl("main", ph$gtable$layout$name, ignore.case = TRUE))

    if (length(title_idx) > 0) {
      for (ti in title_idx) {
      grob <- ph$gtable$grobs[[ti]]
      # If it is a textGrob with gp
      if (!is.null(grob$gp) && !is.null(grob$gp$fontsize)) {
        grob$gp <- grid::gpar(fontsize = title_fontsize)
      } else {
        # Try adjusting in children if it exists
        if (!is.null(grob$children)) {
        for (k in seq_along(grob$children)) {
          child <- grob$children[[k]]
          if (!is.null(child$gp) && !is.null(child$gp$fontsize)) {
          child$gp <- grid::gpar(fontsize = title_fontsize)
          grob$children[[k]] <- child
          } else if (inherits(child, "text")) {
          child$gp <- grid::gpar(fontsize = title_fontsize)
          grob$children[[k]] <- child
          }
        }
        } else {
        # Fallback: try to assign gp directly
        grob$gp <- grid::gpar(fontsize = title_fontsize)
        }
      }
      ph$gtable$grobs[[ti]] <- grob
      }
    }

    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
    dev.off()
    }

  dir.create("~/Dropbox/induction/5_remodeling/multivariate/promoters", recursive = TRUE, showWarnings = FALSE)

  plot_model_heatmap_small(
    lm_prom_fits,
    "Accessibility FC ~ Motif Enrichment (LM, promoters, min 1000 peaks)",
    "~/Dropbox/induction/5_remodeling/multivariate/promoters/accessibility_fc_motif_enrichment_lm_promoters_q0.01.pdf"
  )

  if (!is.null(lm_prom_acc_fits)) {
    plot_model_heatmap_small(
      lm_prom_acc_fits,
      "Accessibility FC ~ Motif Enrichment + Initial Accessibility (LM, promoters, min 1000 peaks)",
      "~/Dropbox/induction/5_remodeling/multivariate/promoters/accessibility_fc_motif_enrichment_plus_initial_accessibility_lm_promoters_q0.01.pdf"
    )
  }



##################################################################
##### fold change accessibility ~ motifs in enhancers  ###########
##################################################################

peak_and_flanking_counts_with_accessibility_enhancers<- peak_and_flanking_counts_with_accessibility %>%
  dplyr::filter( abs (distanceToTSS) >= 4000 )




  # Enhancers: mirror of the promoters analysis but for distal elements
  library(dplyr)
  library(tidyr)
  library(broom)

  # Count peaks per study (enhancers)
  peaks_per_study_enh <- peak_and_flanking_counts_with_accessibility_enhancers %>%
    select(study, peak_id) %>% distinct() %>%
    group_by(study) %>%
    summarise(n_peaks = n(), .groups = "drop")

  peaks_per_study_enh %>% arrange(n_peaks) %>% print(n = Inf, width = Inf)

  # Keep studies with >= 1000 peaks
  valid_studies_enh <- peaks_per_study_enh %>% filter(n_peaks >= 1000) %>% pull(study)

  peak_and_flanking_counts_with_accessibility_enhancers_filtered <- peak_and_flanking_counts_with_accessibility_enhancers %>%
    dplyr::filter(study %in% valid_studies_enh)

  # Build wide data frame for enhancers
  df_enh <- peak_and_flanking_counts_with_accessibility_enhancers_filtered %>%
    select(study, peak_id, central_dinucleotide, motif_enrichment, log2FoldChange,
           dplyr::any_of("log_accessibility")) %>%
    distinct() %>%
    mutate(central_dinucleotide = factor(central_dinucleotide, levels = motif_levels)) %>%
    pivot_wider(
      names_from = central_dinucleotide,
      values_from = motif_enrichment,
      names_prefix = "motif_enrichment_"
    )

  # LM fits for enhancers (no accessibility covariate)
  lm_enh_fits <- lapply(split(df_enh, df_enh$study), function(df_study) {
    curr_study <- unique(df_study$study)

    preds <- grep("^motif_enrichment_", colnames(df_study), value = TRUE)

    # Keep only finite response values
    df_study_clean <- df_study %>% filter(is.finite(log2FoldChange))

    preds_valid <- preds[vapply(df_study_clean[preds], function(x) {
      length(unique(x[!is.na(x)])) > 1
    }, logical(1))]

    if (length(preds_valid) == 0) return(NULL)
    if (nrow(df_study_clean) <= (length(preds_valid) + 1)) return(NULL)

    df_study_clean[preds_valid] <- lapply(df_study_clean[preds_valid], function(x) replace(x, is.na(x), 0))

    preds_escaped <- paste0("`", preds_valid, "`")
    fm <- as.formula(paste("log2FoldChange ~", paste(preds_escaped, collapse = " + ")))

    fit <- tryCatch(lm(fm, data = df_study_clean), error = function(e) NULL)
    if (is.null(fit)) return(NULL)

    broom::tidy(fit) %>%
      filter(term != "(Intercept)") %>%
      mutate(
        term = gsub("`", "", term),
        study = curr_study,
        central_dinucleotide = sub("^motif_enrichment_", "", term)
      )
  }) %>% bind_rows()

  print(lm_enh_fits)
  lm_enh_fits <- lm_enh_fits %>% mutate(q_value = p.adjust(p.value, method = "BH"))
  saveRDS(lm_enh_fits, "~/accessibility_fc~motif_enrichment_lm_enhancers_fits.RDS")

  # LM fits for enhancers controlling for initial accessibility
  fits_list_enh <- lapply(split(df_enh, df_enh$study), function(df_study) {
    if (!"log_accessibility" %in% colnames(df_study)) return(NULL)
    curr_study <- unique(df_study$study)

    preds <- grep("^motif_enrichment_", colnames(df_study), value = TRUE)

    df_clean <- df_study %>%
      filter(is.finite(log2FoldChange), is.finite(log_accessibility))

    if (nrow(df_clean) == 0 || length(unique(df_clean$log_accessibility)) <= 1) return(NULL)

    preds_valid <- preds[vapply(df_clean[preds], function(x) {
      length(unique(x[!is.na(x)])) > 1
    }, logical(1))]

    if (length(preds_valid) == 0) return(NULL)
    if (nrow(df_clean) <= (length(preds_valid) + 1)) return(NULL)

    df_clean[preds_valid] <- lapply(df_clean[preds_valid], function(x) replace(x, is.na(x), 0))

    fm <- as.formula(
      paste(
        "log2FoldChange ~",
        paste(c(paste0("`", preds_valid, "`"), "`log_accessibility`"), collapse = " + ")
      )
    )

    fit <- tryCatch(lm(fm, data = df_clean), error = function(e) NULL)
    if (is.null(fit)) return(NULL)

    broom::tidy(fit) %>%
      filter(term != "(Intercept)") %>%
      mutate(
        term = gsub("`", "", term),
        study = curr_study,
        central_dinucleotide = sub("^motif_enrichment_", "", term)
      )
  })

  fits_list_enh <- Filter(Negate(is.null), fits_list_enh)
  lm_enh_acc_fits <- if (length(fits_list_enh)) bind_rows(fits_list_enh) else NULL

  if (!is.null(lm_enh_acc_fits)) {
    lm_enh_acc_fits <- lm_enh_acc_fits %>% mutate(q_value = p.adjust(p.value, method = "BH"))
    saveRDS(lm_enh_acc_fits, "~/accessibility_fc~motif_enrichment_plus_initial_accessibility_lm_enhancers_fits.RDS")
  }

  # Ensure slope column exists for heatmaps
  if ("estimate" %in% colnames(lm_enh_fits) && !("slope" %in% colnames(lm_enh_fits))) {
    lm_enh_fits <- mutate(lm_enh_fits, slope = estimate)
  }
  if (!is.null(lm_enh_acc_fits) && ("estimate" %in% colnames(lm_enh_acc_fits)) && !("slope" %in% colnames(lm_enh_acc_fits))) {
    lm_enh_acc_fits <- mutate(lm_enh_acc_fits, slope = estimate)
  }

  # Plot heatmaps for enhancers (compact version used for promoters)
  dir.create("~/Dropbox/induction/5_remodeling/multivariate/enhancers", recursive = TRUE, showWarnings = FALSE)

  plot_model_heatmap_small(
    lm_enh_fits,
    "Accessibility FC ~ Motif Enrichment (LM, enhancers, min 1000 peaks)",
    "~/Dropbox/induction/5_remodeling/multivariate/enhancers/accessibility_fc_motif_enrichment_lm_enhancers_q0.01.pdf"
  )

  if (!is.null(lm_enh_acc_fits)) {
    plot_model_heatmap_small(
      lm_enh_acc_fits,
      "Accessibility FC ~ Motif Enrichment + Initial Accessibility (LM, enhancers, min 1000 peaks)",
      "~/Dropbox/induction/5_remodeling/multivariate/enhancers/accessibility_fc_motif_enrichment_plus_initial_accessibility_lm_enhancers_q0.01.pdf"
    )
  }



##### combine the following heatmaps into a single stacked heatmap with shared color scale:
### accessibility_fc_motif_enrichment_lm_promoters_q0.01.pdf, accessibility_fc_motif_enrichment_plus_initial_accessibility_lm_promoters_q0.01.pdf,
### accessibility_fc_motif_enrichment_lm_enhancers_q0.01.pdf, accessibility_fc_motif_enrichment_plus_initial_accessibility_lm_enhancers_q0.01.pdf
# short names for stacking and ordering
lm_prom_fits <- lm_prom_fits %>%
  mutate(region = "Promoters",
         study_region = paste(study, region, sep = " | "))
lm_prom_acc_fits <- lm_prom_acc_fits %>%
  mutate(region = "Promoters+Acc",
         study_region = paste(study, region, sep = " | "))
lm_enh_fits <- lm_enh_fits %>%
  mutate(region = "Enhancers",
         study_region = paste(study, region, sep = " | "))
lm_enh_acc_fits <- lm_enh_acc_fits %>%
  mutate(region = "Enhancers+Acc",
         study_region = paste(study, region, sep = " | "))
lm_combined_regions <- bind_rows(lm_prom_fits, lm_prom_acc_fits, lm_enh_fits, lm_enh_acc_fits)
# Build wide matrix stacking study x region as rows
df_slope_combined_regions <- lm_combined_regions %>%
  select(study, region, study_region, central_dinucleotide, slope) %>%
  group_by(study, region, study_region, central_dinucleotide) %>%
  summarise(slope = mean(as.numeric(unlist(slope)), na.rm = TRUE), .groups = "drop") %>%
  mutate(slope = ifelse(is.nan(slope), NA_real_, slope))
mat_df <- df_slope_combined_regions %>%
  pivot_wider(names_from = central_dinucleotide, values_from = slope) %>%
  as.data.frame()
rownames(mat_df) <- mat_df$study_region
mat <- as.matrix(mat_df[ , setdiff(colnames(mat_df), "study_region"), drop = FALSE])
# Determine row order: keep same study order (orden) and put Promoters, Promoters+Acc, Enhancers, Enhancers+Acc for each study
studies_present <- unique(lm_combined_regions$study)
base_order <- orden[orden %in% studies_present]
desired_rows <- c(
  paste(base_order, "Promoters", sep = " | "),
  paste(base_order, "Promoters+Acc", sep = " | "),
  paste(base_order, "Enhancers", sep = " | "),
  paste(base_order, "Enhancers+Acc", sep = " | ")
)
rows <- intersect(desired_rows, rownames(mat))
cols <- intersect(motif_levels, colnames(mat))


if (length(rows) == 0 || length(cols) == 0) {
  warning("No rows or columns to plot for combined regions heatmap")
} else {
  mat <- mat[rows, cols, drop = FALSE]

  # Build significance matrix (single symbol per study_region x motif)
  df_sig <- lm_combined_regions %>%
    dplyr::select(study_region, central_dinucleotide, q_value) %>%
    dplyr::group_by(study_region, central_dinucleotide) %>%
    dplyr::summarise(sig = if (any(!is.na(q_value) & q_value < 0.01)) "*" else "", .groups = "drop")

  sig_df <- df_sig %>%
    tidyr::pivot_wider(names_from = central_dinucleotide, values_from = sig) %>%
    as.data.frame()

  if (nrow(sig_df) == 0) {
    sig <- matrix("", nrow = nrow(mat), ncol = ncol(mat),
                  dimnames = list(rownames(mat), colnames(mat)))
  } else {
    rownames(sig_df) <- sig_df$study_region
    sig_mat <- as.matrix(sig_df[ , setdiff(colnames(sig_df), "study_region"), drop = FALSE])

    sig <- matrix("", nrow = nrow(mat), ncol = ncol(mat), dimnames = list(rownames(mat), colnames(mat)))
    common_rows <- intersect(rownames(sig_mat), rownames(mat))
    common_cols <- intersect(colnames(sig_mat), colnames(mat))
    if (length(common_rows) > 0 && length(common_cols) > 0) {
      sig[common_rows, common_cols] <- sig_mat[common_rows, common_cols, drop = FALSE]
    }
  }

  # Ensure numeric matrix and handle non-finite values
  mat_numeric <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat_numeric[!is.finite(mat_numeric)] <- NA_real_
  mat <- mat_numeric

  if (all(is.na(mat))) {
    warning("Combined regions matrix contains only NA values; skipping heatmap plot.")
  } else {
    spec <- make_diverging(mat, n = 200)

    out_file <- "~/Dropbox/induction/5_remodeling/multivariate/combined_regions_accessibility_lm_sharedscale.pdf"
    pdf(out_file, height = 15, width = 9, useDingbats = FALSE)
    pheatmap::pheatmap(
      mat,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = spec$colors,
      breaks = spec$breaks,
      main = "Accessibility FC ~ Motif Enrichment (Promoters/Enhancers, LM) [shared scale]",
      display_numbers = sig,
      number_color = "black",
      fontsize_number = 10,
      fontsize_row = 8,
      fontsize_col = 8,
      scale = "none"
    )
    dev.off()
  }
}







##################################################
############### expression #######################
##################################################
########### Main figure 5b and Supplementary figure 7b 


### ebox counts with flanking regions
peak_and_flanking_counts<- readRDS("~/peak_and_flanking_counts.RDS")

## expression info per peak
canonical_eboxes_fc_info<- readRDS("~/proneural/canonical_eboxes_with_expression_accessibility_info.RDS")

expression_info_by_peak<-canonical_eboxes_fc_info %>%
  dplyr::select(peak_id, study, log2FoldChange) %>% distinct()



### show number of distinct peaks per study
peaks_per_study <- expression_info_by_peak %>%
  group_by(study) %>%
  summarise(n_peaks = n(), .groups = "drop")

peaks_per_study %>% arrange(n_peaks) %>% print(n=Inf, width=Inf)

## get studies that have at least 500 peaks
studies_with_1000_peaks <- peaks_per_study %>%
  filter(n_peaks >= 1000) %>%
  pull(study)







expression_studies<-unique(canonical_eboxes_fc_info$study)


peak_and_flanking_counts_filtered <- peak_and_flanking_counts %>%
  dplyr::filter(study %in% expression_studies)


### add expression info to counts, some peaks of peak_and_flanking_counts_filtered
### will not have expression info, so do not use those rows
peak_and_flanking_counts_with_expression <- peak_and_flanking_counts_filtered %>%
  dplyr::left_join(expression_info_by_peak, by = c("peak_id", "study")) %>%
  dplyr::filter(!is.na(log2FoldChange)) %>% distinct()


# remove specific low-peak studies from all datasets (less than 100 peaks)
studies_to_drop <- c(
  "MSGN1 -- EB [Chalamalasetty et al. 2014]",
  "NEUROD1 -- mESC [Pataskar et al. 2016]"
)



##### do rigde regression of expression ~ motif enrichments CAT-CAT + CAT-CAG + CAG-CAG + CAG-CAC + CAC-CAC + CAT-CAC + CAA-CAT + CAA-CAG + CAA-CAC + CAA-CAA
library(glmnet)
library(dplyr)
library(tidyr)

motif_levels <- c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAC-CAC",
                  "CAT-CAC","CAA-CAT","CAA-CAG","CAA-CAC","CAA-CAA")




df_expr <- peak_and_flanking_counts_with_expression %>%
  filter(!study %in% studies_to_drop) %>%
  select(study, peak_id, central_dinucleotide, motif_enrichment, log2FoldChange) %>%
  distinct() %>%
  mutate(central_dinucleotide = factor(central_dinucleotide, levels = motif_levels)) %>%
  pivot_wider(
    names_from = central_dinucleotide,
    values_from = motif_enrichment,
    names_prefix = "motif_enrichment_"
  )


df_expr_acc <- peak_and_flanking_counts_with_expression %>%
  filter(!study %in% studies_to_drop) %>%
  select(study, peak_id, central_dinucleotide, motif_enrichment, log2FoldChange, log_accessibility) %>%
  distinct() %>%
  mutate(central_dinucleotide = factor(central_dinucleotide, levels = motif_levels)) %>%
  pivot_wider(
    names_from = central_dinucleotide,
    values_from = motif_enrichment,
    names_prefix = "motif_enrichment_"
  )




library(dplyr)
library(tidyr)
library(broom)

lm_fits <- lapply(split(df_expr_acc, df_expr_acc$study), function(df_study) {
    curr_study <- unique(df_study$study)
    
    # 1. Identify predictors
    preds <- grep("^motif_enrichment_", colnames(df_study), value = TRUE)
    
    # 2. Clean rows with log2FoldChange NA
    df_study_clean <- df_study %>% tidyr::drop_na(log2FoldChange)
    
    # 3. Variance Check
    preds_validos <- preds[vapply(df_study_clean[preds], function(x) {
        length(unique(x[!is.na(x)])) > 1
    }, logical(1))]
    
    if (length(preds_validos) == 0) return(NULL)
    if (nrow(df_study_clean) <= (length(preds_validos) + 1)) return(NULL)

    # 4. Impute NAs
    df_study_clean[preds_validos] <- lapply(df_study_clean[preds_validos], 
                                            function(x) replace(x, is.na(x), 0))

    # 5. BUILD FORMULA WITH BACKTICKS (The key fix)
    # This transforms 'motif_enrichment_CAA-CAA' into '`motif_enrichment_CAA-CAA`'
    preds_escaped <- paste0("`", preds_validos, "`")
    fm <- as.formula(paste("log2FoldChange ~", paste(preds_escaped, collapse = " + ")))
    
    # 6. Run model
    fit <- tryCatch(lm(fm, data = df_study_clean), error = function(e) NULL)
    
    if (is.null(fit)) return(NULL)

    # 7. Extract results
    res <- broom::tidy(fit) %>%
        dplyr::filter(term != "(Intercept)") %>%
        # Clean the backticks from the term name for the final result
        dplyr::mutate(
            term = gsub("`", "", term),
            study = curr_study,
            central_dinucleotide = sub("^motif_enrichment_", "", term)
        )
    
    return(res)
}) %>% dplyr::bind_rows()

# Now it should show data
print(lm_fits)
# See final result
print(lm_fits)

lm_fits <- lm_fits %>%
    dplyr::mutate(q_value = p.adjust(p.value, method = "BH"))

saveRDS(lm_fits, "~/expression~motif_enrichment_lm_fits.RDS")


######## The same but controlling for initial accessibility
library(dplyr)
library(tidyr)
library(broom)

lm_acc_fits <- lapply(split(df_expr_acc, df_expr_acc$study), function(df_study) {
    curr_study <- unique(df_study$study)

    preds <- grep("^motif_enrichment_", colnames(df_study), value = TRUE)

    df_study_clean <- df_study %>%
        tidyr::drop_na(log2FoldChange, log_accessibility)

    preds_valid <- preds[vapply(df_study_clean[preds], function(x) {
        length(unique(x[!is.na(x)])) > 1
    }, logical(1))]

    if (length(preds_valid) == 0) return(NULL)
    if (nrow(df_study_clean) <= (length(preds_valid) + 1)) return(NULL)

    df_study_clean[preds_valid] <- lapply(df_study_clean[preds_valid], function(x) replace(x, is.na(x), 0))

    preds_escaped <- paste0("`", preds_valid, "`")
    fm <- as.formula(paste("log2FoldChange ~", paste(c(preds_escaped, "`log_accessibility`"), collapse = " + ")))

    fit <- tryCatch(lm(fm, data = df_study_clean), error = function(e) NULL)
    if (is.null(fit)) return(NULL)

    broom::tidy(fit) %>%
        dplyr::filter(term != "(Intercept)") %>%
        dplyr::mutate(
            term = gsub("`", "", term),
            study = curr_study,
            central_dinucleotide = sub("^motif_enrichment_", "", term)
        )
}) %>% dplyr::bind_rows()

print(lm_acc_fits)

lm_acc_fits <- lm_acc_fits %>%
    dplyr::mutate(q_value = p.adjust(p.value, method = "BH"))

saveRDS(lm_acc_fits, "~/expression~motif_enrichment_plus_accessibility_lm_fits.RDS")


#### read all the models for the heatmap
readRDS("~/expression~motif_enrichment_lm_fits.RDS")-> lm_fits
readRDS("~/expression~motif_enrichment_plus_accessibility_lm_fits.RDS")-> lm_acc_fits






############################### HEATMAPS
library(ggplot2)
library(dplyr)
# Prepare and plot heatmaps for the four model result objectss

motif_levels <- c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAC-CAC",
          "CAT-CAC","CAA-CAT","CAA-CAG","CAA-CAC","CAA-CAA")

plot_model_heatmap <- function(df, title, outfile) {
  df <- df %>%
  dplyr::mutate(
    study = factor(study, levels = orden),
    central_dinucleotide = factor(central_dinucleotide, levels = motif_levels)
  )

  # Summarise duplicates: ensure a single numeric slope per study+motif
  df_slope <- df %>%
    dplyr::select(study, central_dinucleotide, slope) %>%
    dplyr::group_by(study, central_dinucleotide) %>%
    dplyr::summarise(
      slope = mean(as.numeric(unlist(slope)), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(slope = ifelse(is.nan(slope), NA_real_, slope))

  mat <- df_slope %>%
  tidyr::pivot_wider(names_from = central_dinucleotide, values_from = slope) %>%
  as.data.frame()

  if (nrow(mat) == 0) {
    warning("No data to plot for: ", title)
    return(invisible(NULL))
  }

  rownames(mat) <- mat$study
  mat <- as.matrix(mat[,-1])
  # Safe row/column selection to avoid 'subscript out of bounds'
  rows <- intersect(orden, rownames(mat))
  cols <- intersect(motif_levels, colnames(mat))
  if (length(rows) == 0 || length(cols) == 0) {
    warning("No rows or columns to plot for: ", title)
    return(invisible(NULL))
  }
  mat <- mat[rows, cols, drop = FALSE]

  # Summarise significance (q_value) into a single symbol per study+motif
  df_sig <- df %>%
    dplyr::select(study, central_dinucleotide, q_value) %>%
    dplyr::group_by(study, central_dinucleotide) %>%
    dplyr::summarise(sig = if (any(!is.na(q_value) & q_value < 0.01)) "*" else "", .groups = "drop")

  sig <- df_sig %>%
  tidyr::pivot_wider(names_from = central_dinucleotide, values_from = sig) %>%
  as.data.frame()
  if (nrow(sig) == 0) {
    # create empty sig matrix matching mat dimensions
    sig <- matrix("", nrow = nrow(mat), ncol = ncol(mat),
                  dimnames = list(rownames(mat), colnames(mat)))
  } else {
    rownames(sig) <- sig$study
    sig <- as.matrix(sig[,-1, drop = FALSE])
    # Ensure same rows/cols as mat (fill missing with "")
    all_rows <- rows
    all_cols <- cols
    # Reorder/expand sig to match mat rows/cols
    sig_full <- matrix("", nrow = length(all_rows), ncol = length(all_cols),
                       dimnames = list(all_rows, all_cols))
    common_rows <- intersect(rownames(sig), all_rows)
    common_cols <- intersect(colnames(sig), all_cols)
    if (length(common_rows) > 0 && length(common_cols) > 0) {
      sig_full[common_rows, common_cols] <- sig[common_rows, common_cols, drop = FALSE]
    }
    sig <- sig_full
  }

  spec <- make_diverging(mat, n = 200)

  pdf(outfile, height = 15, width = 9, useDingbats = FALSE)
  pheatmap::pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = spec$colors,
  breaks = spec$breaks,
  main = title,
  display_numbers = sig,
  number_color = "black",
  fontsize_number = 14,
  fontsize_row = 8,
  fontsize_col = 8,
  scale = "none"
  )
  dev.off()
}




#### the heatmaps filtereign to get studies that have at least 1000 peaks with expression info (studies_with_1000_peaks)

# Heatmaps restricted to studies with ≥1000 peaks having expression info
library(dplyr)

# Filter model results
er_fits_1000      <- expression_ridge_fits %>% filter(study %in% studies_with_1000_peaks)
era_fits_1000     <- expression_ridge_acc_fits %>% filter(study %in% studies_with_1000_peaks)
lm_fits_1000      <- lm_fits %>% filter(study %in% studies_with_1000_peaks)
lm_acc_fits_1000  <- lm_acc_fits %>% filter(study %in% studies_with_1000_peaks)

# Helper to set/restore ordering and plot
plot_min1000 <- function(df, title, outfile) {
  if (!nrow(df)) return(invisible(NULL))
  old_orden <- orden
  orden <- old_orden[old_orden %in% unique(df$study)]
  on.exit(orden <- old_orden, add = TRUE)
  plot_model_heatmap(df, title, outfile)
}


plot_min1000(lm_fits_1000,
             "Expression ~ Motif Enrichment (LM) [min 1000 peaks]",
             "~/Dropbox/induction/7_expression/multivariate/min_1000peaks/expression_motif_enrichment_lm_min1000peaks_q0.01.pdf")

plot_min1000(lm_acc_fits_1000,
             "Expression ~ Motif Enrichment + Accessibility (LM) [min 1000 peaks]",
             "~/Dropbox/induction/7_expression/multivariate/min_1000peaks/expression_motif_enrichment_plus_accessibility_lm_min1000peaks_q0.01.pdf")

#### concatenate lm_fits_1000 and lm_acc_fits_1000 into a single stacked heatmap sharing the same color scale

# short names for stacking and ordering
lm_fits_1000 <- lm_fits_1000 %>%
  mutate(model = "LM",
          study_model = paste(study, model, sep = " | "))
lm_acc_fits_1000 <- lm_acc_fits_1000 %>%
  mutate(model = "LM+Acc",
          study_model = paste(study, model, sep = " | "))

lm_combined <- bind_rows(lm_fits_1000, lm_acc_fits_1000)

# ensure slope column exists
if ("estimate" %in% colnames(lm_combined) && !("slope" %in% colnames(lm_combined))) {
  lm_combined <- dplyr::mutate(lm_combined, slope = estimate)
}

# Build wide matrix stacking study x model as rows
df_slope_combined <- lm_combined %>%
  dplyr::select(study, model, study_model, central_dinucleotide, slope) %>%
  dplyr::group_by(study, model, study_model, central_dinucleotide) %>%
  dplyr::summarise(slope = mean(as.numeric(unlist(slope)), na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(slope = ifelse(is.nan(slope), NA_real_, slope))

mat_df <- df_slope_combined %>%
  tidyr::pivot_wider(names_from = central_dinucleotide, values_from = slope) %>%
  as.data.frame()

if (nrow(mat_df) == 0) {
  warning("No data to plot for combined LM heatmap")
} else {
  rownames(mat_df) <- mat_df$study_model
  mat <- as.matrix(mat_df[ , setdiff(colnames(mat_df), "study_model"), drop = FALSE])

  # Determine row order: keep same study order (orden) and put LM then LM+Acc for each study
  studies_present <- unique(lm_combined$study)
  base_order <- orden[orden %in% studies_present]
  desired_rows <- c(paste(base_order, "LM", sep = " | "), paste(base_order, "LM+Acc", sep = " | "))
  rows <- intersect(desired_rows, rownames(mat))
  cols <- intersect(motif_levels, colnames(mat))
  if (length(rows) == 0 || length(cols) == 0) {
    warning("No rows or columns to plot for combined LM heatmap")
  } else {
    mat <- mat[rows, cols, drop = FALSE]

    # build significance matrix (single symbol per study_model x motif)
    df_sig <- lm_combined %>%
      dplyr::select(study_model, central_dinucleotide, q_value) %>%
      dplyr::group_by(study_model, central_dinucleotide) %>%
      dplyr::summarise(sig = if (any(!is.na(q_value) & q_value < 0.01)) "*" else "", .groups = "drop")

    sig_df <- df_sig %>%
      tidyr::pivot_wider(names_from = central_dinucleotide, values_from = sig) %>%
      as.data.frame()
    if (nrow(sig_df) == 0) {
      sig <- matrix("", nrow = nrow(mat), ncol = ncol(mat), dimnames = list(rownames(mat), colnames(mat)))
    } else {
      rownames(sig_df) <- sig_df$study_model
      sig_mat <- as.matrix(sig_df[ , setdiff(colnames(sig_df), "study_model"), drop = FALSE])
      # build full sig matrix matching mat rows/cols, fill missing with ""
      sig <- matrix("", nrow = nrow(mat), ncol = ncol(mat), dimnames = list(rownames(mat), colnames(mat)))
      common_rows <- intersect(rownames(sig_mat), rownames(mat))
      common_cols <- intersect(colnames(sig_mat), colnames(mat))
      if (length(common_rows) > 0 && length(common_cols) > 0) {
        sig[common_rows, common_cols] <- sig_mat[common_rows, common_cols, drop = FALSE]
      }
    }

    # compute shared color scale across the combined matrix
    # Ensure 'mat' is a numeric matrix (it may have been coerced to character)
    mat <- as.matrix(mat)
    mat_numeric <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat),
                          dimnames = dimnames(mat))
    # Replace non-finite values with NA
    mat_numeric[!is.finite(mat_numeric)] <- NA_real_
    mat <- mat_numeric

    # If the matrix contains no numeric data, skip plotting to avoid cut.default errors
    if (all(is.na(mat))) {
      warning("Combined LM matrix contains only NA values; skipping heatmap plot.")
    } else {
      spec <- make_diverging(mat, n = 200)

      out_file <- "~/Dropbox/induction/7_expression/multivariate/min_1000peaks/expression_motif_enrichment_lm_combined_min1000peaks_q0.01_sharedscale.pdf"
      pdf(out_file, height = 15, width = 9, useDingbats = FALSE)
      pheatmap::pheatmap(
        mat,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        color = spec$colors,
        breaks = spec$breaks,
        main = "Expression ~ Motif Enrichment (LM vs LM+Acc) [min 1000 peaks]",
        display_numbers = sig,
        number_color = "black",
        fontsize_number = 10,
        fontsize_row = 8,
        fontsize_col = 8,
        scale = "none"
      )
      dev.off()
    }
  }
}





########################################################################
############# grouped gene ontologies ##################################
###########################################################################
######################## Main figure 5c




############### I can group the Gene Ontologies, and then I can plot each group together for each gene.
### Plotting the variation in terms of the NES, and weighting it by the p-value—for example, using boxplots.
## And discarding the ontologies that do not fit into any group.

## Actually, better than a boxplot, I'll plot balloons for each ontology on a Y-axis. All these balloons packed within the same gene/category.

### I have only included the groups that have enough ontologies.


# Create a list of grouped ontologies and their associated ontologies
grouped_ontology_list <- list(
    "Cytoskeletal organization and ciliary motility" = c(
        "intermediate filament organization",
        "intermediate filament cytoskeleton organization",
        "intermediate filament-based process",
        "microtubule-based movement",
        "cilium assembly",
        "cilium organization",
        "cilium or flagellum-dependent cell motility",
        "cilium-dependent cell motility",
        "sperm motility",
        "flagellated sperm motility",
        "cilium movement",
        "cilium movement involved in cell motility",
        "axoneme assembly",
        "microtubule bundle formation"
    ),
    "Tissue repair & chemotaxis" = c(
        "leukocyte migration",
        "leukocyte chemotaxis",
        "granulocyte chemotaxis",
        "myeloid leukocyte migration",
        "regulation of chemotaxis",
        "chemotaxis",
        "taxis",
        "wound healing",
        "positive regulation of endothelial cell migration",
        "cellular response to alcohol"
    ),
    "Cardiac and vascular development" = c(
        "cardiac muscle tissue development",
        "endocardial cushion development",
        "cardiac chamber development",
        "cardiac chamber morphogenesis",
        "smooth muscle cell differentiation",
        "mesenchyme morphogenesis",
        "skeletal system morphogenesis",
        "pattern specification process",
        "regionalization",
        "heart valve development",
        "heart valve morphogenesis",
        "artery development",
        "artery morphogenesis",
        "cardiac septum development",
        "semi-lunar valve development",
        "anatomical structure arrangement",
        "cranial nerve morphogenesis",
        "motor neuron axon guidance"
    ),
    "Innate and adaptive immune defense" = c(
        "defense response to bacterium",
        "defense response to virus",
        "response to virus",
        "cellular response to interferon-beta",
        "response to interferon-beta",
        "cell killing",
        "leukocyte mediated immunity",
        "immune effector process",
        "positive regulation of immune effector process",
        "regulation of innate immune response",
        "negative regulation of lymphocyte mediated immunity",
        "negative regulation of leukocyte mediated immunity",
        "negative regulation of immune effector process",
        "lymphocyte mediated immunity",
        "mesodermal cell differentiation",
        "positive regulation of apoptotic signaling pathway",
        "prostaglandin metabolic process",
        "prostanoid metabolic process",
        "antimicrobial humoral immune response mediated by antimicrobial peptide",
        "antimicrobial humoral response",
        "humoral immune response",
        "disruption of anatomical structure in another organism",
        "killing of cells of another organism",
        "disruption of cell in another organism"
    ),
    "Neural development, synaptic organization, and excitability" = c(
        "synaptic transmission, glutamatergic",
        "protein localization to synapse",
        "action potential",
        "potassium ion transport",
        "potassium ion transmembrane transport",
        "regulation of postsynaptic membrane potential",
        "regulation of membrane potential",
        "chemical synaptic transmission, postsynaptic",
        "neuromuscular process",
        "retina development in camera-type eye",
        "segmentation",
        "spinal cord development",
        "neuron fate commitment",
        "learning",
        "associative learning",
        "limbic system development",
        "pallium development",
        "Notch signaling pathway",
        "regulation of trans-synaptic signaling",
        "axon guidance",
        "neuron projection guidance",
        "multicellular organismal response to stress",
        "cell fate commitment",
        "adult behavior",
        "positive regulation of neurotransmitter secretion",
        "positive regulation of neurotransmitter transport",
        "cell-cell adhesion via plasma-membrane adhesion molecules",
        "homophilic cell adhesion via plasma membrane adhesion molecules",
        "axonogenesis",
        "developmental cell growth",
        "regulation of axonogenesis",
        "neuron projection extension",
        "central nervous system neuron differentiation",
        "central nervous system neuron development",
        "forebrain development",
        "telencephalon development",
        "positive regulation of synapse assembly",
        "regulation of synapse assembly",
        "regulation of synapse structure or activity",
        "regulation of synapse organization",
        "neuron migration",
        "synapse assembly"
    ),
    "Skeletal muscle development and contraction" = c(
        "cellular component assembly involved in morphogenesis",
        "cellular anatomical entity morphogenesis",
        "sarcomere organization",
        "myofibril assembly",
        "skeletal muscle contraction",
        "multicellular organismal movement",
        "musculoskeletal movement",
        "striated muscle tissue development",
        "muscle organ development",
        "skeletal muscle organ development",
        "skeletal muscle tissue development",
        "muscle system process",
        "muscle contraction",
        "regulation of muscle system process",
        "muscle tissue development",
        "muscle cell development",
        "striated muscle cell development",
        "striated muscle contraction",
        "muscle cell differentiation",
        "striated muscle cell differentiation"
    ),
    "DNA replication and mitotic chromosome segregation" = c(
        "nuclear DNA replication",
        "cell cycle DNA replication",
        "DNA replication",
        "DNA-templated DNA replication",
        "negative regulation of mitotic nuclear division",
        "negative regulation of nuclear division",
        "nuclear division",
        "mitotic nuclear division",
        "sister chromatid segregation",
        "regulation of chromosome segregation",
        "nuclear chromosome segregation",
        "chromosome segregation",
        "mitotic sister chromatid segregation",
        "metaphase chromosome alignment",
        "mitotic metaphase chromosome alignment",
        "spindle checkpoint signaling",
        "negative regulation of chromosome segregation",
        "negative regulation of metaphase/anaphase transition of cell cycle",
        "negative regulation of chromosome separation",
        "mitotic spindle assembly checkpoint signaling",
        "negative regulation of mitotic sister chromatid separation",
        "negative regulation of mitotic metaphase/anaphase transition",
        "negative regulation of sister chromatid segregation",
        "negative regulation of mitotic sister chromatid segregation",
        "regulation of mitotic sister chromatid segregation",
        "regulation of mitotic sister chromatid separation",
        "mitotic sister chromatid separation",
        "regulation of chromosome separation",
        "chromosome separation"
    ),
    "Ribosome biogenesis and protein translation" = c(
        "cytoplasmic translation",
        "translation at presynapse",
        "translation at synapse",
        "translation at postsynapse",
        "mitochondrial respiratory chain complex assembly",
        "ribosome assembly",
        "mitochondrial translation",
        "mitochondrial gene expression",
        "ribonucleoprotein complex biogenesis",
        "ribosome biogenesis",
        "rRNA processing",
        "rRNA metabolic process",
        "maturation of SSU-rRNA",
        "ribosomal small subunit biogenesis",
        "ribosomal large subunit biogenesis",
        "protein localization to nucleolus",
        "tRNA metabolic process",
        "tRNA processing",
        "tRNA modification"
    )
)

# Convert to dataframe
ontology_df <- do.call(rbind, lapply(names(grouped_ontology_list), function(group) {
    data.frame(
        grouped_ontology = group,
        ontology = grouped_ontology_list[[group]],
        stringsAsFactors = FALSE
    )
}))





#### add to balloon_data_all the grouped ontology
### remove the ontologies that are not in the ontology_df
balloon_data_all_grouped <- balloon_data_all %>%
    inner_join(ontology_df, by = c("Ontology" = "ontology"))










    # Calculate signed log10(padj): sign(NES) * -log10(padj)
    balloon_data_all_grouped$signed_log10padj <- sign(balloon_data_all_grouped$NES) * balloon_data_all_grouped$log10padj

    # Prepare matrix: rows = ontologies, columns = studies, values = signed_log10padj
    ontology_matrix_grouped_signedlogp <- balloon_data_all_grouped %>%
        dplyr::select(Ontology, Study, signed_log10padj, grouped_ontology) %>%
        tidyr::pivot_wider(names_from = Study, values_from = signed_log10padj)

    # Order ontologies as in ontology_balloon_all_gsea.pdf
    ordered_ontologies_grouped_signedlogp <- ordered_ontologies[ordered_ontologies %in% balloon_data_all_grouped$Ontology]

    # Remove columns not needed for heatmap
    ontology_matrix_grouped_signedlogp <- as.data.frame(ontology_matrix_grouped_signedlogp)
    rownames(ontology_matrix_grouped_signedlogp) <- ontology_matrix_grouped_signedlogp$Ontology
    ontology_matrix_grouped_signedlogp$Ontology <- NULL
    ontology_matrix_grouped_signedlogp$grouped_ontology <- NULL

    # Subset and order rows
    ontology_matrix_grouped_signedlogp <- ontology_matrix_grouped_signedlogp[ordered_ontologies_grouped_signedlogp, , drop = FALSE]

    # Convert to matrix
    ontology_matrix_grouped_signedlogp <- as.matrix(ontology_matrix_grouped_signedlogp)

    # Order columns as before
    ontology_matrix_grouped_signedlogp <- ontology_matrix_grouped_signedlogp[, orden]

    # Set color palette with white for NA, using min/max of signed_log10padj
    signedlogp_min <- min(ontology_matrix_grouped_signedlogp, na.rm = TRUE)
    signedlogp_max <- max(ontology_matrix_grouped_signedlogp, na.rm = TRUE)
    breaks_signedlogp <- c(
        seq(signedlogp_min, 0, length.out = 51),
        seq(0, signedlogp_max, length.out = 51)[-1]
    )
    my_palette<- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

    # Plot heatmap for signed log10(padj), NA as white
    pdf("~/Dropbox/induction/7_expression/ontology_grouped_signedlog10padj_heatmap.pdf", width = 15, height = 10)
    pheatmap(
        ontology_matrix_grouped_signedlogp,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        scale = "none",
        show_rownames = FALSE,
        annotation_row = row_anno,
        annotation_colors = anno_colors,
        main = "Signed -log10(adj p) for grouped ontologies (heatmap)",
        color = my_palette,
        breaks = breaks_signedlogp,
        fontsize_col = 10,
        angle_col = 45,
        treeheight_row = 0,
        gaps_row = gaps_row,
        na_col = "white"
    )
    dev.off()












