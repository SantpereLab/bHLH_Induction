
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






  # Modelos lineales con log2FoldChange de accesibilidad en promotores
  # (idénticos en construcción a los de expresión)
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
    
    # 1. Identificar predictores
    preds <- grep("^motif_enrichment_", colnames(df_study), value = TRUE)
    
    # 2. LIMPIEZA CRÍTICA: Mantener solo valores numéricos reales y finitos
    df_study_clean <- df_study %>% 
        dplyr::filter(is.finite(log2FoldChange))
    
    # 3. Validación de varianza en los predictores
    preds_validos <- preds[vapply(df_study_clean[preds], function(x) {
        length(unique(x[!is.na(x)])) > 1
    }, logical(1))]
    
    # Filtros de seguridad
    if (length(preds_validos) == 0) return(NULL)
    if (nrow(df_study_clean) <= (length(preds_validos) + 1)) return(NULL)

    # 4. Imputación de NAs en predictores
    df_study_clean[preds_validos] <- lapply(df_study_clean[preds_validos], 
                                            function(x) replace(x, is.na(x), 0))

    # 5. Construcción de fórmula
    preds_escaped <- paste0("`", preds_validos, "`")
    fm <- as.formula(paste("log2FoldChange ~", paste(preds_escaped, collapse = " + ")))
    
    # 6. Ejecución del modelo
    fit <- tryCatch(lm(fm, data = df_study_clean), error = function(e) NULL)
    
    if (is.null(fit)) return(NULL)

    # 7. Formateo de resultados
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




  # LM en promotores controlando por accesibilidad inicial
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

  # Asegurar columna 'slope' para heatmaps
  if ("estimate" %in% colnames(lm_prom_fits) && !("slope" %in% colnames(lm_prom_fits))) {
    lm_prom_fits <- dplyr::mutate(lm_prom_fits, slope = estimate)
  }
  if (!is.null(lm_prom_acc_fits) &&
      ("estimate" %in% colnames(lm_prom_acc_fits)) &&
      !("slope" %in% colnames(lm_prom_acc_fits))) {
    lm_prom_acc_fits <- dplyr::mutate(lm_prom_acc_fits, slope = estimate)
  }


  # Versión compacta del heatmap para promotores (PDF más pequeño y dinámico según nº de estudios)
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

    # Tamaño dinámico: menos alto si hay menos estudios
    n_rows <- length(rows)
    n_cols <- length(cols)
    pdf_height <- max(4, min(10, 0.28 * n_rows + 3))  # entre 4 y 10
    pdf_width  <- max(5, min(8,  0.35 * n_cols + 2))  # entre 5 y 8

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

    # Buscar y reducir el tamaño de la fuente del título de forma robusta
    title_idx <- which(grepl("title", ph$gtable$layout$name, ignore.case = TRUE) |
               grepl("main", ph$gtable$layout$name, ignore.case = TRUE))

    if (length(title_idx) > 0) {
      for (ti in title_idx) {
      grob <- ph$gtable$grobs[[ti]]
      # Si es un textGrob con gp
      if (!is.null(grob$gp) && !is.null(grob$gp$fontsize)) {
        grob$gp <- grid::gpar(fontsize = title_fontsize)
      } else {
        # Intentar ajustar en children si existe
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
        # Fallback: intentar asignar gp directamente
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
library(dplyr)
library(tidyr)
library(broom)

lm_fits <- lapply(split(df_expr_acc, df_expr_acc$study), function(df_study) {
    curr_study <- unique(df_study$study)
    
    # 1. Identificar predictores
    preds <- grep("^motif_enrichment_", colnames(df_study), value = TRUE)
    
    # 2. Limpiar filas con log2FoldChange NA
    df_study_clean <- df_study %>% tidyr::drop_na(log2FoldChange)
    
    # 3. Check de Varianza
    preds_validos <- preds[vapply(df_study_clean[preds], function(x) {
        length(unique(x[!is.na(x)])) > 1
    }, logical(1))]
    
    if (length(preds_validos) == 0) return(NULL)
    if (nrow(df_study_clean) <= (length(preds_validos) + 1)) return(NULL)

    # 4. Imputar NAs
    df_study_clean[preds_validos] <- lapply(df_study_clean[preds_validos], 
                                            function(x) replace(x, is.na(x), 0))

    # 5. CONSTRUIR FÓRMULA CON BACKTICKS (La corrección clave)
    # Esto transforma 'motif_enrichment_CAA-CAA' en '`motif_enrichment_CAA-CAA`'
    preds_escaped <- paste0("`", preds_validos, "`")
    fm <- as.formula(paste("log2FoldChange ~", paste(preds_escaped, collapse = " + ")))
    
    # 6. Ejecutar modelo
    fit <- tryCatch(lm(fm, data = df_study_clean), error = function(e) NULL)
    
    if (is.null(fit)) return(NULL)

    # 7. Extraer resultados
    res <- broom::tidy(fit) %>%
        dplyr::filter(term != "(Intercept)") %>%
        # Limpiar los backticks del nombre del término para el resultado final
        dplyr::mutate(
            term = gsub("`", "", term),
            study = curr_study,
            central_dinucleotide = sub("^motif_enrichment_", "", term)
        )
    
    return(res)
}) %>% dplyr::bind_rows()

# Ahora sí debería mostrar datos
print(lm_fits)
# Ver resultado final
print(lm_fits)

lm_fits <- lm_fits %>%
    dplyr::mutate(q_value = p.adjust(p.value, method = "BH"))

saveRDS(lm_fits, "~/expression~motif_enrichment_lm_fits.RDS")




######## lo mismo pero controlando por accesibilidad inicial
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




#################### correlation matrix ###################
###########################################################
###### Supplementary figure 7c



setwd("/projects_ng/shared/joel_induction/06_fc")

fc_files <- c(
    "Aydin_2019_Ascl1_eb_12h_DESeq2_results1.csv",
    "Aydin_2019_Neurog2_eb_48h_DESeq2_results1.csv",
    "Barfeld_2017_Myc_lncap_5h_microarray_results2_ENSEMBL.csv",
    "Casey_2018_Ascl1_mesc_DESeq2_results1.csv",
    "Casey_2018_Ascl2_mesc_ESeq2_results1.csv",
    "Casey_2018_Myod1_mesc_DESeq2_results1.csv",
    "Chalamasetty_2014_Msgn1_eb_48h_microarray_results2_ENSEMBL.csv",
    "Conerly_2016_Myf5_mef_DESeq2_results1.csv",
    "Conerly_2016_Myod1_mef_DESeq2_results1.csv",
    "Fong_2012_Myod1_p19_microarray_results2_ENSEMBL.csv",
    "Fong_2012_Myod1_mef_microarray_results2_ENSEMBL.csv",
    "Fong_2012_Neurod2_p19_microarray_results2_ENSEMBL.csv",
    "Fong_2012_Neurod2_mef_microarray_results2_ENSEMBL.csv",
    "Fong_2015_Neurod2_p19_DESeq2_results1.csv",
    "Jaenicke_2016_Myc_imec_DESeq2_results1.csv",
    "Joung_2023_Ascl1_hesc_day7_DESeq2_results1.csv",
    "Jung_2017_Myc_u2os_DESeq2_results1.csv",
    "Lee_2020_Myod1_mef_DESeq2_results1.csv",
    "Li_2019_Twist2_myoblast_DM_DESeq2_results1.csv",
    "Li_2019_Twist2_myoblast_GM_DESeq2_results1.csv",
    "Liu_2014_Ascl2_cd4_microarray_results2_ENSEMBL.csv",
    "Lorenzin_2016_Myc_u2os_DESeq2_results1.csv",
    "Manandhar_2016_Myod1_fibroblast_DESeq2_results1.csv",
    "Matsuda_2016_Neurod1_microglia_Day2_DESeq2_results1.csv",
    "Park_2017_Ascl1_glioblastoma_DESeq2_results1.csv",
    "Pataskar_2016_Neurod1_mesc_DESeq2_results1.csv",
    "Pereira_2024_Neurog2_astrocyte_DESeq2_results1.csv",
    "Smith_2016_Neurog2_mrc5_DESeq2_results1.csv",
    "Vainorius_2023_Ascl1_mesc_DESeq2_results1.csv",
    "Vainorius_2023_Neurog2_mesc_DESeq2_results1.csv",
    "Walz_2014_Myc_u2os_DESeq2_results1.csv",
    "Wang_2023_Ascl1_neuroblastoma_DESeq2_results1.csv",
    "Wapinski_2013_Ascl1_mef_DESeq2_results1.csv",
    "Weber_2016_Hey1_cardiomyocite_DESeq2_results1.csv",
    "Weber_2016_Hey2_cardiomyocite_DESeq2_results1.csv",
    "Woods_2023_Ascl1_neuroblastomaWT31_DESeq2_results1.csv",
    "Woods_2023_Ascl1_neuroblastomaWT51_DESeq2_results1.csv"
)

fc<-data.frame()
i<-1
for (file in fc_files) {
    file_path <- list.files(recursive = TRUE, pattern = paste0("^", file, "$"), full.names = TRUE)
    if (length(file_path) == 0) {
        stop(paste("File not found:", file))
    }
    df <- read.csv(file_path[1])[,c(1,2,3,7)]
    colnames(df) <- c("ensembl","gene_name", "log2FoldChange","padj" )
    
    head(df)
    df$study <- sub("(_DESeq2.*|_microarray.*)", "", file)

    ## join
    if (i == 1) {
        fc <- df
        i <- i + 1
    }
    else {
     fc <- rbind(fc, df)
    }
}



mapping<-c("Aydin_2019_Ascl1_eb_12h" = "ASCL1 -- EB [Aydin et al. 2019]",         
"Aydin_2019_Neurog2_eb_48h"  = "NEUROG2 -- EB [Aydin et al. 2019]",     
"Barfeld_2017_Myc_lncap_5h"  = "MYC -- LNCaP [Barfeld et al. 2017]",          
"Casey_2018_Ascl1_mesc"    = "ASCL1 -- mESC [Casey et al. 2018]",
"Casey_2018_Ascl2_mesc_ESeq2_results1.csv" = "ASCL2 -- mESC [Casey et al. 2018]",
"Casey_2018_Myod1_mesc"  = "MYOD1 -- mESC [Casey et al. 2018]",                  
"Chalamasetty_2014_Msgn1_eb_48h" = "MSGN1 -- EB [Chalamalasetty et al. 2014]",      
"Conerly_2016_Myf5_mef"  = "MYF5 -- MEF [Conerly et al. 2016]",
"Conerly_2016_Myod1_mef"  = "MYOD1 -- MEF [Conerly et al. 2016]",               
"Fong_2012_Myod1_p19"  = "MYOD1 -- P19 [Fong et al. 2012]",                     
"Fong_2012_Myod1_mef"  = "MYOD1 -- MEF [Fong et al. 2012]",                    
"Fong_2012_Neurod2_p19"  = "NEUROD2 -- P19 [Fong et al. 2012]",                  
"Fong_2012_Neurod2_mef" = "NEUROD2 -- MEF [Fong et al. 2012]",                   
"Fong_2015_Neurod2_p19"  = "NEUROD2 -- P19 [Fong et al. 2015]",                   
"Jaenicke_2016_Myc_imec"  = "MYC -- IMEC [Jaenicke et al. 2016]",         
"Joung_2023_Ascl1_hesc_day7"  = "ASCL1 -- hESC [Joung et al. 2023]",             
"Jung_2017_Myc_u2os" = "MYC -- U2OS [Jung et al. 2017]",                      
"Lee_2020_Myod1_mef" = "MYOD1 -- MEF [Lee et al. 2020]",                      
"Li_2019_Twist2_myoblast_DM" = "TWIST2 -- Myoblast (DM) [Li et al. 2019]",              
"Li_2019_Twist2_myoblast_GM" = "TWIST2 -- Myoblast (GM) [Li et al. 2019]",              
"Liu_2014_Ascl2_cd4" = "ASCL2 -- CD4+ T cells [Liu et al. 2023]",                      
"Lorenzin_2016_Myc_u2os" = "MYC -- U2OS [Lorenzin et al. 2016]",                  
"Manandhar_2016_Myod1_fibroblast" = "MYOD1 -- Human fibroblasts [Manandhar et al. 2017]",         
"Matsuda_2016_Neurod1_microglia_Day2" = "NEUROD1 -- Microglia [Matsuda et al. 2016]",     
"Park_2017_Ascl1_glioblastoma" = "ASCL1 -- G523NS [Park et al. 2017]",            
"Pataskar_2016_Neurod1_mesc" = "NEUROD1 -- mESC [Pataskar et al. 2016]",              
"Pereira_2024_Neurog2_astrocyte"  = "NEUROG2 -- Astrocytes [Pereira et al. 2024]",         
"Smith_2016_Neurog2_mrc5" = "NEUROG2 -- MRC5 [Smith et al. 2016]",                 
"Vainorius_2023_Ascl1_mesc" = "ASCL1 -- mESC [Vainorius et al. 2023]",               
"Vainorius_2023_Neurog2_mesc" = "NEUROG2 -- mESC [Vainorius et al. 2023]",             
"Walz_2014_Myc_u2os" = "MYC -- U2OS [Walz et al. 2014]",                      
"Wang_2023_Ascl1_neuroblastoma" = "ASCL1 -- GI-MEN [Wang et al. 2023]",           
"Wapinski_2013_Ascl1_mef" = "ASCL1 -- MEF [Wapinski et al. 2013]",                 
"Weber_2016_Hey1_cardiomyocite" = "HEY1 -- Cardiomiocytes [Weber et al. 2016]",           
"Weber_2016_Hey2_cardiomyocite" = "HEY2 -- Cardiomiocytes [Weber et al. 2016]",           
"Woods_2023_Ascl1_neuroblastomaWT31" = "ASCL1 -- SH-SY5Y [Woods et al. 2023] (1)",      
"Woods_2023_Ascl1_neuroblastomaWT51" = "ASCL1 -- SH-SY5Y [Woods et al. 2023] (2)" )



fc$study <- mapping[fc$study]

### obtain significant genes. genes whose DE is significant in at least one study
significant_genes <- unique(fc$gene_name[fc$padj < 0.05 & !is.na(fc$padj)])






# Find orthologous genes between human and mouse and create a joint correlation heatmap

library(biomaRt)
library(dplyr)
library(tidyr)
library(biomaRt)

library(biomaRt)
library(dplyr)

human <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl"
  # si quieres: mirror = "useast", o version = 110, etc.
)

orthologs <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",                      # símbolo humano genérico
    "mmusculus_homolog_ensembl_gene",         # ID Ensembl ortólogo mouse
    "mmusculus_homolog_associated_gene_name"  # símbolo mouse (MGI)
  ),
  filters = "hgnc_symbol",                    # usamos HGNC SOLO como filtro
  values  = rownames(fc_matrix_human),
  mart    = human
)

orthologs <- orthologs %>%
  dplyr::rename(
    human_ensembl = ensembl_gene_id,
    human_gene    = external_gene_name,
    mouse_ensembl = mmusculus_homolog_ensembl_gene,
    mouse_gene    = mmusculus_homolog_associated_gene_name
  )


# Filter to only genes present in both matrices
orthologs <- orthologs %>%
    filter(human_gene %in% rownames(fc_matrix_human), mouse_gene %in% rownames(fc_matrix_mouse))

# Subset and align matrices by orthologs
fc_matrix_human_ortho <- fc_matrix_human[orthologs$human_gene, , drop = FALSE]
fc_matrix_mouse_ortho <- fc_matrix_mouse[orthologs$mouse_gene, , drop = FALSE]
rownames(fc_matrix_human_ortho) <- orthologs$mouse_gene  # rename to mouse gene names for matching

# Combine matrices (columns = studies, rows = mouse gene names)
fc_matrix_joint <- cbind(fc_matrix_mouse_ortho, fc_matrix_human_ortho)

# Pearson correlation between all pairs of studies (mouse + human)
cor_matrix_joint <- cor(fc_matrix_joint, use = "pairwise.complete.obs", method = "pearson")

pdf("~/Dropbox/induction/7_expression/correlation_joint_human_mouse.pdf", width = 12, height = 10)
pheatmap(
    cor_matrix_joint,
    main = "Pearson correlation between mouse and human studies (orthologous genes)",
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-1, 1, length.out = 101),
    border_color = NA,
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = 45,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    treeheight_row = 0,
    treeheight_col = 0
)
dev.off()


