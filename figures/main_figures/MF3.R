##############################
### Main figure 3a ###########
##############################

plot_selex<-function(factor) {
  neurod2_selex<-read.table(paste0(factor,"_motif_selex.txt"), header = T)
  neurod2_selex_pct<-neurod2_selex[,2:17]/neurod2_selex$total_sequences
  rownames(neurod2_selex_pct)<-c("1","2","3","4")
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CATATG"))]<-"CAT-CAT"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAGCTG"))]<-"CAG-CAG"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CACGTG"))]<-"CAC-CAC"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAATTG"))]<-"CAA-CAA"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAGATG","CATCTG"))]<-"CAT-CAG"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CACATG","CATGTG"))]<-"CAT-CAC"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAAATG","CATTTG"))]<-"CAA-CAT"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAGGTG","CACCTG"))]<-"CAG-CAC"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAACTG","CAGTTG"))]<-"CAA-CAG"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAAGTG","CACTTG"))]<-"CAA-CAC"
  
  #Extract unique names
  nms <- unique(names(neurod2_selex_pct))
  #Sum columns with the same name
  neurod2_selex_pct<-sapply(nms, function(x)  rowSums(neurod2_selex_pct[names(neurod2_selex_pct) %in% x]))
  
  library(reshape2)
  neurod2_selex_pct<-setNames(melt(as.matrix(neurod2_selex_pct)), c('round', 'motif', 'proportion'))
  
  colors<-c("seagreen3", "azure4","yellow", "deeppink", "turquoise2", "chartreuse1","darkgoldenrod4",
                      "red","blue", "black")
                      
  library(ggplot2)
  plot<-ggplot(data = neurod2_selex_pct) + theme_classic() +
    geom_line(aes(x = round, y = proportion, colour = motif, group = motif)) + scale_color_manual(values = colors ) + ggtitle(factor) + theme(legend.position = "none")
  print(plot)
}

plot_selex_yin<-function(factor) {
  neurod2_selex<-as.data.frame(t( read.table(factor, header = T) ) )
  ### If it only contains rounds 3 and 4
  if (dim(neurod2_selex)[1]==3) {
    colnames(neurod2_selex)<-neurod2_selex[1,]
    neurod2_selex<-neurod2_selex[-1,]
    rownames(neurod2_selex)<-c(3:4)
    neurod2_selex<-sapply(neurod2_selex, as.numeric)
    neurod2_selex_pct<-data.matrix(neurod2_selex[,2:17])/as.numeric(neurod2_selex[,1])
    rownames(neurod2_selex_pct)<-c("3","4")
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CATATG"))]<-"CAT-CAT"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAGCTG"))]<-"CAG-CAG"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CACGTG"))]<-"CAC-CAC"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAATTG"))]<-"CAA-CAA"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAGATG","CATCTG"))]<-"CAT-CAG"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CACATG","CATGTG"))]<-"CAT-CAC"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAAATG","CATTTG"))]<-"CAA-CAT"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAGGTG","CACCTG"))]<-"CAG-CAC"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAACTG","CAGTTG"))]<-"CAA-CAG"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAAGTG","CACTTG"))]<-"CAA-CAC"
    
    #Extract unique names
    nms <- unique(colnames(neurod2_selex_pct))
    neurod2_selex_pct<-as.data.frame(neurod2_selex_pct)
    #Sum columns with the same name
    neurod2_selex_pct<-sapply(nms, function(x)  rowSums(neurod2_selex_pct[names(neurod2_selex_pct) %in% x]))
    
    library(reshape2)
    neurod2_selex_pct<-setNames(melt(as.matrix(neurod2_selex_pct)), c('round', 'motif', 'proportion'))
    
    colors<-c("seagreen3", "azure4","yellow", "deeppink", "turquoise2", "chartreuse1","darkgoldenrod4",
                         "red","blue", "black")
                         
    library(ggplot2)
    plot<-ggplot(data = neurod2_selex_pct) + theme_classic() + 
      geom_line(aes(x = round, y = proportion, colour = motif, group = motif)) + scale_color_manual(values = colors ) + ggtitle(factor) + theme(legend.position = "none")
    print(plot)
  }
  else {
    colnames(neurod2_selex)<-neurod2_selex[1,]
    neurod2_selex<-neurod2_selex[-1,]
    rownames(neurod2_selex)<-c(1:4)
    neurod2_selex<-sapply(neurod2_selex, as.numeric)
    neurod2_selex_pct<-data.matrix(neurod2_selex[,2:17])/as.numeric(neurod2_selex[,1])
    rownames(neurod2_selex_pct)<-c("1","2","3","4")
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CATATG"))]<-"CAT-CAT"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAGCTG"))]<-"CAG-CAG"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CACGTG"))]<-"CAC-CAC"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAATTG"))]<-"CAA-CAA"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAGATG","CATCTG"))]<-"CAT-CAG"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CACATG","CATGTG"))]<-"CAT-CAC"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAAATG","CATTTG"))]<-"CAA-CAT"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAGGTG","CACCTG"))]<-"CAG-CAC"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAACTG","CAGTTG"))]<-"CAA-CAG"
    colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAAGTG","CACTTG"))]<-"CAA-CAC"
    
    #Extract unique names
    nms <- unique(colnames(neurod2_selex_pct))
    neurod2_selex_pct<-as.data.frame(neurod2_selex_pct)
    #Sum columns with the same name
    neurod2_selex_pct<-sapply(nms, function(x)  rowSums(neurod2_selex_pct[names(neurod2_selex_pct) %in% x]))
    
    library(reshape2)
    neurod2_selex_pct<-setNames(melt(as.matrix(neurod2_selex_pct)), c('round', 'motif', 'proportion'))
    
    colors<-c("seagreen3", "azure4","yellow", "deeppink", "turquoise2", "chartreuse1","darkgoldenrod4",
                         "red","blue", "black")
                         
    library(ggplot2)
    plot<-ggplot(data = neurod2_selex_pct) + theme_classic() + 
      geom_line(aes(x = round, y = proportion, colour = motif, group = motif)) + scale_color_manual(values = colors ) + ggtitle(factor) + theme(legend.position = "none")
    print(plot)
  }
}




#####################
####### YIN ##########
#####################


# For mycn methyl-selex, they don't find the enriched motif (I checked and they only have cycles 3 and 4).

# To make the logo, they indicate in the supplementary (aaj2239_yin_sm_tables_s1-s6.xlsx) which cycle they use.
# And they compare it with the previous cycle, except when they put a "b". 
# For the ones where my lineplots look wrong, I've seen they use b0.
# That is, against a cycle 0, which I don't know what it is.

# So I discard those from cycle 0 for plotting. 
# Because what I show in the plots are cycles 1 to 4.


setwd("~/htselex/yin")

yin_edbdd<-list.files(path = ".", pattern = "*.motif_counts.txt" )[grep("eDBD",list.files(path = ".", pattern = "*.motif_counts.txt" ))]
yin_fl<-list.files(path = ".", pattern = "*.motif_counts.txt" )[grep("FL",list.files(path = ".", pattern = "*.motif_counts.txt" ))]
yin_fl<-yin_fl[!grepl("eDBD", yin_fl)]

yin_fl_subset<-yin_fl[grep("NEUROD1|NEUROD2|NEUROG2|TWIST2|TWIST1|OLIG2|MSGN1|MESP1|MYOD1|MYF5|ASCL1|ASCL2|TCF4|MYC|MYCN|MITF|HEY1|HEY2|TFE3", yin_fl)]
yin_edbdd_subset<-yin_edbdd[grep("NEUROD1|NEUROD2|NEUROG2|TWIST2|TWIST1|OLIG2|MSGN1|MESP1|MYOD1|MYF5|ASCL1|ASCL2|TCF4|MYC|MYCN|MITF|HEY1|HEY2|TFE3", yin_edbdd)]


yin_fl_max<-yin_fl[grep("MAX", yin_fl)]
yin_edbd_max<-yin_edbdd[grep("MAX", yin_edbdd)]


#(b0 discarded). And also those that only have cycles 3 and 4.

# methyl ht-selex
metil_yin<-unique(c("HEY2*eDBD*TAAGCT40NTGC","HEY2*FL*TAGGAC40NACA","HEY1*FL*TCCTTC40NATT","TFE3*eDBD*TGGTCC40NTTC","MAX*eDBD*TACCGC40NTAC","MAX*FL*TAGCAT40NAGG","OLIG2*eDBD*TTCATT40NTAT","OLIG2*eDBD*TTCATT40NTAT","NEUROG2*eDBD*TTAAAC40NTAT","NEUROD1*FL*TCAGCC40NTTC","MSGN1*eDBD*TATAGC40NGAC","TCF4*eDBD*TTGTCG40NACC","ASCL1*eDBD*TTACGT40NGGA","ASCL2*eDBD*TTGCTG40NCTA","ASCL2*FL*TGCCCG40NCTG","MYOD1*eDBD*TAAGTT40NGAT","MYOD1*eDBD*TAAGTT40NGAT","MYOD1*FL*TGCGAT40NGAT","MYOD1*FL*TGCGAT40NGAT") )
metil_yin<-paste0(metil_yin,".motif_counts.txt")
# ht-selex
nometil_yin<-unique(c("HEY2*eDBD*TCGAGG40NGTT","HEY2*FL*TGCCAA40NTAG","HEY1*FL*TTGCTT40NCCA","TFE3*eDBD*TGGTTA40NTAC","MAX*eDBD*TCTAGG40NCTG","MAX*FL*TTCTGG40NTGC","OLIG2*eDBD*TACTAG40NGGA","OLIG2*eDBD*TACTAG40NGGA","NEUROG2*eDBD*TTTGTA40NTCC","NEUROD1*FL*TTGTAC40NGTA","MSGN1*eDBD*TTCTTG40NTGC","TCF4*eDBD*TGTATA40NGGT","ASCL1*eDBD*TGCGCA40NCAA","ASCL2*eDBD*TAACTC40NGTC","ASCL2*FL*TTCTCC40NAGA","ASCL2*FL*TTCTCC40NAGA","MYOD1*eDBD*TGTCTG40NAGG","MYOD1*eDBD*TGTCTG40NAGG","MYOD1*FL*TACTAA40NATT","MYOD1*FL*TACTAA40NATT") )
nometil_yin<-paste0(nometil_yin,".motif_counts.txt")





######## methylation selected factors ########

# Plot only for HEY1, HEY2, and MAX

selected_factors <- c("HEY1", "HEY2", "MAX")
# Plot methyl and non-methyl side-by-side for selected factors

library(gridExtra)
# Refactor plot_selex_yin to return the ggplot object instead of printing
plot_selex_yin_return <- function(factor) {
  neurod2_selex <- as.data.frame(t(read.table(factor, header = T)))
  if (dim(neurod2_selex)[1] == 3) {
    colnames(neurod2_selex) <- neurod2_selex[1, ]
    neurod2_selex <- neurod2_selex[-1, ]
    rownames(neurod2_selex) <- c(3:4)
    neurod2_selex <- sapply(neurod2_selex, as.numeric)
    neurod2_selex_pct <- data.matrix(neurod2_selex[, 2:17]) / as.numeric(neurod2_selex[, 1])
    rownames(neurod2_selex_pct) <- c("3", "4")
  } else {
    colnames(neurod2_selex) <- neurod2_selex[1, ]
    neurod2_selex <- neurod2_selex[-1, ]
    rownames(neurod2_selex) <- c(1:4)
    neurod2_selex <- sapply(neurod2_selex, as.numeric)
    neurod2_selex_pct <- data.matrix(neurod2_selex[, 2:17]) / as.numeric(neurod2_selex[, 1])
    rownames(neurod2_selex_pct) <- c("1", "2", "3", "4")
  }
  # Rename motifs
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CATATG"))] <- "CAT-CAT"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAGCTG"))] <- "CAG-CAG"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CACGTG"))] <- "CAC-CAC"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAATTG"))] <- "CAA-CAA"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAGATG", "CATCTG"))] <- "CAT-CAG"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CACATG", "CATGTG"))] <- "CAT-CAC"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAAATG", "CATTTG"))] <- "CAA-CAT"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAGGTG", "CACCTG"))] <- "CAG-CAC"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAACTG", "CAGTTG"))] <- "CAA-CAG"
  colnames(neurod2_selex_pct)[which(colnames(neurod2_selex_pct) %in% c("CAAGTG", "CACTTG"))] <- "CAA-CAC"
  # Sum columns with same name
  nms <- unique(colnames(neurod2_selex_pct))
  neurod2_selex_pct <- as.data.frame(neurod2_selex_pct)
  neurod2_selex_pct <- sapply(nms, function(x) rowSums(neurod2_selex_pct[names(neurod2_selex_pct) %in% x]))
  library(reshape2)
  neurod2_selex_pct <- setNames(melt(as.matrix(neurod2_selex_pct)), c('round', 'motif', 'proportion'))
  colors <- c("seagreen3", "azure4", "yellow2", "hotpink1", "hotpink1", "chartreuse1", "goldenrod2", "firebrick2", "blue", "black")
  library(ggplot2)
  plot <- ggplot(data = neurod2_selex_pct) + theme_classic() +
    geom_line(aes(x = round, y = proportion, colour = motif, group = motif)) +
    scale_color_manual(values = colors) +
    ggtitle(factor) +
    theme(legend.position = "right")
  return(plot)
}




pdf("~/Dropbox/induction/4_methylation/methyl_selex_legend.pdf", height = 3, width = 14, useDingbats = F )
for (sel_factor in selected_factors) {
  metil_file <- metil_yin[grepl(sel_factor, metil_yin)]
  nometil_file <- nometil_yin[grepl(sel_factor, nometil_yin)]
  if (length(metil_file) > 0 & length(nometil_file) > 0) {
    p_left <- plot_selex_yin_return(nometil_file[1]) + ggtitle("HT-SELEX")
    p_right <- plot_selex_yin_return(metil_file[1]) + ggtitle("Methyl-HT-SELEX")
    gridExtra::grid.arrange(
      p_left, p_right,
      ncol = 2,
      top = sel_factor
    )
  }
}
dev.off()






#######################################################
#################### Main figure 3b ###################
#######################################################


########### a plot summarizing the mean and dispersion.


combineee <- combineee %>% select(peakid_new, meth, quantile,study) %>% distinct()


combineee <- combineee %>%
  mutate(quantile = case_when(
    quantile <= 25 ~ 4,
    quantile > 25 & quantile <= 50 ~ 3,
    quantile > 50 & quantile <= 75 ~ 2,
    quantile > 75 ~ 1 ))

# Adjust meth for studies where some values are greater than 1
combine_corrected <- combineee %>%
  filter(meth != ".") %>% # Remove rows where meth is "."
  mutate(meth = as.numeric(meth)) 

combine_corrected2<-data.frame()
for (estudio in unique(combine_corrected$study) ) {
  rafa <- combine_corrected %>% filter(study == estudio) 
  if (max(rafa$meth) < 2) { rafa$meth <- rafa$meth * 100 }
  combine_corrected2<-rbind(combine_corrected2, rafa)
}

summary_stats <- combine_corrected2 %>%
  mutate(meth = as.numeric(meth)) %>% # Ensure 'meth' is numeric
  group_by(study, quantile) %>%
  summarise(
    mean_meth = mean(meth, na.rm = TRUE),
    sd_meth = sd(meth, na.rm = TRUE),
    .groups = "drop" # Avoid grouped output warning
  )

study_mapping <- readRDS("study_mapping.RDS")
summary_stats$study <- study_mapping[match(summary_stats$study, study_mapping$study), "renamed_study"]

orden<-readRDS("orden.RDS")
orden<-orden[which(orden %in% summary_stats$study)]
summary_stats$study <- factor(summary_stats$study, levels = rev(orden) )
summary_stats <- summary_stats %>% filter(!is.na(study))

viridis_soft <- rev(viridis(100, option = "G")[1:100] )
library(ggplot2)




### the same but plot a heatmap of the median methylation

combineee <- combinee %>% select(peakid_new, meth, quantile,study) %>% distinct()



# Adjust meth for studies where some values are greater than 1
combine_corrected <- combineee %>%
  filter(meth != ".") %>% # Remove rows where meth is "."
  mutate(meth = as.numeric(meth))
combine_corrected2<-data.frame()
for (estudio in unique(combine_corrected$study) ) {
  rafa <- combine_corrected %>% filter(study == estudio) 
  if (max(rafa$meth) < 2) { rafa$meth <- rafa$meth * 100 }
  combine_corrected2<-rbind(combine_corrected2, rafa)
}



summary_stats <- combine_corrected2 %>%
  mutate(meth = as.numeric(meth)) %>% # Ensure 'meth' is numeric
  group_by(study, quantile) %>%
  summarise(
    median_meth = median(meth, na.rm = TRUE),
    .groups = "drop" # Avoid grouped output warning
  )


orden<-readRDS("orden.RDS")
orden<-orden[which(orden %in% summary_stats$study)]
summary_stats$study <- factor(summary_stats$study, levels = rev(orden) )
summary_stats <- summary_stats %>% filter(!is.na(study))
viridis_soft <- rev(viridis(100, option = "G")[1:100] )
library(ggplot2)



pdf("~/Dropbox/induction/4_methylation/summary_median_meth_quantiles_heatmap.pdf",
    useDingbats = F, height = 7, width = 7)
# Create a heatmap
ggplot(summary_stats, aes(x = quantile, y = study, fill = median_meth)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = viridis_soft, name = "Median % mCpG") +
  theme_minimal() +
  labs(
    title = "Heatmap of Median Methylation Levels",
    x = "Quantile",
    y = "Study"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
dev.off()








#################### Main fig. 3c ####################
######################################################


########## plot only myc and CAT-CAC and CAC-CAC in quantile 4. 

canonical_eboxes_meth_atac<-readRDS("canonical_eboxes_meth_atac.RDS")

canonical_eboxes_meth_atac_myc <- canonical_eboxes_meth_atac[
  grepl("myc|hey1|hey2", canonical_eboxes_meth_atac$study, ignore.case = TRUE),
]

# Convert quantiles to 1-4 (1 = most accessible, 4 = least accessible)
canonical_eboxes_meth_atac_myc$quantile <- dplyr::case_when(
  canonical_eboxes_meth_atac_myc$quantile <= 25 ~ 4,
  canonical_eboxes_meth_atac_myc$quantile > 25 & canonical_eboxes_meth_atac_myc$quantile <= 50 ~ 3,
  canonical_eboxes_meth_atac_myc$quantile > 50 & canonical_eboxes_meth_atac_myc$quantile <= 75 ~ 2,
  canonical_eboxes_meth_atac_myc$quantile > 75 ~ 1
)

canonical_eboxes_meth_atac_4th_quantile<-canonical_eboxes_meth_atac_myc %>% filter(quantile==4)

canonical_eboxes_meth_atac_4th_quantile$meth_window <- factor(canonical_eboxes_meth_atac_4th_quantile$meth_window, levels = c(
  "0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"
))

# get number of peaks (number of unique peakid_new) in each meth window for each study
npeaks_meth_window<- canonical_eboxes_meth_atac_4th_quantile %>%
  select(study, meth_window, peakid_new) %>% distinct() %>% group_by(study,meth_window) %>% summarise(npeaks=n())
 
tabla<-canonical_eboxes_meth_atac_4th_quantile %>% select(central_dinucleotide,study,dist_summit,meth_window)


tabla <- tabla %>%
  group_by(study, central_dinucleotide, dist_summit,meth_window) %>%
  summarise(ocurrencias = n() )


library(ggplot2)
library(ggridges)
library(tidyr)


tabla <- tabla %>%
  mutate(bin = floor(dist_summit / 20) * 20) %>% # Create bins of 20 bp
  group_by(study, central_dinucleotide,meth_window, bin) %>%
  summarise(ocurrencias = sum(ocurrencias, na.rm = TRUE), .groups = "drop") %>%
  rename(dist_summit = bin) 

tabla$dist_summit<- factor(tabla$dist_summit, levels = seq(-400, 400, by = 20)) # Convert to factor for correct order

tabla <- tabla %>%
  group_by(study, central_dinucleotide) %>% # Group by necessary columns
  tidyr::complete(dist_summit,meth_window, fill = list(ocurrencias = 0)) %>% # Fill missing values
  ungroup() # Ungroup after completion



# correct occurences by the number of peaks in the meth windows
tabla <- tabla %>%
  left_join(npeaks_meth_window, by = c("study", "meth_window")) %>%
  mutate(ocurrencias = ocurrencias / npeaks) %>%
  select(-npeaks)


tabla$dist_summit <- as.numeric(as.character(tabla$dist_summit)) # Convert dist_summit to numeric for filtering
tabla <- tabla %>% filter(abs(dist_summit)<=200)

maxii<-tabla %>% group_by(study) %>% summarise(maximo=max(ocurrencias))

tabla<-merge(maxii, tabla, by="study")

tabla$ocurrencias_escaladas<-tabla$ocurrencias/tabla$maximo



# Remove rows with NA central_dinucleotide
tabla <- tabla %>% filter(!is.na(central_dinucleotide) )


study_mapping <- readRDS("study_mapping.RDS")

tabla <- tabla %>%
  mutate(study = study_mapping[match(study, study_mapping$study), "renamed_study"])




orden<-readRDS("orden.RDS")


tabla <- tabla %>%
  mutate(study = factor(study, levels = orden))


tabla<-tabla %>% filter(central_dinucleotide %in% c("CAT-CAC","CAC-CAC")  )



# Define colors for the central dinucleotides
dinucleotide_colors <- c(
   "CAC-CAC" = "goldenrod2", "CAT-CAC" = "darkorange1"
)

library(patchwork)

# Create a list to store plots for each study
plot_list <- list()

for (estudio in unique(tabla$study)) {
  tabla_estudio <- tabla %>% filter(study == estudio)

  gg <- ggplot(tabla_estudio, aes(x = dist_summit, y = ocurrencias_escaladas, group = meth_window, fill = central_dinucleotide)) +
    geom_area(alpha = 0.7, color = "black") +
    facet_grid(meth_window ~ central_dinucleotide, scales = "fixed") +
    scale_fill_manual(values = dinucleotide_colors) +
    theme_void() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.y = element_text(angle = 0, size = 10),
      panel.spacing.y = unit(-4, "lines"),
      legend.position = "bottom"
    ) +
    labs(
      title = estudio,
      x = "Distance to summit (binned by 20 bp)",
      y = "Frequency",
      fill = "Central dinucleotide"
    ) + ylim(0, 1)
  plot_list[[estudio]] <- gg
}

# Combine all plots side by side
combined_plot <- wrap_plots(plot_list, nrow = 1, guides = "collect")

pdf("~/Dropbox/induction/4_methylation/caccac_catcac_myc_4th_quantile_meth_windows.pdf", height = 13, width = 5 * length(plot_list), useDingbats = F)
  print(combined_plot)
dev.off()
