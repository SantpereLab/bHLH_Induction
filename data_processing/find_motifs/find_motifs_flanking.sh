
cd ~/proneural/chip-seq



## find flanking mucleotides in mouse

echo '
find_noncanonical_mouse <- function(peaks, peaks_file, motif) {
    require(dplyr)
    require(valr)
    require(pracma)
    
    system(paste0("bedtools getfasta -fi ~/genomes/mm10.fa -bed ", peaks_file, " | grep -v chr > ", peaks_file, ".fa")) # nolint
    seqs <- read.table(paste0(peaks_file, ".fa"))
    seqs$V1 <- toupper(seqs$V1)
    
    ebox <- data.frame()
    for (i in seq_along(seqs$V1)) {
        ebox_chr <- peaks$chr[i]
        matches <- refindall(seqs$V1[i], motif)
        if (length(matches) == 0) next
        
        for (j in seq_along(matches)) {
            ebox <- rbind(ebox, data.frame(
                chrom = ebox_chr,
                start = peaks$start[i] + matches[j],
                end = peaks$start[i] + matches[j] + 8,
                ebox_seq = substr(seqs$V1[i], matches[j], matches[j] + 7)
            ))
        }
    }
    
    colnames(peaks)[1] <- "chrom"
    ebox <- ebox %>% mutate(start = as.numeric(start), end = as.numeric(end))
    inter <- valr::bed_intersect(peaks, ebox) %>%
        select(chrom, start = start.x, end = end.x, summit = summit.x, peak_length = peak_length.x, peak_id = peak_id.x, ebox_start = start.y, ebox_seq = ebox_seq.y, real_start = real_start.x, real_end = real_end.x)
    
    return(inter)
}

find_noncanonical_mouse_extended <- function(peaks, peakfile, motif) {
    require(dplyr)
    
    peaks$summit<-peaks$V10+peaks$V2
    peaks<-peaks %>% select(V1,V2,V3,summit)
    colnames(peaks)<-c("chr","start","end","summit")

    # add natural length of the peak
    peaks$real_start <- peaks$start
    peaks$real_end <- peaks$end
    peaks <- peaks %>% filter(chr != "chrM")
    npeaks <- nrow(peaks)
    peaks <- peaks %>% mutate(
        peak_id = paste0("peak_", row_number()),
        peak_length = end - start,
        start = summit - 400,
        end = summit + 400
    )
    
    motifs <- do.call(rbind, lapply(seq(1, npeaks, by = 2000), function(i) {
        peaks1 <- peaks[i:min(i + 1999, npeaks), ]
        options(scipen = 999)
        write.table(peaks1, file = paste0("rep_", ceiling(i / 2000), "_", peakfile), sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
        find_noncanonical_mouse(peaks = peaks1, peaks_file = paste0("rep_", ceiling(i / 2000), "_", peakfile), motif)
    } ) )

    # convert the sequence of E-boxes + flanking nucleotide. Into strand-oriented half-sites + the flanking nucleotides.
    # in other words, maintain the first 4 nucleotides. and reverse-complement the last 4 nucleotides. and separate by "-"
    


    motifs$flanking <- NA
    # CATATG eboxes
    motifs$flanking[motifs$ebox_seq == "ACATATGA"] <- "ACAT-TCAT"
    motifs$flanking[motifs$ebox_seq == "ACATATGT"] <- "ACAT-ACAT"
    motifs$flanking[motifs$ebox_seq == "ACATATGG"] <- "ACAT-CCAT"
    motifs$flanking[motifs$ebox_seq == "ACATATGC"] <- "ACAT-GCAT"
    
    motifs$flanking[motifs$ebox_seq == "TCATATGA"] <- "TCAT-TCAT"
    motifs$flanking[motifs$ebox_seq == "TCATATGT"] <- "TCAT-ACAT"
    motifs$flanking[motifs$ebox_seq == "TCATATGG"] <- "TCAT-CCAT"
    motifs$flanking[motifs$ebox_seq == "TCATATGC"] <- "TCAT-GCAT"

    motifs$flanking[motifs$ebox_seq == "GCATATGA"] <- "GCAT-TCAT"
    motifs$flanking[motifs$ebox_seq == "GCATATGT"] <- "GCAT-ACAT"
    motifs$flanking[motifs$ebox_seq == "GCATATGG"] <- "GCAT-CCAT"
    motifs$flanking[motifs$ebox_seq == "GCATATGC"] <- "GCAT-GCAT"

    motifs$flanking[motifs$ebox_seq == "CCATATGA"] <- "CCAT-TCAT"
    motifs$flanking[motifs$ebox_seq == "CCATATGT"] <- "CCAT-ACAT"
    motifs$flanking[motifs$ebox_seq == "CCATATGG"] <- "CCAT-CCAT"
    motifs$flanking[motifs$ebox_seq == "CCATATGC"] <- "CCAT-GCAT"

    # CAGATG eboxes
    motifs$flanking[motifs$ebox_seq == "ACAGATGA"] <- "ACAG-TCAT"
    motifs$flanking[motifs$ebox_seq == "ACAGATGT"] <- "ACAG-ACAT"
    motifs$flanking[motifs$ebox_seq == "ACAGATGG"] <- "ACAG-CCAT"
    motifs$flanking[motifs$ebox_seq == "ACAGATGC"] <- "ACAG-GCAT"

    motifs$flanking[motifs$ebox_seq == "TCAGATGA"] <- "TCAG-TCAT"
    motifs$flanking[motifs$ebox_seq == "TCAGATGT"] <- "TCAG-ACAT"
    motifs$flanking[motifs$ebox_seq == "TCAGATGG"] <- "TCAG-CCAT"
    motifs$flanking[motifs$ebox_seq == "TCAGATGC"] <- "TCAG-GCAT"

    motifs$flanking[motifs$ebox_seq == "GCAGATGA"] <- "GCAG-TCAT"
    motifs$flanking[motifs$ebox_seq == "GCAGATGT"] <- "GCAG-ACAT"
    motifs$flanking[motifs$ebox_seq == "GCAGATGG"] <- "GCAG-CCAT"
    motifs$flanking[motifs$ebox_seq == "GCAGATGC"] <- "GCAG-GCAT"

    motifs$flanking[motifs$ebox_seq == "CCAGATGA"] <- "CCAG-TCAT"
    motifs$flanking[motifs$ebox_seq == "CCAGATGT"] <- "CCAG-ACAT"
    motifs$flanking[motifs$ebox_seq == "CCAGATGG"] <- "CCAG-CCAT"
    motifs$flanking[motifs$ebox_seq == "CCAGATGC"] <- "CCAG-GCAT"

    # CATCTG eboxes
    motifs$flanking[motifs$ebox_seq == "ACATCTGA"] <- "ACAT-TCAG"
    motifs$flanking[motifs$ebox_seq == "ACATCTGT"] <- "ACAT-ACAG"
    motifs$flanking[motifs$ebox_seq == "ACATCTGG"] <- "ACAT-CCAG"
    motifs$flanking[motifs$ebox_seq == "ACATCTGC"] <- "ACAT-GCAG"

    motifs$flanking[motifs$ebox_seq == "TCATCTGA"] <- "TCAT-TCAG"
    motifs$flanking[motifs$ebox_seq == "TCATCTGT"] <- "TCAT-ACAG"
    motifs$flanking[motifs$ebox_seq == "TCATCTGG"] <- "TCAT-CCAG"
    motifs$flanking[motifs$ebox_seq == "TCATCTGC"] <- "TCAT-GCAG"

    motifs$flanking[motifs$ebox_seq == "GCATCTGA"] <- "GCAT-TCAG"
    motifs$flanking[motifs$ebox_seq == "GCATCTGT"] <- "GCAT-ACAG"
    motifs$flanking[motifs$ebox_seq == "GCATCTGG"] <- "GCAT-CCAG"
    motifs$flanking[motifs$ebox_seq == "GCATCTGC"] <- "GCAT-GCAG"

    motifs$flanking[motifs$ebox_seq == "CCATCTGA"] <- "CCAT-TCAG"
    motifs$flanking[motifs$ebox_seq == "CCATCTGT"] <- "CCAT-ACAG"
    motifs$flanking[motifs$ebox_seq == "CCATCTGG"] <- "CCAT-CCAG"
    motifs$flanking[motifs$ebox_seq == "CCATCTGC"] <- "CCAT-GCAG"

    # CAGCTG eboxes
    motifs$flanking[motifs$ebox_seq == "ACAGCTGA"] <- "ACAG-TCAG"
    motifs$flanking[motifs$ebox_seq == "ACAGCTGT"] <- "ACAG-ACAG"
    motifs$flanking[motifs$ebox_seq == "ACAGCTGG"] <- "ACAG-CCAG"
    motifs$flanking[motifs$ebox_seq == "ACAGCTGC"] <- "ACAG-GCAG"

    motifs$flanking[motifs$ebox_seq == "TCAGCTGA"] <- "TCAG-TCAG"
    motifs$flanking[motifs$ebox_seq == "TCAGCTGT"] <- "TCAG-ACAG"
    motifs$flanking[motifs$ebox_seq == "TCAGCTGG"] <- "TCAG-CCAG"
    motifs$flanking[motifs$ebox_seq == "TCAGCTGC"] <- "TCAG-GCAG"

    motifs$flanking[motifs$ebox_seq == "GCAGCTGA"] <- "GCAG-TCAG"
    motifs$flanking[motifs$ebox_seq == "GCAGCTGT"] <- "GCAG-ACAG"
    motifs$flanking[motifs$ebox_seq == "GCAGCTGG"] <- "GCAG-CCAG"
    motifs$flanking[motifs$ebox_seq == "GCAGCTGC"] <- "GCAG-GCAG"

    motifs$flanking[motifs$ebox_seq == "CCAGCTGA"] <- "CCAG-TCAG"
    motifs$flanking[motifs$ebox_seq == "CCAGCTGT"] <- "CCAG-ACAG"
    motifs$flanking[motifs$ebox_seq == "CCAGCTGG"] <- "CCAG-CCAG"
    motifs$flanking[motifs$ebox_seq == "CCAGCTGC"] <- "CCAG-GCAG"


    # CAGGTG eboxes
    motifs$flanking[motifs$ebox_seq == "ACAGGTGA"] <- "ACAG-TCAC"
    motifs$flanking[motifs$ebox_seq == "ACAGGTGT"] <- "ACAG-ACAC"
    motifs$flanking[motifs$ebox_seq == "ACAGGTGG"] <- "ACAG-CCAC"
    motifs$flanking[motifs$ebox_seq == "ACAGGTGC"] <- "ACAG-GCAC"

    motifs$flanking[motifs$ebox_seq == "TCAGGTGA"] <- "TCAG-TCAC"
    motifs$flanking[motifs$ebox_seq == "TCAGGTGT"] <- "TCAG-ACAC"
    motifs$flanking[motifs$ebox_seq == "TCAGGTGG"] <- "TCAG-CCAC"
    motifs$flanking[motifs$ebox_seq == "TCAGGTGC"] <- "TCAG-GCAC"

    motifs$flanking[motifs$ebox_seq == "GCAGGTGA"] <- "GCAG-TCAC"
    motifs$flanking[motifs$ebox_seq == "GCAGGTGT"] <- "GCAG-ACAC"
    motifs$flanking[motifs$ebox_seq == "GCAGGTGG"] <- "GCAG-CCAC"
    motifs$flanking[motifs$ebox_seq == "GCAGGTGC"] <- "GCAG-GCAC"

    motifs$flanking[motifs$ebox_seq == "CCAGGTGA"] <- "CCAG-TCAC"
    motifs$flanking[motifs$ebox_seq == "CCAGGTGT"] <- "CCAG-ACAC"
    motifs$flanking[motifs$ebox_seq == "CCAGGTGG"] <- "CCAG-CCAC"
    motifs$flanking[motifs$ebox_seq == "CCAGGTGC"] <- "CCAG-GCAC"

    # CACCTG eboxes
    motifs$flanking[motifs$ebox_seq == "ACACCTGA"] <- "ACAC-TCAG"
    motifs$flanking[motifs$ebox_seq == "ACACCTGT"] <- "ACAC-ACAG"
    motifs$flanking[motifs$ebox_seq == "ACACCTGG"] <- "ACAC-CCAG"
    motifs$flanking[motifs$ebox_seq == "ACACCTGC"] <- "ACAC-GCAG"

    motifs$flanking[motifs$ebox_seq == "TCACCTGA"] <- "TCAC-TCAG"
    motifs$flanking[motifs$ebox_seq == "TCACCTGT"] <- "TCAC-ACAG"
    motifs$flanking[motifs$ebox_seq == "TCACCTGG"] <- "TCAC-CCAG"
    motifs$flanking[motifs$ebox_seq == "TCACCTGC"] <- "TCAC-GCAG"

    motifs$flanking[motifs$ebox_seq == "GCACCTGA"] <- "GCAC-TCAG"
    motifs$flanking[motifs$ebox_seq == "GCACCTGT"] <- "GCAC-ACAG"
    motifs$flanking[motifs$ebox_seq == "GCACCTGG"] <- "GCAC-CCAG"
    motifs$flanking[motifs$ebox_seq == "GCACCTGC"] <- "GCAC-GCAG"

    motifs$flanking[motifs$ebox_seq == "CCACCTGA"] <- "CCAC-TCAG"
    motifs$flanking[motifs$ebox_seq == "CCACCTGT"] <- "CCAC-ACAG"
    motifs$flanking[motifs$ebox_seq == "CCACCTGG"] <- "CCAC-CCAG"
    motifs$flanking[motifs$ebox_seq == "CCACCTGC"] <- "CCAC-GCAG"

    # CATGTG eboxes
    motifs$flanking[motifs$ebox_seq == "ACATGTGA"] <- "ACAT-TCAC"
    motifs$flanking[motifs$ebox_seq == "ACATGTGT"] <- "ACAT-ACAC"
    motifs$flanking[motifs$ebox_seq == "ACATGTGG"] <- "ACAT-CCAC"
    motifs$flanking[motifs$ebox_seq == "ACATGTGC"] <- "ACAT-GCAC"

    motifs$flanking[motifs$ebox_seq == "TCATGTGA"] <- "TCAT-TCAC"
    motifs$flanking[motifs$ebox_seq == "TCATGTGT"] <- "TCAT-ACAC"
    motifs$flanking[motifs$ebox_seq == "TCATGTGG"] <- "TCAT-CCAC"
    motifs$flanking[motifs$ebox_seq == "TCATGTGC"] <- "TCAT-GCAC"

    motifs$flanking[motifs$ebox_seq == "GCATGTGA"] <- "GCAT-TCAC"
    motifs$flanking[motifs$ebox_seq == "GCATGTGT"] <- "GCAT-ACAC"
    motifs$flanking[motifs$ebox_seq == "GCATGTGG"] <- "GCAT-CCAC"
    motifs$flanking[motifs$ebox_seq == "GCATGTGC"] <- "GCAT-GCAC"

    motifs$flanking[motifs$ebox_seq == "CCATGTGA"] <- "CCAT-TCAC"
    motifs$flanking[motifs$ebox_seq == "CCATGTGT"] <- "CCAT-ACAC"
    motifs$flanking[motifs$ebox_seq == "CCATGTGG"] <- "CCAT-CCAC"
    motifs$flanking[motifs$ebox_seq == "CCATGTGC"] <- "CCAT-GCAC"

    # CACATG eboxes

    motifs$flanking[motifs$ebox_seq == "ACACATGA"] <- "ACAC-TCAT"
    motifs$flanking[motifs$ebox_seq == "ACACATGT"] <- "ACAC-ACAT"
    motifs$flanking[motifs$ebox_seq == "ACACATGG"] <- "ACAC-CCAT"
    motifs$flanking[motifs$ebox_seq == "ACACATGC"] <- "ACAC-GCAT"

    motifs$flanking[motifs$ebox_seq == "TCACATGA"] <- "TCAC-TCAT"
    motifs$flanking[motifs$ebox_seq == "TCACATGT"] <- "TCAC-ACAT"
    motifs$flanking[motifs$ebox_seq == "TCACATGG"] <- "TCAC-CCAT"
    motifs$flanking[motifs$ebox_seq == "TCACATGC"] <- "TCAC-GCAT"

    motifs$flanking[motifs$ebox_seq == "GCACATGA"] <- "GCAC-TCAT"
    motifs$flanking[motifs$ebox_seq == "GCACATGT"] <- "GCAC-ACAT"
    motifs$flanking[motifs$ebox_seq == "GCACATGG"] <- "GCAC-CCAT"
    motifs$flanking[motifs$ebox_seq == "GCACATGC"] <- "GCAC-GCAT"

    motifs$flanking[motifs$ebox_seq == "CCACATGA"] <- "CCAC-TCAT"
    motifs$flanking[motifs$ebox_seq == "CCACATGT"] <- "CCAC-ACAT"
    motifs$flanking[motifs$ebox_seq == "CCACATGG"] <- "CCAC-CCAT"
    motifs$flanking[motifs$ebox_seq == "CCACATGC"] <- "CCAC-GCAT"

    # CACGTG eboxes
    motifs$flanking[motifs$ebox_seq == "ACACGTGA"] <- "ACAC-TCAC"
    motifs$flanking[motifs$ebox_seq == "ACACGTGT"] <- "ACAC-ACAC"
    motifs$flanking[motifs$ebox_seq == "ACACGTGG"] <- "ACAC-CCAC"
    motifs$flanking[motifs$ebox_seq == "ACACGTGC"] <- "ACAC-GCAC"

    motifs$flanking[motifs$ebox_seq == "TCACGTGA"] <- "TCAC-TCAC"
    motifs$flanking[motifs$ebox_seq == "TCACGTGT"] <- "TCAC-ACAC"
    motifs$flanking[motifs$ebox_seq == "TCACGTGG"] <- "TCAC-CCAC"
    motifs$flanking[motifs$ebox_seq == "TCACGTGC"] <- "TCAC-GCAC"

    motifs$flanking[motifs$ebox_seq == "GCACGTGA"] <- "GCAC-TCAC"
    motifs$flanking[motifs$ebox_seq == "GCACGTGT"] <- "GCAC-ACAC"
    motifs$flanking[motifs$ebox_seq == "GCACGTGG"] <- "GCAC-CCAC"
    motifs$flanking[motifs$ebox_seq == "GCACGTGC"] <- "GCAC-GCAC"

    motifs$flanking[motifs$ebox_seq == "CCACGTGA"] <- "CCAC-TCAC"
    motifs$flanking[motifs$ebox_seq == "CCACGTGT"] <- "CCAC-ACAC"
    motifs$flanking[motifs$ebox_seq == "CCACGTGG"] <- "CCAC-CCAC"
    motifs$flanking[motifs$ebox_seq == "CCACGTGC"] <- "CCAC-GCAC"

    # CAATTG eboxes
    motifs$flanking[motifs$ebox_seq == "ACAATTGA"] <- "ACAA-TCAA"
    motifs$flanking[motifs$ebox_seq == "ACAATTGT"] <- "ACAA-ACAA"
    motifs$flanking[motifs$ebox_seq == "ACAATTGG"] <- "ACAA-CCAA"
    motifs$flanking[motifs$ebox_seq == "ACAATTGC"] <- "ACAA-GCAA"

    motifs$flanking[motifs$ebox_seq == "TCAATTGA"] <- "TCAA-TCAA"
    motifs$flanking[motifs$ebox_seq == "TCAATTGT"] <- "TCAA-ACAA"
    motifs$flanking[motifs$ebox_seq == "TCAATTGG"] <- "TCAA-CCAA"
    motifs$flanking[motifs$ebox_seq == "TCAATTGC"] <- "TCAA-GCAA"

    motifs$flanking[motifs$ebox_seq == "GCAATTGA"] <- "GCAA-TCAA"
    motifs$flanking[motifs$ebox_seq == "GCAATTGT"] <- "GCAA-ACAA"
    motifs$flanking[motifs$ebox_seq == "GCAATTGG"] <- "GCAA-CCAA"
    motifs$flanking[motifs$ebox_seq == "GCAATTGC"] <- "GCAA-GCAA"

    motifs$flanking[motifs$ebox_seq == "CCAATTGA"] <- "CCAA-TCAA"
    motifs$flanking[motifs$ebox_seq == "CCAATTGT"] <- "CCAA-ACAA"
    motifs$flanking[motifs$ebox_seq == "CCAATTGG"] <- "CCAA-CCAA"
    motifs$flanking[motifs$ebox_seq == "CCAATTGC"] <- "CCAA-GCAA"

    # CAAATG eboxes
    motifs$flanking[motifs$ebox_seq == "ACAAATGA"] <- "ACAA-TCAT"
    motifs$flanking[motifs$ebox_seq == "ACAAATGT"] <- "ACAA-ACAT"
    motifs$flanking[motifs$ebox_seq == "ACAAATGG"] <- "ACAA-CCAT"
    motifs$flanking[motifs$ebox_seq == "ACAAATGC"] <- "ACAA-GCAT"

    motifs$flanking[motifs$ebox_seq == "TCAAATGA"] <- "TCAA-TCAT"
    motifs$flanking[motifs$ebox_seq == "TCAAATGT"] <- "TCAA-ACAT"
    motifs$flanking[motifs$ebox_seq == "TCAAATGG"] <- "TCAA-CCAT"
    motifs$flanking[motifs$ebox_seq == "TCAAATGC"] <- "TCAA-GCAT"

    motifs$flanking[motifs$ebox_seq == "GCAAATGA"] <- "GCAA-TCAT"
    motifs$flanking[motifs$ebox_seq == "GCAAATGT"] <- "GCAA-ACAT"
    motifs$flanking[motifs$ebox_seq == "GCAAATGG"] <- "GCAA-CCAT"
    motifs$flanking[motifs$ebox_seq == "GCAAATGC"] <- "GCAA-GCAT"

    motifs$flanking[motifs$ebox_seq == "CCAAATGA"] <- "CCAA-TCAT"
    motifs$flanking[motifs$ebox_seq == "CCAAATGT"] <- "CCAA-ACAT"
    motifs$flanking[motifs$ebox_seq == "CCAAATGG"] <- "CCAA-CCAT"
    motifs$flanking[motifs$ebox_seq == "CCAAATGC"] <- "CCAA-GCAT"

    # CAACTG eboxes
    motifs$flanking[motifs$ebox_seq == "ACAACTGA"] <- "ACAA-TCAG"
    motifs$flanking[motifs$ebox_seq == "ACAACTGT"] <- "ACAA-ACAG"
    motifs$flanking[motifs$ebox_seq == "ACAACTGG"] <- "ACAA-CCAG"
    motifs$flanking[motifs$ebox_seq == "ACAACTGC"] <- "ACAA-GCAG"

    motifs$flanking[motifs$ebox_seq == "TCAACTGA"] <- "TCAA-TCAG"
    motifs$flanking[motifs$ebox_seq == "TCAACTGT"] <- "TCAA-ACAG"
    motifs$flanking[motifs$ebox_seq == "TCAACTGG"] <- "TCAA-CCAG"
    motifs$flanking[motifs$ebox_seq == "TCAACTGC"] <- "TCAA-GCAG"

    motifs$flanking[motifs$ebox_seq == "GCAACTGA"] <- "GCAA-TCAG"
    motifs$flanking[motifs$ebox_seq == "GCAACTGT"] <- "GCAA-ACAG"
    motifs$flanking[motifs$ebox_seq == "GCAACTGG"] <- "GCAA-CCAG"
    motifs$flanking[motifs$ebox_seq == "GCAACTGC"] <- "GCAA-GCAG"

    motifs$flanking[motifs$ebox_seq == "CCAACTGA"] <- "CCAA-TCAG"
    motifs$flanking[motifs$ebox_seq == "CCAACTGT"] <- "CCAA-ACAG"
    motifs$flanking[motifs$ebox_seq == "CCAACTGG"] <- "CCAA-CCAG"
    motifs$flanking[motifs$ebox_seq == "CCAACTGC"] <- "CCAA-GCAG"


    # CAAGTG eboxes
    motifs$flanking[motifs$ebox_seq == "ACAAGTGA"] <- "ACAA-TCAC"
    motifs$flanking[motifs$ebox_seq == "ACAAGTGT"] <- "ACAA-ACAC"
    motifs$flanking[motifs$ebox_seq == "ACAAGTGG"] <- "ACAA-CCAC"
    motifs$flanking[motifs$ebox_seq == "ACAAGTGC"] <- "ACAA-GCAC"

    motifs$flanking[motifs$ebox_seq == "TCAAGTGA"] <- "TCAA-TCAC"
    motifs$flanking[motifs$ebox_seq == "TCAAGTGT"] <- "TCAA-ACAC"
    motifs$flanking[motifs$ebox_seq == "TCAAGTGG"] <- "TCAA-CCAC"
    motifs$flanking[motifs$ebox_seq == "TCAAGTGC"] <- "TCAA-GCAC"

    motifs$flanking[motifs$ebox_seq == "GCAAGTGA"] <- "GCAA-TCAC"
    motifs$flanking[motifs$ebox_seq == "GCAAGTGT"] <- "GCAA-ACAC"
    motifs$flanking[motifs$ebox_seq == "GCAAGTGG"] <- "GCAA-CCAC"
    motifs$flanking[motifs$ebox_seq == "GCAAGTGC"] <- "GCAA-GCAC"

    motifs$flanking[motifs$ebox_seq == "CCAAGTGA"] <- "CCAA-TCAC"
    motifs$flanking[motifs$ebox_seq == "CCAAGTGT"] <- "CCAA-ACAC"
    motifs$flanking[motifs$ebox_seq == "CCAAGTGG"] <- "CCAA-CCAC"
    motifs$flanking[motifs$ebox_seq == "CCAAGTGC"] <- "CCAA-GCAC"

    # CAGTTG eboxes
    motifs$flanking[motifs$ebox_seq == "ACAGTTGA"] <- "ACAG-TCAA"
    motifs$flanking[motifs$ebox_seq == "ACAGTTGT"] <- "ACAG-ACAA"
    motifs$flanking[motifs$ebox_seq == "ACAGTTGG"] <- "ACAG-CCAA"
    motifs$flanking[motifs$ebox_seq == "ACAGTTGC"] <- "ACAG-GCAA"

    motifs$flanking[motifs$ebox_seq == "TCAGTTGA"] <- "TCAG-TCAA"
    motifs$flanking[motifs$ebox_seq == "TCAGTTGT"] <- "TCAG-ACAA"
    motifs$flanking[motifs$ebox_seq == "TCAGTTGG"] <- "TCAG-CCAA"
    motifs$flanking[motifs$ebox_seq == "TCAGTTGC"] <- "TCAG-GCAA"

    motifs$flanking[motifs$ebox_seq == "GCAGTTGA"] <- "GCAG-TCAA"
    motifs$flanking[motifs$ebox_seq == "GCAGTTGT"] <- "GCAG-ACAA"
    motifs$flanking[motifs$ebox_seq == "GCAGTTGG"] <- "GCAG-CCAA"
    motifs$flanking[motifs$ebox_seq == "GCAGTTGC"] <- "GCAG-GCAA"

    motifs$flanking[motifs$ebox_seq == "CCAGTTGA"] <- "CCAG-TCAA"
    motifs$flanking[motifs$ebox_seq == "CCAGTTGT"] <- "CCAG-ACAA"
    motifs$flanking[motifs$ebox_seq == "CCAGTTGG"] <- "CCAG-CCAA"
    motifs$flanking[motifs$ebox_seq == "CCAGTTGC"] <- "CCAG-GCAA"

    # CACTTG eboxes
    motifs$flanking[motifs$ebox_seq == "ACACTTGA"] <- "ACAC-TCAA"
    motifs$flanking[motifs$ebox_seq == "ACACTTGT"] <- "ACAC-ACAA"
    motifs$flanking[motifs$ebox_seq == "ACACTTGG"] <- "ACAC-CCAA"
    motifs$flanking[motifs$ebox_seq == "ACACTTGC"] <- "ACAC-GCAA"

    motifs$flanking[motifs$ebox_seq == "TCACTTGA"] <- "TCAC-TCAA"
    motifs$flanking[motifs$ebox_seq == "TCACTTGT"] <- "TCAC-ACAA"
    motifs$flanking[motifs$ebox_seq == "TCACTTGG"] <- "TCAC-CCAA"
    motifs$flanking[motifs$ebox_seq == "TCACTTGC"] <- "TCAC-GCAA"

    motifs$flanking[motifs$ebox_seq == "GCACTTGA"] <- "GCAC-TCAA"
    motifs$flanking[motifs$ebox_seq == "GCACTTGT"] <- "GCAC-ACAA"
    motifs$flanking[motifs$ebox_seq == "GCACTTGG"] <- "GCAC-CCAA"
    motifs$flanking[motifs$ebox_seq == "GCACTTGC"] <- "GCAC-GCAA"

    motifs$flanking[motifs$ebox_seq == "CCACTTGA"] <- "CCAC-TCAA"
    motifs$flanking[motifs$ebox_seq == "CCACTTGT"] <- "CCAC-ACAA"
    motifs$flanking[motifs$ebox_seq == "CCACTTGG"] <- "CCAC-CCAA"
    motifs$flanking[motifs$ebox_seq == "CCACTTGC"] <- "CCAC-GCAA"


    # CATTTG eboxes
    motifs$flanking[motifs$ebox_seq == "ACATTTGA"] <- "ACAT-TCAA"
    motifs$flanking[motifs$ebox_seq == "ACATTTGT"] <- "ACAT-ACAA"
    motifs$flanking[motifs$ebox_seq == "ACATTTGG"] <- "ACAT-CCAA"
    motifs$flanking[motifs$ebox_seq == "ACATTTGC"] <- "ACAT-GCAA"

    motifs$flanking[motifs$ebox_seq == "TCATTTGA"] <- "TCAT-TCAA"
    motifs$flanking[motifs$ebox_seq == "TCATTTGT"] <- "TCAT-ACAA"
    motifs$flanking[motifs$ebox_seq == "TCATTTGG"] <- "TCAT-CCAA"
    motifs$flanking[motifs$ebox_seq == "TCATTTGC"] <- "TCAT-GCAA"

    motifs$flanking[motifs$ebox_seq == "GCATTTGA"] <- "GCAT-TCAA"
    motifs$flanking[motifs$ebox_seq == "GCATTTGT"] <- "GCAT-ACAA"
    motifs$flanking[motifs$ebox_seq == "GCATTTGG"] <- "GCAT-CCAA"
    motifs$flanking[motifs$ebox_seq == "GCATTTGC"] <- "GCAT-GCAA"

    motifs$flanking[motifs$ebox_seq == "CCATTTGA"] <- "CCAT-TCAA"
    motifs$flanking[motifs$ebox_seq == "CCATTTGT"] <- "CCAT-ACAA"
    motifs$flanking[motifs$ebox_seq == "CCATTTGG"] <- "CCAT-CCAA"
    motifs$flanking[motifs$ebox_seq == "CCATTTGC"] <- "CCAT-GCAA"

    return(motifs)
}


args <- commandArgs(trailingOnly = TRUE)
file <- args[1]

motif <- c(".CA..TG." )
peaks <- read.table(file, header = FALSE)
all_eboxes <- find_noncanonical_mouse_extended(peaks, file, motif)

estudio<-gsub("_confidence.bed", "", file)
all_eboxes$study<-estudio
saveRDS(all_eboxes,paste0(estudio,"_flanking.RDS") )
' > find_eboxes_mouse_flanking.R




# the same script in human


echo '
find_noncanonical_human <- function(peaks, peaks_file, motif) {
    require(dplyr)
    require(valr)
    require(pracma)
    
    system(paste0("bedtools getfasta -fi ~/genomes/hg38.fa -bed ", peaks_file, " | grep -v chr > ", peaks_file, ".fa")) # nolint
    seqs <- read.table(paste0(peaks_file, ".fa"))
    seqs$V1 <- toupper(seqs$V1)
    
    ebox <- data.frame()
    for (i in seq_along(seqs$V1)) {
        ebox_chr <- peaks$chr[i]
        matches <- refindall(seqs$V1[i], motif)
        if (length(matches) == 0) next
        
        for (j in seq_along(matches)) {
            ebox <- rbind(ebox, data.frame(
                chrom = ebox_chr,
                start = peaks$start[i] + matches[j],
                end = peaks$start[i] + matches[j] + 8,
                ebox_seq = substr(seqs$V1[i], matches[j], matches[j] + 7)
            ))
        }
    }
    
    colnames(peaks)[1] <- "chrom"
    ebox <- ebox %>% mutate(start = as.numeric(start), end = as.numeric(end))
    inter <- valr::bed_intersect(peaks, ebox) %>%
        select(chrom, start = start.x, end = end.x, summit = summit.x, peak_length = peak_length.x, peak_id = peak_id.x, ebox_start = start.y, ebox_seq = ebox_seq.y, real_start = real_start.x, real_end = real_end.x)
    
    return(inter)
}

find_noncanonical_human_extended <- function(peaks, peakfile, motif) {
    require(dplyr)
    
    peaks$summit<-peaks$V10+peaks$V2
    peaks<-peaks %>% select(V1,V2,V3,summit)
    colnames(peaks)<-c("chr","start","end","summit")

    # add natural length of the peak
    peaks$real_start <- peaks$start
    peaks$real_end <- peaks$end
    peaks <- peaks %>% filter(chr != "chrM")
    npeaks <- nrow(peaks)
    peaks <- peaks %>% mutate(
        peak_id = paste0("peak_", row_number()),
        peak_length = end - start,
        start = summit - 400,
        end = summit + 400
    )
    
    motifs <- do.call(rbind, lapply(seq(1, npeaks, by = 2000), function(i) {
        peaks1 <- peaks[i:min(i + 1999, npeaks), ]
        options(scipen = 999)
        write.table(peaks1, file = paste0("rep_", ceiling(i / 2000), "_", peakfile), sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
        find_noncanonical_human(peaks = peaks1, peaks_file = paste0("rep_", ceiling(i / 2000), "_", peakfile), motif)
    } ) )

    motifs$flanking <- NA
    # Copy the same ebox_seq to flanking assignments as in mouse
    # If you want to use the same mapping, you can copy-paste the block from mouse here
    # If you want to use a different mapping for human, adjust accordingly

    # --- BEGIN COPY FROM MOUSE ---
    motifs$flanking[motifs$ebox_seq == "ACATATGA"] <- "ACAT-TCAT"
    motifs$flanking[motifs$ebox_seq == "ACATATGT"] <- "ACAT-ACAT"
    motifs$flanking[motifs$ebox_seq == "ACATATGG"] <- "ACAT-CCAT"
    motifs$flanking[motifs$ebox_seq == "ACATATGC"] <- "ACAT-GCAT"
    motifs$flanking[motifs$ebox_seq == "TCATATGA"] <- "TCAT-TCAT"
    motifs$flanking[motifs$ebox_seq == "TCATATGT"] <- "TCAT-ACAT"
    motifs$flanking[motifs$ebox_seq == "TCATATGG"] <- "TCAT-CCAT"
    motifs$flanking[motifs$ebox_seq == "TCATATGC"] <- "TCAT-GCAT"
    motifs$flanking[motifs$ebox_seq == "GCATATGA"] <- "GCAT-TCAT"
    motifs$flanking[motifs$ebox_seq == "GCATATGT"] <- "GCAT-ACAT"
    motifs$flanking[motifs$ebox_seq == "GCATATGG"] <- "GCAT-CCAT"
    motifs$flanking[motifs$ebox_seq == "GCATATGC"] <- "GCAT-GCAT"
    motifs$flanking[motifs$ebox_seq == "CCATATGA"] <- "CCAT-TCAT"
    motifs$flanking[motifs$ebox_seq == "CCATATGT"] <- "CCAT-ACAT"
    motifs$flanking[motifs$ebox_seq == "CCATATGG"] <- "CCAT-CCAT"
    motifs$flanking[motifs$ebox_seq == "CCATATGC"] <- "CCAT-GCAT"
    # ... (copy all the rest of the assignments from the mouse block above)
    # --- END COPY FROM MOUSE ---

    return(motifs)
}

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]

motif <- c(".CA..TG." )
peaks <- read.table(file, header = FALSE)
all_eboxes <- find_noncanonical_human_extended(peaks, file, motif)

estudio<-gsub("_confidence.bed", "", file)
all_eboxes$study<-estudio
saveRDS(all_eboxes,paste0(estudio,"_flanking.RDS") )
' > find_eboxes_human_flanking.R



module load R
module load BEDTools

cd human

for file in *confidence.bed
do
  sbatch --partition=normal --wrap="Rscript ../find_eboxes_human_flanking.R $file"
done


rm rep_*bed
rm rep_*fa

cd ../mouse

for file in *confidence.bed
do
  sbatch --partition=normal --wrap="Rscript ../find_eboxes_mouse_flanking.R $file"
done


rm rep_*bed
rm rep_*fa



