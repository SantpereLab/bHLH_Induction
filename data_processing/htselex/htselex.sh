cd p
mkdir htselex
cd htselex



# From a table I downloaded from the ENA, extract the accessions for the bHLH factors
# First, extract all bHLH from the Lambert table
cat lambert_table.txt | awk '{if ($4=="bHLH") print $3}' > human_bHLH.txt

################## YIN 2017
###############################
#############################

cd n/htselex
mkdir yin
cd yin

# Remove those that are bisulfite sequencing
rm *Bis*

# And also remove the "Nor" files, because they are bisulfite.
rm *Nor*

# I'm not using the bisulfite ones because I don't understand them and they are very complex.
# It's better to use the methyl-HT-SELEX files because they are simpler, and as stated in a supplementary [file/table], they work the same for bHLH factors.
# Information on which one is methyl-HT-SELEX and which is not can be found on the website.

grep -f ../human_bHLH.txt filereport_read_run_PRJEB3289_tsv.txt | cut -f8 > yin_ftps_bhlh.txt
while read line; do sbatch --wrap=" wget $line ";done<yin_ftps_bhlh.txt
sbatch --wrap="gunzip *gz"


# Split this operation into the different IDs
ls *fastq | tr "_" "\t" | tr "." "\t" | awk '{print $1"*"$2"*"$5}' | sort | uniq > ids_regexpr.txt
cat ids_regexpr.txt | egrep 'NEUROD1|NEUROD2|NEUROG2|TWIST2|TWIST1|OLIG2|MSGN1|MESP1|MYOD1|MYF5|ASCL1|ASCL2|TCF4|MYC|MYCN|MITF|HEY1|HEY2|TFE3' > ids_regexpr_subset.txt

cat ids_regexpr.txt | egrep 'MAX' > ids_regexpr_subset_max.txt

cat yin_ftps_bhlh.txt | tr "/" "\t" | tr "_" "\t" | cut -f6 | sort | uniq > bhlh_in_yin.txt

while read line; do mkdir $line;done<ids_regexpr.txt

while read line
do
  echo $line
  for fastq in ${line}*fastq
  do
    cd $line
    cat ../$fastq | wc -l | awk '{print $1/4 }' >> ${fastq}_wc.txt
    for motif in CAAATG CAACTG CAAGTG CAATTG CACATG CACCTG CACGTG CACTTG CAGATG CAGCTG CAGGTG CAGTTG CATATG CATCTG CATGTG CATTTG
    do
      grep $motif ../$fastq | wc -l >> ${fastq}_motifs.txt
    done
    cd ..
  done
done<ids_regexpr_subset.txt



# Make the combined file

while read line
do
  echo $line
  cd $line
  for motif in CAAATG CAACTG CAAGTG CAATTG CACATG CACCTG CACGTG CACTTG CAGATG CAGCTG CAGGTG CAGTTG CATATG CATCTG CATGTG CATTTG
  do
    echo $motif >> hexanucleotides.txt
  done
  cat <(echo "motif 1 2 3 4" | tr " " "\t") <(paste <(echo "total_sequence") *_1_*wc.txt *_2_*wc.txt *_3_*wc.txt *_4_*wc.txt) <(paste hexanucleotides.txt *_1_*.fastq_motifs.txt *_2_*.fastq_motifs.txt *_3_*.fastq_motifs.txt *_4_*.fastq_motifs.txt) > $line.motif_counts.txt
  cd ..
done<ids_regexpr_subset.txt


