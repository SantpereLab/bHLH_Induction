cd p
mkdir htselex
cd htselex

#de una tabla que me he descargado del ENA, sacar los accessions de los bHLH factors
#primero sacar todos los bHLH de la tabla del lambert
cat lambert_table.txt | awk '{if ($4=="bHLH") print $3}' > human_bHLH.txt






################## YIN 2017
###############################
#############################

cd n/htselex
mkdir yin
cd yin
#removear las que son bisulphite sequencing
rm *Bis*
#y los Nor tambien removearlos, porque son bis.
rm *Nor*

# Los bisulphite no los uso pq no los entiendo y son muy complejos.
# Mejor usar los metil-htselex, pq son mas simples, y como pone en una suplementaria, funcionan igual para los bHLH.
# Cual es metil-htselex y cual no, la info esta en la pagina web.


grep -f ../human_bHLH.txt filereport_read_run_PRJEB3289_tsv.txt | cut -f8 > yin_ftps_bhlh.txt
while read line; do sbatch --wrap=" wget $line ";done<yin_ftps_bhlh.txt
sbatch --wrap="gunzip *gz"


#dividir esta operacion en los distintos ids 
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



#hacer la file conjunta

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


