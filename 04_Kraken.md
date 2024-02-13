# Identify all contaminants present in genome and reads

### Create Kraken2 database

```
#this section needs to be updated when I have removed the contaminants from the B. dahlbomii genome.
ml miniconda3;source activate kraken2
kraken2-build --download-taxonomy --db BacteriaVirusArchaea/
kraken2-build --download-library viral --db BacteriaVirusArchaea/
kraken2-build --download-library bacteria --db BacteriaVirusArchaea/
kraken2-build --download-library archaea --db BacteriaVirusArchaea/
kraken2-build --download-library fungi --db BacteriaVirusArchaea/
kraken2-build --download-library protozoa --db BacteriaVirusArchaea/
kraken2-build --download-library nematoda --db BacteriaVirusArchaea/
#bioawk -c fastx '{print ">"$name"|kraken:taxid|51029\n"$seq}'  ../01_AlignmentNDeseq/SCNgenome.fasta >HglycinesTax.fa
kraken2-build --add-to-library HglycinesTax.fa --db BacteriaVirusArchaea/
kraken2-build --build --db BacteriaVirusArchaea/ --threads 36

for f in ../01_AlignmentNDeseq/*fastq.gz; do ln -s $f;done
```

### RNASEQ
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/12_Kraken
ln -s /work/gif/remkv6/Baum/12_Endoreduplication/02_Kraken/BacteriaVirusArchaea

paste <(ls -1 *R1*fastq) <(ls -1 *R2*fastq) |tr "\t" " " |while read line; do echo "kraken2 -db /work/gif/remkv6/Baum/12_Endoreduplication/02_Kraken/BacteriaVirusArchaea --threads 36 --report "$line".report   --unclassified-out "${line%.*}"unclassified#.fq --classified-out "${line%.*}"classified#.fq --paired "$line" > "${line%.*}"Kraken.out" ;done  |awk '{print $1,$2,$3,$4,$5,$6,$8,$9,$10,$12,$13,$15,$16,$17,$18,$19,$21}' >kraken.sh

echo "ml miniconda3;source activate kraken2;kraken2 -db /work/gif/remkv6/Baum/12_Endoreduplication/02_Kraken/BacteriaVirusArchaea --threads 36 --report B-dahl-3_S1_L001_R2_001.fastq.
report --unclassified-out B-dahl-3_S1_L001_R2_001unclassified#.fq --classified-out B-dahl-3_S1_L001_R2_001classified#.fq --paired B-dahl-3_S1_L001_R1_001.fastq B-dahl-3_S1_L001_R2_001.fastq > B-dahl-3_S1_L001_R2_001Kraken.out" >kraken.sh

# Summarize the report
awk '$1>0 && $3>100' B-dahl-3_S1_L001_R2_001.fastq.report |uniq|sort -k1,1nr  >B-dahl-3_S1_L001_R2_001.fastq.report.summary

```

### PACBIO HIFI
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/12_Kraken
ln -s ../10_FindBuscos/30_537121167.hifi_reads.fastq.gz

ls -1 *hifi_reads.fastq.gz |while read line; do echo "kraken2 -db /work/gif/remkv6/Baum/12_Endoreduplication/02_Kraken/BacteriaVirusArchaea --threads 36 --report "$line".report --gzip-compressed  --unclassified-out "${line%.*}"unclassified#.fq --classified-out "${line%.*}"classified#.fq --paired "$line" > "${line%.*}"Kraken.out" ;done  


kraken2 -db /work/gif/remkv6/Baum/12_Endoreduplication/02_Kraken/BacteriaVirusArchaea --threads 36 --report 30_537121167.hifi_reads.fastq.gz.report --gzip-compressed  --unclassified-out 30_537121167.hifi_reads.fastqunclassified.fq --classified-out 30_537121167.hifi_reads.fastqclassified.fq 30_537121167.hifi_reads.fastq.gz > 30_537121167.hifi_reads.fastqKraken.out

awk '$1>0 && $3>100' BdahlbomiiHiC_R2.fastq.report |uniq|sort -k1,1nr  >BdahlbomiiHiC_R2.fastq.report.summary

```

### HIC
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/12_Kraken
ln -s /work/gif3/data/amytoth/Dovetail_data_2021/Bombus_dahlbomii/HiRise/DTG-OmniC-150_R1_001.fastq BdahlbomiiHiC_R1.fastq
ln -s /work/gif3/data/amytoth/Dovetail_data_2021/Bombus_dahlbomii/HiRise/DTG-OmniC-150_R2_001.fastq BdahlbomiiHiC_R2.fastq

paste <(ls -1 *HiC_R1*fastq) <(ls -1 *HiC_R2*fastq) |tr "\t" " " |while read line; do echo "kraken2 -db /work/gif/remkv6/Baum/12_Endoreduplication/02_Kraken/BacteriaVirusArchaea --threads 36 --report "$line".report   --unclassified-out "${line%.*}"unclassified#.fq --classified-out "${line%.*}"classified#.fq --paired "$line" > "${line%.*}"Kraken.out" ;done  

echo " ml miniconda3;source activate kraken2;kraken2 -db /work/gif/remkv6/Baum/12_Endoreduplication/02_Kraken/BacteriaVirusArchaea --threads 36 --report  BdahlbomiiHiC_R2.fastq.report   --unclassified-out BdahlbomiiHiC_R2unclassified#.fq --classified-out BdahlbomiiHiC_R2classified#.fq --paired BdahlbomiiHiC_R1.fastq BdahlbomiiHiC_R2.fastq > BdahlbomiiHiC_R2Kraken.out" >HiCKraken.sh

awk '$1>0 && $3>100' BdahlbomiiHiC_R2.fastq.report |uniq|sort -k1,1nr  >BdahlbomiiHiC_R2.fastq.report.summary
```

These all came out a little odd because of the use of an older database that included the H. glycines genome.  