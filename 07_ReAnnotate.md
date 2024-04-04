# start annotation on B. dahlbomii with new genome
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/13_ReAnnotate

for f in  ../08_RNA_Alignment_n19/*fastq; do ls $f;done
ln -s ../11_Juicer2/3Ddna/04_3d-dnaDiploidRepCovFinalize/ConcatDahlbomiiGenome.FINAL.fasta

echo "ml star;STAR  --runMode genomeGenerate --genomeDir /work/gif/remkv6/Toth/12_Bombus_dahlbomii/13_ReAnnotate --genomeFastaFiles ConcatDahlbomiiGenome.FINAL.fasta" >dahlMakedb.sh


for f in *R1_001.part*.fastq; do echo "ml star;STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate  --genomeDir /work/gif/remkv6/Toth/12_Bombus_dahlbomii/13_ReAnnotate --outFileNamePrefix  "${f}" --readFilesIn  "$f" "${f};done |sed 's/R1/R2/3' >dahlAlignStar.sh

# merge the bams
ls -lrth *out.bam |awk '$5!="0"{print $9}' >bam.list
echo "ml samtools; samtools merge -@ 24  -b bam.list MergedRNA.bam ;samtools sort -o MergedRNA_sorted.bam -T TEMP  --threads 24   MergedRNA.bam; samtools index MergedRNA_sorted.bam ">Dahl_merge.sh
```

### Call repeats in the genome
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/13_ReAnnotate
ln -s ../04_RepeatModeler/consensi.fa.classified

ml miniconda3; source activate repeatmasker; RepeatMasker -pa 36 -norna  -dir RepeatmaskerOut -small -gff -lib consensi.fa.classified ConcatDahlbomiiGenome.FINAL.fasta
#softmask may need rerun the script due to not knowing the identity of the out file.
ml bedtools2;bedtools maskfasta -fi ConcatDahlbomiiGenome.FINAL.fasta -fo SoftmaskedConcatDahlbomiiGenome.FINAL.fasta -bed ConcatDahlbomiiGenome.FINAL.fasta.out.gff  -soft
```
Left off here on 08/07/23 for alignments, busco, and softmasking to run.

### Run Braker
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/07_braker

ln -s ../04_RepeatModeler/SoftmaskedB_dahlbomiiGenome.FINAL.fasta
ln -s ../08_RNA_Alignment_n19/MergedRNA_sorted.bam
ln -s ../08_RNA_Alignment_n19/MergedRNA_sorted.bam.bai

git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augustus;
cp -rf /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/augustus-3.3.2-rtcnsefyulxnfscwgpwi5tc7civqwvsq/bin/ .
cd scripts
cp -rf  /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/braker-2.1.2-75wblifp2zieps5rf7tzp7ajcwvzo2oz/cfg/ .
cd ../config/species/
cp -rf ../../../ConcatDahlbomiiGenome.FINAL_Passerformes_Busco/run_hymenoptera_odb10/augustus_output/retraining_parameters/BUSCO_ConcatDahlbomiiGenome.FINAL_Passerformes_Busco .


ml bamtools;ml genemark-et/4.69-xeer4uv; ml braker/3.0.2;braker.pl --genome=SoftmaskedConcatDahlbomiiGenome.FINAL.fasta  --species=BUSCO_ConcatDahlbomiiGenome.FINAL_Passerformes_Busco --bam=MergedRNA_sorted.bam  --AUGUSTUS_CONFIG_PATH=/work/gif/remkv6/Toth/12_Bombus_dahlbomii/13_ReAnnotate/Augustus/config --overwrite --useexisting




ml braker;ml augustus;cat braker.gtf |gtf2gff.pl --gff3 --out=B_dahlbomiiBraker.gff3
```
### Annotation Stats

Old Stats
```
awk '$3=="gene"{print $5-$4}' B_dahlbomiiBraker.gff3|summary.sh
Total:  92,266,192
Count:  16,534
Mean:   5,580
Median: 2,044
Min:    200
Max:    299,836

awk '$3=="mRNA"{print $5-$4}' B_dahlbomiiBraker.gff3|summary.sh
Total:  124,794,620
Count:  20,407
Mean:   6,115
Median: 2,272
Min:    200
Max:    299,836

awk '$3=="CDS"{print $5-$4}' B_dahlbomiiBraker.gff3|summary.sh
Total:  27,581,393
Count:  112,990
Mean:   244
Median: 176
Min:    2
Max:    14,333
```
New Stats
```
 awk '$3=="gene"{print $5-$4}' B_dahlbomiiBraker.gff3|summary.sh
Total:  102,114,211
Count:  15,767
Mean:   6,476
Median: 2,394
Min:    142
Max:    330,148
awk '$3=="mRNA"{print $5-$4}' B_dahlbomiiBraker.gff3|summary.sh
Total:  138,653,209
Count:  19,937
Mean:   6,954
Median: 2,543
Min:    142
Max:    330,148
awk '$3=="CDS"{print $5-$4}' B_dahlbomiiBraker.gff3|summary.sh
Total:  29,061,897
Count:  118,519
Mean:   245
Median: 173
Min:    2
Max:    14,333
```


### Busco on Annotation
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/02_Busco/05_BrakerAnnotation
gffread B_dahlbomiiBraker.gff3 -g ../ConcatDahlbomiiGenome.FINAL.fasta -y Dahlbomii_Proteins.fasta -x Dahlbomii_Transcripts


sh runBuscoProteinsToth.sh Dahlbomii_Proteins.fasta
#runBusco
#############################################################
#!/bin/bash
#runBusco.sh
#Here is how to run this script: sh runBusco.sh Genome.fasta
Genome="$1"

ml purge;ml miniconda3;source activate busco5_env ;
busco -i ${Genome} \
-o ${Genome%.*}_Busco_eukaryota \
-m prot \
-l eukaryota_odb10 \
-c 35 \
-f

busco -i ${Genome} \
-o ${Genome%.*}_Busco_metazoa \
-m prot \
-l metazoa_odb10 \
-c 35 \
-f

busco -i ${Genome} \
-o ${Genome%.*}_Busco_hymenoptera \
-m prot \
-l hymenoptera_odb10 \
-c 35 \
-f
##############################################################
```
Old scores
```
        --------------------------------------------------
        |Results from dataset eukaryota_odb10             |
        --------------------------------------------------
        |C:88.3%[S:65.9%,D:22.4%],F:3.1%,M:8.6%,n:255     |
        |225    Complete BUSCOs (C)                       |
        |168    Complete and single-copy BUSCOs (S)       |
        |57     Complete and duplicated BUSCOs (D)        |
        |8      Fragmented BUSCOs (F)                     |
        |22     Missing BUSCOs (M)                        |
        |255    Total BUSCO groups searched               |
        --------------------------------------------------

        --------------------------------------------------
        |Results from dataset metazoa_odb10               |
        --------------------------------------------------
        |C:89.7%[S:65.7%,D:24.0%],F:1.4%,M:8.9%,n:954     |
        |856    Complete BUSCOs (C)                       |
        |627    Complete and single-copy BUSCOs (S)       |
        |229    Complete and duplicated BUSCOs (D)        |
        |13     Fragmented BUSCOs (F)                     |
        |85     Missing BUSCOs (M)                        |
        |954    Total BUSCO groups searched               |
        --------------------------------------------------

        --------------------------------------------------
        |Results from dataset hymenoptera_odb10           |
        --------------------------------------------------
        |C:89.3%[S:63.5%,D:25.8%],F:1.5%,M:9.2%,n:5991    |
        |5351   Complete BUSCOs (C)                       |
        |3804   Complete and single-copy BUSCOs (S)       |
        |1547   Complete and duplicated BUSCOs (D)        |
        |88     Fragmented BUSCOs (F)                     |
        |552    Missing BUSCOs (M)                        |
        |5991   Total BUSCO groups searched               |
        --------------------------------------------------
```
New Scores!
```
        --------------------------------------------------
        |Results from dataset eukaryota_odb10             |
        --------------------------------------------------
        |C:95.7%[S:72.2%,D:23.5%],F:2.4%,M:1.9%,n:255     |
        |244    Complete BUSCOs (C)                       |
        |184    Complete and single-copy BUSCOs (S)       |
        |60     Complete and duplicated BUSCOs (D)        |
        |6      Fragmented BUSCOs (F)                     |
        |5      Missing BUSCOs (M)                        |
        |255    Total BUSCO groups searched               |
        --------------------------------------------------
        --------------------------------------------------
        |Results from dataset metazoa_odb10               |
        --------------------------------------------------
        |C:96.8%[S:71.9%,D:24.9%],F:1.3%,M:1.9%,n:954     |
        |924    Complete BUSCOs (C)                       |
        |686    Complete and single-copy BUSCOs (S)       |
        |238    Complete and duplicated BUSCOs (D)        |
        |12     Fragmented BUSCOs (F)                     |
        |18     Missing BUSCOs (M)                        |
        |954    Total BUSCO groups searched               |
        --------------------------------------------------
        --------------------------------------------------
        |Results from dataset hymenoptera_odb10           |
        --------------------------------------------------
        |C:95.1%[S:67.8%,D:27.3%],F:1.4%,M:3.5%,n:5991    |
        |5697   Complete BUSCOs (C)                       |
        |4062   Complete and single-copy BUSCOs (S)       |
        |1635   Complete and duplicated BUSCOs (D)        |
        |84     Fragmented BUSCOs (F)                     |
        |210    Missing BUSCOs (M)                        |
        |5991   Total BUSCO groups searched               |
        --------------------------------------------------

```





