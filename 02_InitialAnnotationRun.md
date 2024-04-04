# Annotate the B. dahlbomii genome

### Fastqc of RNAseq
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/05_Fastqc
for f in /work/gif/archiveNova/2022_TothBeeAssemblies/2023_B_dahlbomii_RNA/*gz; do ln -s $f;done
for f in *gz; do echo "ml fastqc;  fastqc "$f ; done >fastqc.sh

```

### RNA-seq alignment
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/06_RNA_Alignment

for f in /work/gif/archiveNova/2022_TothBeeAssemblies/2023_B_dahlbomii_RNA/*gz; do ln -s $f;done

ln -s ../02_Busco/04_RetryN19_Nobreaks/B_dahlbomiiGenome.FINAL.fasta

echo "ml star;STAR  --runMode genomeGenerate --genomeDir /work/gif/remkv6/Toth/12_Bombus_dahlbomii/08_RNA_Alignment_n19 --genomeFastaFiles B_dahlbomiiGenome.FINAL.fasta" >dahlMakedb.sh

# timed out after 7 days trying alignment with split fastq
~/common_scripts/fastq-splitter.pl --n-parts 10 B-dahl-3_S1_L001_R1_001.fastq
~/common_scripts/fastq-splitter.pl --n-parts 10 B-dahl-3_S1_L001_R2_001.fastq


for f in *R1_001.part*.fastq; do echo "ml star;STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate  --genomeDir /work/gif/remkv6/Toth/12_Bombus_dahlbomii/08_RNA_Alignment_n19 --outFileNamePrefix  "${f}" --readFilesIn  "$f" "${f};done |sed 's/R1/R2/3' >dahlAlignStar.sh

# merge the bams
ls -lrth *out.bam |awk '$5!="0"{print $9}' >bam.list
echo "ml samtools; samtools merge -@ 24  -b bam.list MergedRNA.bam ;samtools sort -o MergedRNA_sorted.bam -T TEMP  --threads 24   MergedRNA.bam; samtools index MergedRNA_sorted.bam ">Dahl_merge.sh
```

### Call repeats in the genome
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/04_RepeatModeler

 ln -s ../02_Busco/04_RetryN19_Nobreaks/B_dahlbomiiGenome.FINAL.fasta

ml miniconda3; source activate repeatmodeler2;
BuildDatabase -name dahlbomii B_dahlbomiiGenome.FINAL.fasta
ml miniconda3; source activate repeatmodeler2; RepeatModeler -database dahlbomii -pa 36
ln -s */consensi.fa.classified
ml miniconda3; source activate repeatmasker; RepeatMasker -pa 36 -norna  -dir RepeatmaskerOut -small -gff -lib consensi.fa.classified B_dahlbomiiGenome.FINAL.fasta
#softmask
ml bedtools2;bedtools maskfasta -fi B_dahlbomiiGenome.FINAL.fasta -fo SoftmaskedB_dahlbomiiGenome.FINAL.fasta -bed RepeatmaskerOut//B_dahlbomiiGenome.FINAL.fasta.out.gff  -soft
```


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
 cp -rf ../../../../02_Busco/04_RetryN19_Nobreaks/B_dahlbomiiGenome.FINAL_Busco/run_hymenoptera_odb10/augustus_output/retraining_parameters/BUSCO_B_dahlbomiiGenome.FINAL_Busco .


ml bamtools;ml genemark-et/4.69-xeer4uv;ml augustus/3.4.0-py39-gvoxpx6; ml braker/3.0.2;braker.pl --genome=SoftmaskedB_dahlbomiiGenome.FINAL.fasta --softmasking --species=BUSCO_B_dahlbomiiGenome.FINAL_Busco --bam=MergedRNA_sorted.bam  --AUGUSTUS_CONFIG_PATH=/work/gif/remkv6/Toth/12_Bombus_dahlbomii/07_braker/Augustus/config/ --overwrite --useexisting




ml braker;ml augustus;cat braker.gtf |gtf2gff.pl --gff3 --out=B_dahlbomiiBraker.gff3
```
### Annotation Stats
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
### Busco on Annotation
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/02_Busco/05_BrakerAnnotation
#remove terminator in sequence
sed 's/\*//g' braker.aa >B_dahlbomiiBraker_Proteins.fasta
ln -s ../../07_braker/braker/B_dahlbomiiBraker_Proteins.fasta

sh runBuscoProteinsToth.sh B_dahlbomiiBraker_Proteins.fasta
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

Results
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
