# Scaffolding B. dahlbomii genome with juicer and Hi-C

### Setup juicer
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/
cp -rf /work/gif/00_JuicerSkeleton 01_Juicer

/work/gif/remkv6/Toth/12_Bombus_dahlbomii/01_Juicer/references
#rename scaffolds so they are more managable
awk '/^>/{print ">Scaffold_" ++i; next}{print}' /work/gif2/data/amytoth/Dovetail_data_2021/Bombus_dahlbomii/HiRise/jasmine-iow2045-mb-hirise-5830x__09-30-2021__hic_output.fasta >B_dahlbomiiGenome.fasta
ml bwa;bwa index B_dahlbomiiGenome.fasta
cd ../
bioawk -c fastx '{print $name"\t"length($seq)}' references/B_dahlbomiiGenome.fasta >chrom.sizes
cd fastq/
ln -s /work/gif2/data/amytoth/Dovetail_data_2021/Bombus_dahlbomii/HiRise/DTG-OmniC-150_R1_001.fastq Bdahlbomii_R1.fastq
ln -s /work/gif2/data/amytoth/Dovetail_data_2021/Bombus_dahlbomii/HiRise/DTG-OmniC-150_R2_001.fastq Bdahlbomii_R2.fastq
cd ../splits
split -a 3 -l 90000000 -d --additional-suffix=_R1.fastq ../fastq/Bdahlbomii_R1.fastq &
split -a 3 -l 90000000 -d --additional-suffix=_R2.fastq ../fastq/Bdahlbomii_R2.fastq &


echo "ml jdk; ml bwa; ml juicer/1.5.7; juicer.sh -d /work/gif/remkv6/Toth/12_Bombus_dahlbomii/01_Juicer -p chrom.sizes -s none -z references/B_dahlbomiiGenome.fasta  -q short -Q 2:00:00 -l medium -L 12:00:00 -t 36" >Dahljuicer.sh
```

### Run 3ddna
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/01_Juicer/3Ddna/01_3d-dnaDefault
ln -s ../../references/B_dahlbomiiGenome.fasta
ln -s ../../aligned/merged_nodups.txt


module load miniconda3/4.3.30-qdauveb;source activate 3d-dna;module load jdk;module load parallel;cd /work/gif/remkv6/Toth/12_Bombus_dahlbomii/01_Juicer/3Ddna/01_3d-dnaDefault ;bash run-asm-pipeline.sh B_dahlbomiiGenome.fasta merged_nodups.txt

#Finalize -- unused
 echo "module load miniconda3/4.3.30-qdauveb;source activate 3d-dna;module load jdk;module load parallel;cd /work/gif/remkv6/Toth/12_Bombus_dahlbomii/01_Juicer/3Ddna/04_3d-dnaDiploidRepCovFinalize ;bash run-asm-pipeline-post-review.sh --sort-output -s seal -i 500 -r B_dahlbomiiGenome.0.review.assembly /work/gif/remkv6/Toth/12_Bombus_dahlbomii/01_Juicer/references/B_dahlbomiiGenome.fasta /work/gif/remkv6/Toth/12_Bombus_dahlbomii/01_Juicer/aligned/merged_nodups.txt ">3Ddna.sh

 #Finalize with N=18
 module load miniconda3/4.3.30-qdauveb;source activate 3d-dna;module load jdk;module load parallel;cd /work/gif/remkv6/Toth/12_Bombus_dahlbomii/01_Juicer/3Ddna/05_3ddnaFinalizeN18 ;bash run-asm-pipeline-post-review.sh --sort-output -s seal -i 500 -r B_dahlbomiiGenome.0.N18.assembly.assembly /work/gif/remkv6/Toth/12_Bombus_dahlbomii/01_Juicer/references/B_dahlbomiiGenome.fasta /work/gif/remkv6/Toth/12_Bombus_dahlbomii/01_Juicer/aligned/merged_nodups.txt


```
### Assembly stats
Received assembly
```
---------------- Information for assembly 'B_dahlbomiiGenome.fasta' ----------------


                                         Number of scaffolds         40
                                     Total size of scaffolds  265320182
                                            Longest scaffold   32487547
                                           Shortest scaffold      19293
                                 Number of scaffolds > 1K nt         40 100.0%
                                Number of scaffolds > 10K nt         40 100.0%
                               Number of scaffolds > 100K nt         32  80.0%
                                 Number of scaffolds > 1M nt         24  60.0%
                                Number of scaffolds > 10M nt         13  32.5%
                                          Mean scaffold size    6633005
                                        Median scaffold size    1783693
                                         N50 scaffold length   15367549
                                          L50 scaffold count          7
                                         n90 scaffold length    7505742
                                          L90 scaffold count         15
                                                 scaffold %A      31.67
                                                 scaffold %C      18.30
                                                 scaffold %G      18.31
                                                 scaffold %T      31.72
                                                 scaffold %N       0.00
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs      70.6%
              Percentage of assembly in unscaffolded contigs      29.4%
                      Average number of contigs per scaffold        1.4
Average length of break (>25 Ns) between contigs in scaffold        100

                                           Number of contigs         57
                              Number of contigs in scaffolds         29
                          Number of contigs not in scaffolds         28
                                       Total size of contigs  265318482
                                              Longest contig   19374274
                                             Shortest contig      19293
                                   Number of contigs > 1K nt         57 100.0%
                                  Number of contigs > 10K nt         57 100.0%
                                 Number of contigs > 100K nt         49  86.0%
                                   Number of contigs > 1M nt         38  66.7%
                                  Number of contigs > 10M nt         13  22.8%
                                            Mean contig size    4654710
                                          Median contig size    2035224
                                           N50 contig length   12007047
                                            L50 contig count          9
                                           n90 contig length    3239930
                                            L90 contig count         25
                                                   contig %A      31.67
                                                   contig %C      18.30
                                                   contig %G      18.31
                                                   contig %T      31.72
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0
```

Scaffolded N-18 assembly
```
---------------- Information for assembly 'B_dahlbomiiGenome.FINAL.fasta' ----------------


                                         Number of scaffolds         18
                                     Total size of scaffolds  265332182
                                            Longest scaffold   32487547
                                           Shortest scaffold      33372
                                 Number of scaffolds > 1K nt         18 100.0%
                                Number of scaffolds > 10K nt         18 100.0%
                               Number of scaffolds > 100K nt         17  94.4%
                                 Number of scaffolds > 1M nt         17  94.4%
                                Number of scaffolds > 10M nt         14  77.8%
                                          Mean scaffold size   14740677
                                        Median scaffold size   15367549
                                         N50 scaffold length   17156755
                                          L50 scaffold count          7
                                         n90 scaffold length   10784557
                                          L90 scaffold count         14
                                                 scaffold %A      31.75
                                                 scaffold %C      18.35
                                                 scaffold %G      18.26
                                                 scaffold %T      31.64
                                                 scaffold %N       0.01
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs      92.7%
              Percentage of assembly in unscaffolded contigs       7.3%
                      Average number of contigs per scaffold        3.3
Average length of break (>25 Ns) between contigs in scaffold        334

                                           Number of contigs         59
                              Number of contigs in scaffolds         57
                          Number of contigs not in scaffolds          2
                                       Total size of contigs  265318482
                                              Longest contig   19374274
                                             Shortest contig       5000
                                   Number of contigs > 1K nt         59 100.0%
                                  Number of contigs > 10K nt         58  98.3%
                                 Number of contigs > 100K nt         50  84.7%
                                   Number of contigs > 1M nt         37  62.7%
                                  Number of contigs > 10M nt         13  22.0%
                                            Mean contig size    4496923
                                          Median contig size    1877921
                                           N50 contig length   12007047
                                            L50 contig count          9
                                           n90 contig length    3239930
                                            L90 contig count         25
                                                   contig %A      31.75
                                                   contig %C      18.35
                                                   contig %G      18.26
                                                   contig %T      31.64
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0
```
n=19 Assembly
```
---------------- Information for assembly '../04_RetryN19_Nobreaks/B_dahlbomiiGenome.FINAL.fasta' ----------------


                                         Number of scaffolds         20
                                     Total size of scaffolds  265330182
                                            Longest scaffold   32487547
                                           Shortest scaffold      33372
                                 Number of scaffolds > 1K nt         20 100.0%
                                Number of scaffolds > 10K nt         20 100.0%
                               Number of scaffolds > 100K nt         19  95.0%
                                 Number of scaffolds > 1M nt         18  90.0%
                                Number of scaffolds > 10M nt         14  70.0%
                                          Mean scaffold size   13266509
                                        Median scaffold size   13970555
                                         N50 scaffold length   15794180
                                          L50 scaffold count          7
                                         n90 scaffold length   11734535
                                          L90 scaffold count         14
                                                 scaffold %A      31.75
                                                 scaffold %C      18.35
                                                 scaffold %G      18.26
                                                 scaffold %T      31.64
                                                 scaffold %N       0.00
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs      83.9%
              Percentage of assembly in unscaffolded contigs      16.1%
                      Average number of contigs per scaffold        2.9
Average length of break (>25 Ns) between contigs in scaffold        316

                                           Number of contigs         57
                              Number of contigs in scaffolds         51
                          Number of contigs not in scaffolds          6
                                       Total size of contigs  265318482
                                              Longest contig   19374274
                                             Shortest contig      19293
                                   Number of contigs > 1K nt         57 100.0%
                                  Number of contigs > 10K nt         57 100.0%
                                 Number of contigs > 100K nt         49  86.0%
                                   Number of contigs > 1M nt         38  66.7%
                                  Number of contigs > 10M nt         13  22.8%
                                            Mean contig size    4654710
                                          Median contig size    2035224
                                           N50 contig length   12007047
                                            L50 contig count          9
                                           n90 contig length    3239930
                                            L90 contig count         25
                                                   contig %A      31.75
                                                   contig %C      18.35
                                                   contig %G      18.26
                                                   contig %T      31.64
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0

```
### Run Busco
Received Assembly
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/02_Busco/01_ReceivedAssembly
ln -s ../../01_Juicer/references/B_dahlbomiiGenome.fasta

echo  "ml miniconda3;source activate busco5_env ; busco -i B_dahlbomiiGenome.fasta -o B_dahlbomiiGenome.Busco -m geno --auto-lineage-euk --long --augustus -c 35 -f"  >busco.sh

--------------------------------------------------
|Results from generic domain eukaryota_odb10      |
--------------------------------------------------
|C:89.4%[S:89.4%,D:0.0%],F:1.2%,M:9.4%,n:255      |
|228    Complete BUSCOs (C)                       |
|228    Complete and single-copy BUSCOs (S)       |
|0      Complete and duplicated BUSCOs (D)        |
|3      Fragmented BUSCOs (F)                     |
|24     Missing BUSCOs (M)                        |
|255    Total BUSCO groups searched               |
--------------------------------------------------

--------------------------------------------------
|Results from dataset hymenoptera_odb10           |
--------------------------------------------------
|C:92.0%[S:91.7%,D:0.3%],F:0.3%,M:7.7%,n:5991     |
|5513   Complete BUSCOs (C)                       |
|5495   Complete and single-copy BUSCOs (S)       |
|18     Complete and duplicated BUSCOs (D)        |
|17     Fragmented BUSCOs (F)                     |
|461    Missing BUSCOs (M)                        |
|5991   Total BUSCO groups searched               |
--------------------------------------------------

```
Scaffolded assembly, n=19
```
--------------------------------------------------
|Results from generic domain eukaryota_odb10      |
--------------------------------------------------
|C:89.0%[S:89.0%,D:0.0%],F:1.6%,M:9.4%,n:255      |
|227    Complete BUSCOs (C)                       |
|227    Complete and single-copy BUSCOs (S)       |
|0      Complete and duplicated BUSCOs (D)        |
|4      Fragmented BUSCOs (F)                     |
|24     Missing BUSCOs (M)                        |
|255    Total BUSCO groups searched               |
--------------------------------------------------

--------------------------------------------------
|Results from dataset hymenoptera_odb10           |
--------------------------------------------------
|C:91.9%[S:91.6%,D:0.3%],F:0.2%,M:7.9%,n:5991     |
|5508   Complete BUSCOs (C)                       |
|5488   Complete and single-copy BUSCOs (S)       |
|20     Complete and duplicated BUSCOs (D)        |
|14     Fragmented BUSCOs (F)                     |
|469    Missing BUSCOs (M)                        |
|5991   Total BUSCO groups searched               |
--------------------------------------------------

```
