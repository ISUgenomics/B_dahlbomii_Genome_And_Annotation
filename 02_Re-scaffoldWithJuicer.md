# Concatenate missing parts of the genome to existing genome and rescaffold with juicer/3ddna
 ```
 p -rf /work/gif/00_JuicerSkeleton/ 11_Juicer2
cd 11_Juicer2/references/
cat ../../10_FindBuscos/02_Flye/assembly.fasta ../../07_braker/SoftmaskedB_dahlbomiiGenome.FINAL.fasta >ConcatDahlbomiiGenome.fasta
bwa index ConcatDahlbomiiGenome.fasta
cd ../
bioawk -c fastx '{print $name"\t"length($seq)}' references/ConcatDahlbomiiGenome.fasta >chrom.sizes
cd fastq/
ln -s ../../01_Juicer/fastq/Bdahlbomii_R1.fastq
ln -s ../../01_Juicer/fastq/Bdahlbomii_R2.fastq
cd ../splits/
for f in ../../01_Juicer/splits/*fastq; do ln -s $f;done
 ml bwa; ml juicer/1.5.7; juicer.sh -d /work/gif/remkv6/Toth/12_Bombus_dahlbomii/11_Juicer2 -p chrom.sizes -s none -z references/ConcatDahlbomiiGenome.fasta  -q short -Q 2:00:00 -l medium -L 12:00:00 -t 8
```
### 3ddna
```
echo "module load miniconda3/4.3.30-qdauveb;source activate 3d-dna;module load jdk;module load parallel;cd /work/gif/remkv6/Toth/12_Bombus_dahlbomii/11_Juicer2/3Ddna/01_3d-dnaDefault; bash run-asm-pipeline.sh ConcatDahlbomiiGenome.fasta merged_nodups.txt" >3ddna.sh

#Most of these new contigs clustered together into 3 spots, not sure why they were not assembled in the intial assembly. 

echo "module load miniconda3/4.3.30-qdauveb;source activate 3d-dna;module load jdk;module load parallel;cd  /work/gif/remkv6/Toth/12_Bombus_dahlbomii/11_Juicer2/3Ddna/04_3d-dnaDiploidRepCovFinalize; bash run-asm-pipeline-post-review.sh --sort-output -i 500 -r ConcatDahlbomiiGenome.0.review.assembly ConcatDahlbomiiGenome.fasta merged_nodups.txt" >3ddna.sh
```

### Evaluate the new assembly
```
Eukaryota_odb10
S:98.82%, 252
D:0.00%, 0
F:0.00%, 0
I:0.00%, 0
M:1.18%, 3
N:255

Hymenoptera_odb10
S:98.23%, 5885
D:0.28%, 17
F:0.80%, 48
I:0.00%, 0
M:0.68%, 41
N:5991

These are significantly better than the original busco estimates

                                         Number of scaffolds         32
                                     Total size of scaffolds  275147538
                                            Longest scaffold   32487547
                                           Shortest scaffold       5000
                                 Number of scaffolds > 1K nt         32 100.0%
                                Number of scaffolds > 10K nt         30  93.8%
                               Number of scaffolds > 100K nt         23  71.9%
                                 Number of scaffolds > 1M nt         18  56.2%
                                Number of scaffolds > 10M nt         15  46.9%
                                          Mean scaffold size    8598361
                                        Median scaffold size    7505742
                                         N50 scaffold length   17299723
                                          L50 scaffold count          7
                                         n90 scaffold length   11609214
                                          L90 scaffold count         15
                                                 scaffold %A      31.66
                                                 scaffold %C      18.32
                                                 scaffold %G      18.38
                                                 scaffold %T      31.61
                                                 scaffold %N       0.04
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs      86.8%
              Percentage of assembly in unscaffolded contigs      13.2%
                      Average number of contigs per scaffold        8.1
Average length of break (>25 Ns) between contigs in scaffold        469

                                           Number of contigs        258
                              Number of contigs in scaffolds        241
                          Number of contigs not in scaffolds         17
                                       Total size of contigs  275041338
                                              Longest contig   19374274
                                             Shortest contig       5000
                                   Number of contigs > 1K nt        258 100.0%
                                  Number of contigs > 10K nt        256  99.2%
                                 Number of contigs > 100K nt         67  26.0%
                                   Number of contigs > 1M nt         38  14.7%
                                  Number of contigs > 10M nt         13   5.0%
                                            Mean contig size    1066052
                                          Median contig size      45903
                                           N50 contig length   11383760
                                            L50 contig count         10
                                           n90 contig length    1877921
                                            L90 contig count         29
                                                   contig %A      31.67
                                                   contig %C      18.33
                                                   contig %G      18.38
                                                   contig %T      31.62
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0
```
Most of the extra extraneous contigs are likley contamination, as they did not have any signal in the juicbox plot.   

HiC plot from from round 1 + new scaffolds added (tiny bottom right)
![Before Modification](https://github.com/ISUgenomics/B_dahlbomii_Genome_And_Annotation/blob/main/Assets/Round2Before.png)
HiC plot after editing 
![After Modification](https://github.com/ISUgenomics/B_dahlbomii_Genome_And_Annotation/blob/main/Assets/Round2After.png)



