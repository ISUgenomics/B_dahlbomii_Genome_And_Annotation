# fix contamination found by ncbi and recalculate stats

Sequence name,  length,            span(s),                            apparent source
HiC_scaffold_2     24216648       24097917..24216648       prok:high GC Gram+
HiC_scaffold_6     18407461       9764792..9764872           adaptor:NGB00972.1-not_cleaned



```

# lican used the nocontam genome, so need to remove the genes that were on the mis-scaffolded contaminant on chromosome 2
/work/gif3/masonbrink/Toth/01_dahlbomii/01_FunctionalAnnotation
less SortedFunctionalAnnotationBdahlbomii.gff3 |awk -F"\t" '{if( $1=="HiC_scaffold_2" &&  $5>24097916){} else {print} }'   >BifidoRemovedSortedFunctionalAnnotationBdahlbomii.gff3



cp -rf /work/LAS/amytoth-lab/2023_B_dahlbomii_Genome/FinalFilesOfImportance/NoContaminants/* .

# this one still has the contaminating scaffolds within the genome, they need removed. 
ln -s ../01_FunctionalAnnotation/BifidoRemovedSortedFunctionalAnnotationBdahlbomii.gff3

# starting number of genes
awk '$3=="gene"' /work/LAS/amytoth-lab/2023_B_dahlbomii_Genome/FinalFilesOfImportance/NoContaminants/NoContaminantB_dahlbomiiBraker.gff3 |wc
  15113  136017  941696

# after eliminating the bifidobacterium insertion
grep ">" FinalB_dahlbomiiGenome.fasta |sed 's/>//g' |grep -w -f - BifidoRemovedSortedFunctionalAnnotationBdahlbomii.gff3 >FixedNoContaminantB_dahlbomiiBraker.gff3

less FixedNoContaminantB_dahlbomiiBraker.gff3 |awk '$3=="gene"' |wc
  15028  253815 2445328
```

Fix the genome
```
/work/gif3/masonbrink/Toth/01_dahlbomii/02_FindDahlbomiiContam
bioawk -c fastx '{print $name":1-"length($seq)}' FinalB_dahlbomiiGenome.fasta >ChromosomeCoordinates2Extract.list
samtools faidx FinalB_dahlbomiiGenome.fasta -r ChromosomeCoordinates2Extract.list >BifidoRemovedB_dahlbomiiGenome.fasta

# are the coordinates fixed when fasta lengths are assessed from the new file? yes
bioawk -c fastx '{print $name":1-"length($seq)}' BifidoRemovedB_dahlbomiiGenome.fasta |less

# softmask the supposed illumina adapter


# coordinates to mask
HiC_scaffold_6  9764792 9764872

vi maskCoords.bed

ml bedtools2;bedtools maskfasta -fi BifidoRemovedB_dahlbomiiGenome.fasta -fo IlluminaMaskedBifidoRemovedB_dahlbomiiGenome.fasta -bed maskCoords.bed 

#change to samtools format
vi maskCoords.bed 
# did it get masked? yes
samtools faidx BifidoRemovedB_dahlbomiiGenome.fasta -r maskCoords.bed                        
>HiC_scaffold_6:9764792-9764872
ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGATCTCTCTCTTTTCCTCCT
CCCTCCGTTGTTGTTGTTGAG
samtools faidx IlluminaMaskedBifidoRemovedB_dahlbomiiGenome.fasta -r maskCoords.bed
>HiC_scaffold_6:9764792-9764872
ANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNN

```

Get the proteins and transcripts
```
gffread -g BifidoRemovedB_dahlbomiiGenome.fasta FixedNoContaminantB_dahlbomiiBraker.gff3 -y B_dahlbomiiProteins.fasta -x B_dahlbomiiTranscripts.fasta

```


### Recompute statistics that have changed

gene characteristics
```
 awk '$3=="gene"{print $5-$4}' FunctionalAnnotationB_dahlbomii.gff3|summary.sh
Total:  101,173,916
Count:  15,028
Mean:   6,732
Median: 2,552
Min:    200
Max:    330,148

awk '$3=="mRNA"{print $5-$4}' FunctionalAnnotationB_dahlbomii.gff3|summary.sh
Total:  137,712,600
Count:  19,197
Mean:   7,173
Median: 2,677
Min:    200
Max:    330,148

awk '$3=="CDS"{print $5-$4}' FunctionalAnnotationB_dahlbomii.gff3|summary.sh
Total:  28,166,128
Count:  117,484
Mean:   239
Median: 173
Min:    2
Max:    14,333
```

Genome assembly stats
```
 new_Assemblathon.pl B_dahlbomiiGenome.fasta

---------------- Information for assembly 'B_dahlbomiiGenome.fasta' ----------------


                                         Number of scaffolds         20
                                     Total size of scaffolds  274050013
                                            Longest scaffold   32487547
                                           Shortest scaffold       5000
                                 Number of scaffolds > 1K nt         20 100.0%
                                Number of scaffolds > 10K nt         19  95.0%
                               Number of scaffolds > 100K nt         19  95.0%
                                 Number of scaffolds > 1M nt         18  90.0%
                                Number of scaffolds > 10M nt         15  75.0%
                                          Mean scaffold size   13702501
                                        Median scaffold size   13970555
                                         N50 scaffold length   17299723
                                          L50 scaffold count          7
                                         n90 scaffold length   11609214
                                          L90 scaffold count         15
                                                 scaffold %A      31.69
                                                 scaffold %C      18.29
                                                 scaffold %G      18.34
                                                 scaffold %T      31.64
                                                 scaffold %N       0.04
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs      87.1%
              Percentage of assembly in unscaffolded contigs      12.9%
                      Average number of contigs per scaffold       12.2
Average length of break (>25 Ns) between contigs in scaffold        470

                                           Number of contigs        245
                              Number of contigs in scaffolds        240
                          Number of contigs not in scaffolds          5
                                       Total size of contigs  273944233
                                              Longest contig   19374274
                                             Shortest contig          1
                                   Number of contigs > 1K nt        244  99.6%
                                  Number of contigs > 10K nt        243  99.2%
                                 Number of contigs > 100K nt         63  25.7%
                                   Number of contigs > 1M nt         38  15.5%
                                  Number of contigs > 10M nt         13   5.3%
                                            Mean contig size    1118140
                                          Median contig size      43663
                                           N50 contig length   11383760
                                            L50 contig count         10
                                           n90 contig length    2069253
                                            L90 contig count         28
                                                   contig %A      31.70
                                                   contig %C      18.29
                                                   contig %G      18.35
                                                   contig %T      31.66
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0
```
How many of our mRNAs were functionally annotated by a hit in NCBI NR or SWISS-prot?
```
# How many genes and mRNAs were functionally annotated
awk '$3=="mRNA"' FunctionalAnnotationB_dahlbomii.gff3 |grep -c "Note="
15180

#total mRNAs present
awk '$3=="mRNA"' FunctionalAnnotationB_dahlbomii.gff3 |wc
  19197  340320 3679236

#proportion of mRNAs that have a functional annotation from swissprot or NR
15,180/19,197=79.07%

#number of genes that gained a functional annotation via having an annotated mRNA
less FunctionalAnnotationB_dahlbomii.gff3 |awk '$3=="gene"' |grep -c "Note="
11086

# number of genes
less FunctionalAnnotationB_dahlbomii.gff3 |awk '$3=="gene"' |wc
  15028  253815 2445328

11086/15189= 73.0%


# number of transposon associated genes
 grep -i -e "transposase" -e "transcriptase" -e "retroelement" -e "helitron" -e "reverse" -e "helitron" -e "transposon" FunctionalAnnotationB_dahlbomii.gff3 |awk '$3=="gene"' |wc
     33     592    5861

```

How many alternatively spliced transcripts per gene do we have
```
How many transcripts per gene?
awk '$3=="mRNA"' FunctionalAnnotationB_dahlbomii.gff3 |awk '{print $9}' |sed 's/Parent=/\t/g' |sed 's/;/\t/g' |cut -f 3 |sort |uniq -c |awk '{print $1}' |sort |uniq -c |sort -k2,2n
  
  11929 1
   2328 2
    563 3
    143 4
     48 5
     11 6
      4 7
      1 8
      1 9
```
BUSCO on the final genome
```
Final genome buscos 
hymenoptera_odb10
        C:97.7%[S:97.4%,D:0.3%],F:0.5%,M:1.8%,n:5991
        5856    Complete BUSCOs (C)
        5836    Complete and single-copy BUSCOs (S)
        20      Complete and duplicated BUSCOs (D)
        29      Fragmented BUSCOs (F)
        106     Missing BUSCOs (M)
        5991    Total BUSCO groups searched


# run busco for eukaryota to put in table S1
/work/gif3/masonbrink/Toth/01_dahlbomii/02_FindDahlbomiiContam/FinalFiles
 echo "sh runBuscoEukaryota.sh B_dahlbomiiGenome.fasta" >EukBusco.sh

runBuscoEukaryota.sh
 ########################################################
 #!/bin/bash
#runBusco.sh
#Here is how to run this script: sh runBusco.sh Genome.fasta
Genome="$1"

ml miniconda3;source activate busco5_env ;
busco -i ${Genome} \
-o ${Genome%.*}_Eukaryota_Busco \
-m geno \
-l eukaryota_odb10 \
--long \
--augustus \
-c 35 \
-f
############################################################
```

How many gaps are in the genome?
```
#Counts per scaffold.  Four chromosome scale sequences had zero gaps!, meaning they were completely assembled as contigs
 less B_dahlbomiiGenome.fasta |tr "\n" "\t" |sed 's/>/\n>/g' |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'NR>1' |sed 's/N*N/\n/g' |awk -F'\t' '/^>/ {if (name != "") print name, count; name = substr($0, 2, index($0, "\t") - 1); count = 0; next} {count++} END {if (name != "") print name, count}'
HiC_scaffold_1   1
HiC_scaffold_2   4
HiC_scaffold_3   2
HiC_scaffold_4   0
HiC_scaffold_5   0
HiC_scaffold_6   69
HiC_scaffold_7   1
HiC_scaffold_8   47
HiC_scaffold_9   2
HiC_scaffold_10  1
HiC_scaffold_11  2
HiC_scaffold_12  2
HiC_scaffold_13  0
HiC_scaffold_14  5
HiC_scaffold_15  76
HiC_scaffold_16  2
HiC_scaffold_17  9
HiC_scaffold_18  0
HiC_scaffold_20  0
HiC_scaffold_32  0

#How many gaps total?
less B_dahlbomiiGenome.fasta |tr "\n" "\t" |sed 's/>/\n>/g' |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'NR>1' |sed 's/N*N/\n/g' |awk -F'\t' '/^>/ {if (name != "") print name, count; name = substr($0, 2, index($0, "\t") - 1); count = 0; next} {count++} END {if (name != "") print name, count}' |awk '{print $2}' |summary.sh
Total:  223
Count:  20
Mean:   11
Median: 2
Min:    0
Max:    76

```

Number of N's per 100kbp
```
274050013 Genome length

Number of N's in the genome
less B_dahlbomiiGenome.fasta |grep -v ">" |sed 's/N/N\n/g' |grep -c "N"
105580
105580/1000 =106

247050013/1000 = 247,050

106/247050= 4.29e-4 
4.29e-4 * 100 = 0.0429%

```

