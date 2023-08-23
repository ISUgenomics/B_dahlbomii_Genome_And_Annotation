# Identify contaminants in the genome assembly of B. dahlbomii



### Are these little contigs just rare assembled haplotypes? Check for low input read mapping rates and use for blobtools

Map long reads
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/14_Blobtools

ln -s ../10_FindBuscos/30_537121167.hifi_reads.fastq.gz
ln -s ../13_ReAnnotate/ConcatDahlbomiiGenome.FINAL.fasta


echo "sh runMinimapNbamSort.sh 30_537121167.hifi_reads.fastq.gz ConcatDahlbomiiGenome.FINAL.fasta" >ReadMapping.sh

#runMinimap.sh
##############################################################################
#!/bin/bash
query=$1
target=$2
outname="${query%.*}_${target%.*}_minimap2.sam"
module load minimap2
minimap2 -x asm5 -a -t 36 $target $query > ${outname}

ml samtools;samtools view --threads 36 -b -o ${outname%.*}.bam ${outname}
samtools sort  -o ${outname%.*}_sorted.bam -T TEMP --threads 36 ${outname%.*}.bam
samtools index ${outname%.*}_sorted.bam
##############################################################################

```

Megablast to NCBI NT
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/14_Blobtools

fasta-splitter.pl --n-parts 20 ConcatDahlbomiiGenome.FINAL.fasta

for f in *part*fasta; do echo "sh runMegablast2nt.sh "$f;done >blasts.sh

#!/bin/bash
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
#tar -zxvf taxdb.tar.gz

#module load blast-plus
FASTA="$1"
blastn \
-task megablast \n
-query ${FASTA} \
-db  /work/LAS/BioDatabase/BLASTdb/NCBI/Archives/Current/nt \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-culling_limit 10 \
-num_threads 36 \
-evalue 1e-3 \
-out ${FASTA%.**}.vs.nt.cul10.1e3.megablast.out
```



Blobtools
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/14_Blobtools

#9 of the pseudomolecules were taking a long time. Since we know that those are Dahlbomii, I will just substitute the blast hits to bombus from scaffold_5
cat <(head -n 9 AllBlasts.tab) AllBlasts.tab >FixedAllBlasts.tab


module load singularity;module load blobtools2;
singularity shell /opt/rit/singularity/images/blobtools2/2.2.0/blobtools2.simg

cp -rf /work/gif3/masonbrink/Serb/08_ContigElimination/01_Blobtools/taxdump/ .
ln -s ../13_ReAnnotate/ConcatDahlbomiiGenome.FINAL_Passerformes_Busco/run_hymenoptera_odb10/full_table.tsv

blobtools create --fasta ConcatDahlbomiiGenome.FINAL.fasta --cov 30_537121167.hifi_reads.fastq_ConcatDahlbomiiGenome.FINAL_minimap2_sorted.bam --busco full_table.tsv --hits FixedAllBlasts.tab  --taxdump taxdump  test


#connect to novaDTN through my pc's terminal
/work/gif/remkv6/Olsen/Bison/08_Blobtools/03_Blobenate
module load singularity;module load blobtools2;
singularity shell /opt/rit/singularity/images/blobtools2/2.2.0/blobtools2.simg

blobtools view --interactive test


Elimination of contigs that had zero coverage and were not one of the following phyla: Mollusca, Brachiopoda, Hemichordata, No-hit, or undef.   
3169 contigs total evaluated, 2,923 pass these filters
Total contaminated contigs 3169-2923= 246 contigs


#Took default circle picture of PGA assembly 2, mapped nanopore reads, blastn to NCBI nt
```
![Blobtools Blob plot](./Assets/Cales.blob.circle.png)
![Blobtools Snail plot](./Assets/Cales.snail.png)




### Create lists of contigs/scaffolds that meet the qualifications above
```
less CalesBlob.csv |sed 's/,/\t/g' |sed 's/"//g' |sed 's/^\t//g' |awk '$4<10{print $6}' >LowCoverageScaffolds.list
less CalesBlob.csv |sed 's/,/\t/g' |sed 's/"//g' |sed 's/^\t//g' |awk '$5!="Mollusca" ' |awk  '$5!="Echinodermata"' |awk '$5!="Brachiopoda"' |awk '$5!="Hemichordata"' | awk '$5!="undef" ' |awk '$5!="no-hit"' |awk 'NR>1{print $6}' > ContaminatingScaffolds.list


```

### Create lists of contigs that are just repeats that were not incorporated
```
/work/gif/remkv6/Serb/02_ContigElimination/03_RepeatOnlyContigs
cp /work/LAS/serb-lab/2021_Serb_Pterimorphia/03_Cales_genome/ContainedRemovedCalesGenome.fasta.out.EDTA.gff.gz .
cp /work/LAS/serb-lab/2021_Serb_Pterimorphia/03_Cales_genome/ContainedRemovedCalesGenome.fasta.out.RMod.gff.gz .

for f in *; do gunzip $f;done

#use bedtools to merge the repeat coordinates
cat ContainedRemovedCalesGenome.fasta.out.EDTA.gff ContainedRemovedCalesGenome.fasta.out.RMod.gff |sort -k1,1V -k4,5n |bedtools merge -i - -d 100 >MergedRepeats.bed

#get scaffold lengths (note this is the contam removed genome)
ml bioawk; bioawk -c fastx '{print $name,length($seq)}' ContainedRemovedCalesGenome.fasta >Scaffold.lengths

# get cumulative total of repeats on each scaffold
sed 's/_/\t/2' MergedRepeats.bed |awk '$2>18' |sed 's/\t/_/1' |awk '{print $1,$3-$2}' |awk -v A=0 '{if($1==A) {b=b+$2} else {print A,b; A=$1;b=$2}}' >RepeatTotal.tab

# get the scaffold lengths in the same order as the RepeatTotal.tab file
awk '{print $1}'  RepeatTotal.tab |while read line; do grep -w $line Scaffold.lengths ;done >RepScaffoldsLengthsInOrder.tab

#How many scaffolds are 90% repeat or greater
paste RepScaffoldsLengthsInOrder.tab <(awk 'NR>1' RepeatTotal.tab) |awk '{print $1,$4/$2}' |awk '$2>.89' |wc
  1201    2402   29981

# I actually dont use this, as I merge repeat coordinates with the little contig to big contig mapping coordinates
paste RepScaffoldsLengthsInOrder.tab <(awk 'NR>1' RepeatTotal.tab) |awk '{print $1,$4/$2}' |awk '$2>.8999 {print $1}' >RepeatScaffolds.list
```

### How many contigs have a combined total of containment to pseudomolecules and repeats of greater than 90%
```

# merge the mapping files with the repeat files
less LittleContigs_BigContigs_minimap2.paf |cut -f 1,3,4 |cat - ../03_RepeatOnlyContigs/MergedRepeats.bed |sed 's/_/\t/2'  |awk '$2>18' |sed 's/\t/_/1' |tr " " "\t" |sort -k1,1V -k2,3n |bedtools merge -i - -d 100 >ContainedRepetitiveContigsTotal.bed


#this is adding up the cumulative total (repeats + containment)for each scaffold
sed 's/_/\t/2' ContainedRepetitiveContigsTotal.bed |awk '$2>18' |sed 's/\t/_/1' |awk '{print $1,$3-$2}' |awk -v A=0 '{if($1==A) {b=b+$2} else {print A,b; A=$1;b=$2}}' |awk 'NR>1' >ContainedAndRepeatsTotal.tab

# get scaffold lengths
bioawk -c fastx '{print $name,length($seq)}' ../C_ales_Genome.FINAL.fasta >Scaffold.lengths

# get the scaffold lengths in the same order as the containedAndRepeatsTotal.tab file
awk '{print $1}'  ContainedAndRepeatsTotal.tab |while read line; do grep -w $line Scaffold.lengths ;done >ContainedRepetitiveContigsTotalLengthsInOrder.tab

#how many contigs are contained and/or repetitive at 90% or greater
paste ContainedRepetitiveContigsTotalLengthsInOrder.tab <(awk 'NR>1' ContainedAndRepeatsTotal.tab) |awk '{print $1,$4/$2}' |awk '$2>.8999' |wc
  3528    7056   91306

paste ContainedRepetitiveContigsTotalLengthsInOrder.tab <(awk 'NR>1' ContainedAndRepeatsTotal.tab) |awk '{print $1,$4/$2}' |awk '$2>.8999 {print $1}' >ContainedRepetitiveContigs2Remove.list

#combine contaminating contigs, low coverage contigs, and 90% contained/repeat contigs into a list
cat ContainedRepetitiveContigs2Remove.list ../02_Blobtools/ContaminatingScaffolds.list  ../02_Blobtools/LowCoverageScaffolds.list |sort|uniq>CombinedAllContigs2Remove.list

#How many contigs will be removed using this list
wc CombinedAllContigs2Remove.list
3718  3718 65985 CombinedAllContigs2Remove.list

#Create the end fasta file for subseqeuent use
cat CombinedAllContigs2Remove.list <(awk '{print $1}' Scaffold.lengths )|sort|uniq -c |awk '$1==1{print $2}' |cdbyank ../03_RepeatOnlyContigs/ContainedRemovedCalesGenome.fasta.cidx >CalesGenomeCleanScaffolds.fasta

```

### Results
```
~/common_scripts/new_Assemblathon.pl CalesGenomeCleanScaffolds.fasta

---------------- Information for assembly 'CalesGenomeCleanScaffolds.fasta' ----------------


                                         Number of scaffolds       1049
                                     Total size of scaffolds 2237305496
                                            Longest scaffold  240195547
                                           Shortest scaffold       1000
                                 Number of scaffolds > 1K nt       1047  99.8%
                                Number of scaffolds > 10K nt        536  51.1%
                               Number of scaffolds > 100K nt         24   2.3%
                                 Number of scaffolds > 1M nt         18   1.7%
                                Number of scaffolds > 10M nt         18   1.7%
                                          Mean scaffold size    2132798
                                        Median scaffold size      10655
                                         N50 scaffold length  121434358
                                          L50 scaffold count          7
                                         n90 scaffold length   84817633
                                          L90 scaffold count         16
                                                 scaffold %A      32.08
                                                 scaffold %C      17.87
                                                 scaffold %G      17.87
                                                 scaffold %T      32.09
                                                 scaffold %N       0.09
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs      99.4%
              Percentage of assembly in unscaffolded contigs       0.6%
                      Average number of contigs per scaffold        7.8
Average length of break (>25 Ns) between contigs in scaffold        281

                                           Number of contigs       8191
                              Number of contigs in scaffolds       7380
                          Number of contigs not in scaffolds        811
                                       Total size of contigs 2235298496
                                              Longest contig    7630288
                                             Shortest contig         28
                                   Number of contigs > 1K nt       8108  99.0%
                                  Number of contigs > 10K nt       7046  86.0%
                                 Number of contigs > 100K nt       2905  35.5%
                                   Number of contigs > 1M nt        651   7.9%
                                  Number of contigs > 10M nt          0   0.0%
                                            Mean contig size     272897
                                          Median contig size      50855
                                           N50 contig length    1207478
                                            L50 contig count        508
                                           n90 contig length     133914
                                            L90 contig count       2476
                                                   contig %A      32.11
                                                   contig %C      17.89
                                                   contig %G      17.88
                                                   contig %T      32.12
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0
```

![Blobtools Snail Plot](https://github.com/ISUgenomics/B_dahlbomii_Genome_And_Annotation/blob/main/Notebook_Masonbrink/Assets/BdahlbomiiSnail.png)
