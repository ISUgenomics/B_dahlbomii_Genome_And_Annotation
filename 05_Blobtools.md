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


```
![Blobtools Blob plot](https://github.com/ISUgenomics/B_dahlbomii_Genome_And_Annotation/blob/main/Assets/PostAssemblyBlobplot.png)
![Blobtools Snail plot](https://github.com/ISUgenomics/B_dahlbomii_Genome_And_Annotation/blob/main/Assets/PostAssemblySnail.png)






Final plot after filtering contaminating scaffolds
![Blobtools Snail Plot](https://github.com/ISUgenomics/B_dahlbomii_Genome_And_Annotation/blob/main/Assets/BdahlbomiiSnail.png)
