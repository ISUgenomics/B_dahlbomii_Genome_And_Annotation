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