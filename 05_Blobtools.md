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
Megablast ouput of top best hit
```
HiC_scaffold_1  30195   8.396e+05       HiC_scaffold_5  OU342930.1      77.816  1464352 254517  44864   1576486 2979461 12199451        13654846        0.0     Bombus terrestris       Eukaryota       Bombus terrestris genome assembly, chromosome: 10
HiC_scaffold_10 30195   8.337e+05       HiC_scaffold_5  OU342930.1      77.755  1465226 254105  45171   1487581 2890829 12199451        13654816        0.0     Bombus terrestris       Eukaryota       Bombus terrestris genome assembly, chromosome: 10
HiC_scaffold_11 65598   1.236e+05       HiC_scaffold_11 HG995280.1      96.275  75837   1995    297     7573037 7648548 7502796 7578127 0.0     Bombus pascuorum        Eukaryota       Bombus pascuorum genome assembly, chromosome: 13
HiC_scaffold_12 65598   1.569e+05       HiC_scaffold_12 HG995277.1      97.139  93349   1981    193     8110585 8203593 7135049 7042051 0.0     Bombus pascuorum        Eukaryota       Bombus pascuorum genome assembly, chromosome: 10
HiC_scaffold_13 30195   8.330e+05       HiC_scaffold_5  OU342930.1      77.744  1464740 254697  45184   1492643 2895399 12199451        13654872        0.0     Bombus terrestris       Eukaryota       Bombus terrestris genome assembly, chromosome: 10
HiC_scaffold_14 30195   8.323e+05       HiC_scaffold_5  OU342930.1      77.735  1464960 254678  44806   1467391 2870452 12199453        13654816        0.0     Bombus terrestris       Eukaryota       Bombus terrestris genome assembly, chromosome: 10
HiC_scaffold_15 65598   1.397e+05       HiC_scaffold_15 HG995279.1      97.155  83013   1850    162     5203366 5286170 5916411 5999119 0.0     Bombus pascuorum        Eukaryota       Bombus pascuorum genome assembly, chromosome: 12
HiC_scaffold_16 30201   2.684e+05       HiC_scaffold_16 OU443163.1      79.008  420201  67792   12639   7096155 7499262 2048536 2465412 0.0     Bombus sylvestris       Eukaryota       Bombus sylvestris genome assembly, chromosome: 23
HiC_scaffold_17 30195   8.316e+05       HiC_scaffold_5  OU342930.1      77.725  1465161 254901  44829   1313971 2717408 12199453        13654875        0.0     Bombus terrestris       Eukaryota       Bombus terrestris genome assembly, chromosome: 10
HiC_scaffold_18 30195   97522   HiC_scaffold_18 OU342933.1      79.851  140541  22765   2366    3697458 3833166 5531383 5671202 0.0     Bombus terrestris       Eukaryota       Bombus terrestris genome assembly, chromosome: 13
HiC_scaffold_19 77635   6613    HiC_scaffold_19 CP062939.1      89.308  5387    421     108     77385   82679   2491536 2496859 0.0     Bifidobacterium subtile Bacteria        Bifidobacterium subtile strain KCTC 3272 chromosome, complete genome
HiC_scaffold_2  30195   8.389e+05       HiC_scaffold_5  OU342930.1      77.811  1464372 254301  44859   1556250 2958940 12199451        13654875        0.0     Bombus terrestris       Eukaryota       Bombus terrestris genome assembly, chromosome: 10
HiC_scaffold_20 30194   3363    HiC_scaffold_20 OV883986.1      70.383  22838   5315    1245    153891  175970  175833  153687  0.0     Bombus pratorum Eukaryota       Bombus pratorum genome assembly, chromosome: 4
HiC_scaffold_21 2575261 3297    HiC_scaffold_21 MK770119.1      76.065  6622    1422    137     38002   44539   40324   33782   0.0     Pantoea phage vB_PagS_AAS21     Viruses Pantoea phage vB_PagS_AAS21, complete genome
HiC_scaffold_22 77635   3522    HiC_scaffold_22 CP062939.1      76.201  7034    1464    180     26126   33050   2177205 2170273 0.0     Bifidobacterium subtile Bacteria        Bifidobacterium subtile strain KCTC 3272 chromosome, complete genome
HiC_scaffold_23 1603886 1463    HiC_scaffold_23 CP062948.1      82.360  1695    287     11      61687   63373   382523  380833  0.0     Bifidobacterium lemurum Bacteria        Bifidobacterium lemurum strain DSM 28807 chromosome, complete genome
HiC_scaffold_24 33905;1254439   8684    HiC_scaffold_24 LR698979.1      82.969  9753    1526    100     6404    16096   1767098 1776775 0.0     Bifidobacterium thermophilum;Bifidobacterium thermophilum RBL67 Bacteria        Bifidobacterium thermophilum isolate MGYG-HGUT-02334 genome assembly, chromosome: 1
HiC_scaffold_25 2170413 1070    HiC_scaffold_25 BK058804.1      79.293  1555    302     19      29677   31221   27374   25830   0.0     Caudoviricetes sp.      Viruses MAG TPA_asm: Siphoviridae sp. isolate ctatL78, partial genome
HiC_scaffold_26 77635   1330    HiC_scaffold_26 CP062939.1      84.996  1313    193     4       14172   15482   1575111 1576421 0.0     Bifidobacterium subtile Bacteria        Bifidobacterium subtile strain KCTC 3272 chromosome, complete genome
HiC_scaffold_27 2020965 1424    HiC_scaffold_27 CP071591.1      87.066  1268    154     9       55504   56766   1997930 1996668 0.0     Bifidobacterium imperatoris     Bacteria        Bifidobacterium imperatoris strain JCM 32708 chromosome
HiC_scaffold_28 1682    1216    HiC_scaffold_28 CP062951.1      78.239  1976    374     49      37505   39447   2033015 2031063 0.0     Bifidobacterium longum subsp. infantis  Bacteria        Bifidobacterium longum subsp. infantis strain JCM 11347 chromosome, complete genome
HiC_scaffold_29 2170413 665     HiC_scaffold_29 MN855933.1      80.154  912     164     15      33439   34346   4539    5437    0.0     Caudoviricetes sp.      Viruses MAG: Siphoviridae sp. isolate 66, complete genome
HiC_scaffold_3  30195   8.369e+05       HiC_scaffold_5  OU342930.1      77.829  1458484 252913  44367   1587011 2984197 12199451        13648782        0.0     Bombus terrestris       Eukaryota       Bombus terrestris genome assembly, chromosome: 10
HiC_scaffold_30 42906   494     HiC_scaffold_30 MT039150.1      89.922  387     33      6       30562   30945   3442    3825    5.63e-133       Serratia entomophila    Bacteria        Serratia entomophila strain AGR_345 plasmid unnamed3, complete sequence
HiC_scaffold_31 30191   1731    HiC_scaffold_31 OU427027.1      73.665  5111    1074    247     1       4974    108828  113803  0.0     Bombus hypnorum Eukaryota       Bombus hypnorum genome assembly, chromosome: 8
HiC_scaffold_32 65598   4348    HiC_scaffold_32 HG995279.1      82.419  5011    853     28      1       4999    3703618 3698624 0.0     Bombus pascuorum        Eukaryota       Bombus pascuorum genome assembly, chromosome: 12
HiC_scaffold_4  30195   8.357e+05       HiC_scaffold_5  OU342930.1      77.777  1465064 253960  44961   1515797 2918912 12199451        13654846        0.0     Bombus terrestris       Eukaryota       Bombus terrestris genome assembly, chromosome: 10
HiC_scaffold_5  30195   8.396e+05       HiC_scaffold_5  OU342930.1      77.816  1464352 254517  44864   1576486 2979461 12199451        13654846        0.0     Bombus terrestris       Eukaryota       Bombus terrestris genome assembly, chromosome: 10
HiC_scaffold_6  30195   8.343e+05       HiC_scaffold_5  OU342930.1      77.753  1464541 255202  44976   1506560 2909665 12199451        13654816        0.0     Bombus terrestris       Eukaryota       Bombus terrestris genome assembly, chromosome: 10
HiC_scaffold_7  30201   3.281e+05       HiC_scaffold_7  OU443163.1      81.014  423439  71257   4891    16360989        16777126        2070630 2492233 0.0     Bombus sylvestris       Eukaryota       Bombus sylvestris genome assembly, chromosome: 23
HiC_scaffold_8  30195   2.485e+05       HiC_scaffold_8  OU342938.1      77.163  436580  93923   4822    4642707 5076045 985559  1419600 0.0     Bombus terrestris       Eukaryota       Bombus terrestris genome assembly, chromosome: 18
HiC_scaffold_9  30195   2.083e+05       HiC_scaffold_9  OU342933.1      81.990  248578  41778   2209    2395811 2642811 5634571 5387409 0.0     Bombus terrestris       Eukaryota       Bombus terrestris genome assembly, chromosome: 13
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


Output from blobtols
```
"sel","_id","gc","length","30_537121167.hifi_reads.fastq_ConcatDahlbomiiGenome.FINAL_minimap2_sorted_cov","bestsumorder_phylum","id"
"","0","0.3688","32487547","27.5395","Arthropoda","HiC_scaffold_1"
"","1","0.3525","24216648","25.1137","Arthropoda","HiC_scaffold_2"
"","2","0.3563","19734131","25.1221","Arthropoda","HiC_scaffold_3"
"","3","0.3661","19374274","26.1943","Arthropoda","HiC_scaffold_4"
"","4","0.3644","18741811","24.9687","Arthropoda","HiC_scaffold_5"
"","5","0.3518","18407461","25.2911","Arthropoda","HiC_scaffold_6"
"","6","0.3762","17299723","25.66","Arthropoda","HiC_scaffold_7"
"","7","0.3667","16961761","26.6204","Arthropoda","HiC_scaffold_8"
"","8","0.3699","15367549","27.113","Arthropoda","HiC_scaffold_9"
"","9","0.3669","13970555","27.5592","Arthropoda","HiC_scaffold_10"
"","10","0.3754","13034363","27.5943","Arthropoda","HiC_scaffold_11"
"","11","0.3882","12441523","28.7737","Arthropoda","HiC_scaffold_12"
"","12","0.3673","12007047","26.8592","Arthropoda","HiC_scaffold_13"
"","13","0.3658","11734535","31.3478","Arthropoda","HiC_scaffold_14"
"","14","0.3861","11609214","26.4285","Arthropoda","HiC_scaffold_15"
"","15","0.361","7505742","24.687","Arthropoda","HiC_scaffold_16"
"","16","0.352","5244765","86.2438","Arthropoda","HiC_scaffold_17"
"","17","0.3798","3846188","23.5803","Arthropoda","HiC_scaffold_18"
"","18","0.5466","228834","5.1966","Actinobacteria","HiC_scaffold_19"
"","19","0.4487","178907","7.4422","Arthropoda","HiC_scaffold_20"
"","20","0.393","113038","21.8358","Uroviricota","HiC_scaffold_21"
"","21","0.5598","108254","4.436","Actinobacteria","HiC_scaffold_22"
"","22","0.5339","103634","4.5048","Actinobacteria","HiC_scaffold_23"
"","23","0.5584","91059","5.1736","Actinobacteria","HiC_scaffold_24"
"","24","0.5537","86719","3.8782","Uroviricota","HiC_scaffold_25"
"","25","0.5381","70448","4.5279","Actinobacteria","HiC_scaffold_26"
"","26","0.5313","59575","4.5306","Actinobacteria","HiC_scaffold_27"
"","27","0.5983","39470","4.1672","Actinobacteria","HiC_scaffold_28"
"","28","0.4656","39391","44.1131","Uroviricota","HiC_scaffold_29"
"","29","0.413","33372","1.5051","Proteobacteria","HiC_scaffold_30"
"","30","0.435","5000","0","Arthropoda","HiC_scaffold_31"
"","31","0.276","5000","8.9968","Arthropoda","HiC_scaffold_32"
```

Kept scaffolds 1-18,20 and 32, as the others were contaminants or lacked coverage (scaffold_31). 

```
GC	Length	Coverage	bestsumorder_phylum	ID
0.393	113038	21.8358	Uroviricota	HiC_scaffold_21
0.5537	86719	3.8782	Uroviricota	HiC_scaffold_25
0.4656	39391	44.1131	Uroviricota	HiC_scaffold_29

3 contigs at a length of 239,148

GC	Length	Coverage	bestsumorder_phylum	ID
0.5466	228834	5.1966	Actinobacteria	HiC_scaffold_19
0.5598	108254	4.436	Actinobacteria	HiC_scaffold_22
0.5339	103634	4.5048	Actinobacteria	HiC_scaffold_23
0.5584	91059	5.1736	Actinobacteria	HiC_scaffold_24
0.5381	70448	4.5279	Actinobacteria	HiC_scaffold_26
0.5313	59575	4.5306	Actinobacteria	HiC_scaffold_27
0.5983	39470	4.1672	Actinobacteria	HiC_scaffold_28
0.413	33372	1.5051	Proteobacteria	HiC_scaffold_30

8 contigs at a length of 734,646


Bifidobacterium accidental integration into HiC_scaffold_2

HiC_scaffold_2:24097917-24216648 
awk '/^[^>]/ {total += length($0); gc += gsub(/[GC]/, "", $0)} END {printf "GC percent: %.2f%%\n", (gc / total) * 100}' BacterialInsertion.fasta
GC percent: 54.91%

So 9 contigs totaling 853,377 bp were of bacterial origin
```
Final Table of Contamination

| ID                    | GC     | Length | Coverage | Best Sum Order Phylum | Best Blast Hit Genus Species       |
|-----------------------|--------|--------|----------|-----------------------|-----------------------------------|
| HiC_scaffold_21       | 0.393  | 113038 | 21.8358  | Uroviricota           | Pantoea phage vB_PagS_AAS21      |
| HiC_scaffold_25       | 0.5537 | 86719  | 3.8782   | Uroviricota           | Caudoviricetes                    |
| HiC_scaffold_29       | 0.4656 | 39391  | 44.1131  | Uroviricota           | Caudoviricetes                    |
| HiC_scaffold_19       | 0.5466 | 228834 | 5.1966   | Actinobacteria        | Bifidobacterium subtile           |
| HiC_scaffold_22       | 0.5598 | 108254 | 4.436    | Actinobacteria        | Bifidobacterium subtile           |
| HiC_scaffold_23       | 0.5339 | 103634 | 4.5048   | Actinobacteria        | Bifidobacterium lemurum           |
| HiC_scaffold_24       | 0.5584 | 91059  | 5.1736   | Actinobacteria        | Bifidobacterium thermophilum      |
| HiC_scaffold_26       | 0.5381 | 70448  | 4.5279   | Actinobacteria        | Bifidobacterium subtile           |
| HiC_scaffold_27       | 0.5313 | 59575  | 4.5306   | Actinobacteria        | Bifidobacterium imperatoris       |
| HiC_scaffold_28       | 0.5983 | 39470  | 4.1672   | Actinobacteria        | Bifidobacterium longum            |
| HiC_scaffold_30       | 0.413  | 33372  | 1.5051   | Proteobacteria        | Serratia entomophila              |
| HiC_scaffold_2_split  | 0.5491 | 118731 | N/A      | Actinobacteria        | Bifidobacterium bifidum           |


