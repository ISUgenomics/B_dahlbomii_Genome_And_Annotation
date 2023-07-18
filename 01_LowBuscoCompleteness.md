#  Why is the busco such a low score?

# map reads
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/10_FindBuscos

ln -s ../07_braker/SoftmaskedB_dahlbomiiGenome.FINAL.fasta
ln -s /work/gif3/data/amytoth/Dovetail_data_2021/Bombus_dahlbomii/PacBio_FASTQ/30_537121167.hifi_reads.fastq.gz

sh runMinimap.sh 30_537121167.hifi_reads.fastq.gz SoftmaskedB_dahlbomiiGenome.FINAL.fasta
```

Results
```
# reads available
zcat 30_537121167.hifi_reads.fastq.gz |wc -l |awk '{print $1/4}'                                                641888
#Reads mapped with 5% mismatch allowed
less 30_537121167.hifi_reads.fastq_SoftmaskedB_dahlbomiiGenome.FINAL_minimap2.paf |cut -f 1 |uniq|wc -l
623870


623870/641888 = 97.19%
```

### Assemble unmapped reads
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/10_FindBuscos
less 30_537121167.hifi_reads.fastq_SoftmaskedB_dahlbomiiGenome.FINAL_minimap2.paf |cut -f 1 |uniq >MappedReads.list
zcat 30_537121167.hifi_reads.fastq.gz |awk 'NR%4==1 {print substr($1,2)}' >TotalReads.list
cat TotalReads.list MappedReads.list |sort |uniq -c |awk '$1==1 {print $2}' >UnmappedReads.list
seqtk subseq 30_537121167.hifi_reads.fastq.gz UnmappedReads.list >UnmappedReads.fastq

mkdir 02_Flye
ln -s UnmappedReads.fastq
ml miniconda3; source activate flye;flye  --pacbio-hifi UnmappedReads.fastq -o /work/gif/remkv6/Toth/12_Bombus_dahlbomii/10_FindBuscos/02_Flye -t 36
```
results
```
new_Assemblathon.pl assembly.fasta

---------------- Information for assembly 'assembly.fasta' ----------------


                                         Number of scaffolds        197
                                     Total size of scaffolds    9722856
                                            Longest scaffold     341424
                                           Shortest scaffold      10537
                                 Number of scaffolds > 1K nt        197 100.0%
                                Number of scaffolds > 10K nt        197 100.0%
                               Number of scaffolds > 100K nt         17   8.6%
                                 Number of scaffolds > 1M nt          0   0.0%
                                Number of scaffolds > 10M nt          0   0.0%
                                          Mean scaffold size      49355
                                        Median scaffold size      33800
                                         N50 scaffold length      68968
                                          L50 scaffold count         46
                                         n90 scaffold length      23543
                                          L90 scaffold count        144
                                                 scaffold %A      30.23
                                                 scaffold %C      19.77
                                                 scaffold %G      19.68
                                                 scaffold %T      30.32
                                                 scaffold %N       0.00
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0
```

### Run megablast on assembly
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/10_FindBuscos/03_Megablast
ln -s ../02_Flye/assembly.fasta
cp ~/common_scripts/runMegablast2nt.sh .

echo "sh runMegablast2nt.sh assembly.fasta" >Dhalblast.sh

# how many were bombus
less assembly.vs.nt.cul10.1e3.megablast.out |sort -k1,1 -u |grep  "Bombus" |wc
    185    4255   31677

# how many were something else (bifidobacteria and caudoviricetes)
less assembly.vs.nt.cul10.1e3.megablast.out |sort -k1,1 -u |grep -v "Bombus" |wc
     12     298    2365
```

### Run compleasm on the assembly
```
sh runCompleasm.sh assembly.fasta eukaryota

#!/bin/bash
#runCompleasm.sh
#Here is how to run this script: sh runCompleasm Genome.fasta Lineage
Genome="$1"
Lineage="$2"
ml py-pandas/2.0.1-py310-ujlqazv
compleasm.py download ${Lineage}
compleasm.py run \
-t 36 \
-l ${Lineage} \
-a  ${Genome} \
-o ${Lineage}_${Genome%.*}_Compleasm



S:10.20%, 26
D:0.00%, 0
F:0.00%, 0
I:0.00%, 0
M:89.80%, 229
N:255

```

# Run compleasm to see if the busco scores are accurate
```
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/10_FindBuscos/01_MiniBusco
ml py-pandas/2.0.1-py310-ujlqazv
wget https://github.com/huangnengCSU/compleasm/releases/download/v0.2.2/compleasm-0.2.2_x64-linux.tar.bz2
tar -jxvf compleasm-0.2.2_x64-linux.tar.bz2

./compleasm_kit/compleasm.py download eukaryota
compleasm_kit/compleasm.py run -t16 -l eukaryota -a ../SoftmaskedB_dahlbomiiGenome.FINAL.fasta -o EukOut



S:90.20%, 230
D:0.00%, 0
F:0.39%, 1
I:0.00%, 0
M:9.41%, 24
N:255

Yeah, the busco scores are correct.  And not even a fragment present...

```

### Compare compleasm between unmappable reads assembly and final dahlbomii assembly
```
#How many are not missing in unmappable reads assembly that are missing from the dahlbomii final  assembly?
less full_table.tsv |awk 'NR>1 && $2!="Missing"' |cut -f 1 |cat - <(less ../../../../01_MiniBusco/EukOut/eukaryota_odb10/full_table.tsv |awk 'NR>1 && $2=="Missing" {print $1}' ) |sort|uniq -c |awk '$1==2' |wc
     21      42     449
21/24 missing were found this way.
```