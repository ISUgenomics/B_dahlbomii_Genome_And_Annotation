#  Why is the busco such a low score?

# map long reads
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

#mapping percentage
623870/641888 = 97.19%

641888- 623870 = 18,018 reads 
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


What would be the complete percent if these contigs were added to the genome?
(21+228)/255 = 97.6% complete buscos
Which contigs have these?
/work/gif/remkv6/Toth/12_Bombus_dahlbomii/10_FindBuscos/03_Megablast
less ../04_Compleasm/02_UnmappedGenome/eukaryota_assembly_Compleasm/eukaryota_odb10/full_table.tsv |awk 'NR>1 && $2!="Missing" {print $3}' |grep -w -f - <(less assembly.vs.nt.cul10.1e3.megablast.out |sort -k1,1 -u |grep  "Bombus") |less

contig_112      65598   28439   contig_112      HG995279.1      96.965  17004   405     57      2       16965   12613214        12596282        0.0     Bombus pascuorum        Eukaryota       Bombus pascuorum genome assembly, chromosome: 12
contig_118      30201   30180   contig_118      OU443161.1      95.158  19250   711     90      4036    23187   3087523 3068397 0.0     Bombus sylvestris       Eukaryota       Bombus sylvestris genome assembly, chromosome: 21
contig_119      207624  23848   contig_119      HG995132.1      94.647  15525   595     101     14697   30068   12894757        12879316        0.0     Bombus campestris       Eukaryota       Bombus campestris genome assembly, chromosome: 7
contig_125      30201   14990   contig_125      OU443152.1      96.732  9026    247     16      78202   87221   381491  372508  0.0     Bombus sylvestris       Eukaryota       Bombus sylvestris genome assembly, chromosome: 13
contig_126      65598   45061   contig_126      HG995279.1      95.765  28125   889     129     23426   51454   11728999        11756917        0.0     Bombus pascuorum        Eukaryota       Bombus pascuorum genome assembly, chromosome: 12
contig_128      65598   45995   contig_128      HG995272.1      96.656  27845   642     72      23921   51696   17506619        17478995        0.0     Bombus pascuorum        Eukaryota       Bombus pascuorum genome assembly, chromosome: 5
contig_141      65598   16829   contig_141      HG995279.1      93.323  11547   530     61      65766   77236   12484544        12473163        0.0     Bombus pascuorum        Eukaryota       Bombus pascuorum genome assembly, chromosome: 12
contig_152      85660   53555   contig_152      HG995195.1      96.299  32802   897     142     9967    42599   2852129 2884782 0.0     Bombus hortorum Eukaryota       Bombus hortorum genome assembly, chromosome: 8
contig_154      65598   47958   contig_154      HG995279.1      94.515  31520   1003    211     127274  158461  10198026        10229151        0.0     Bombus pascuorum        Eukaryota       Bombus pascuorum genome assembly, chromosome: 12
contig_159      30195   58013   contig_159      OU342925.1      94.848  37482   1384    208     1       37306   4793110 4756000 0.0     Bombus terrestris       Eukaryota       Bombus terrestris genome assembly, chromosome: 5
contig_163      85660   48944   contig_163      HG995197.1      95.411  30875   1177    73      25896   56632   11856111        11825339        0.0     Bombus hortorum Eukaryota       Bombus hortorum genome assembly, chromosome: 10
contig_169      65598   54979   contig_169      HG995279.1      96.581  33342   840     84      46708   79945   8040835 8073980 0.0     Bombus pascuorum        Eukaryota       Bombus pascuorum genome assembly, chromosome: 12
contig_17       85660   50196   contig_17       HG995195.1      95.133  32070   1152    169     67361   99180   2411106 2443016 0.0     Bombus hortorum Eukaryota       Bombus hortorum genome assembly, chromosome: 8
contig_175      30201   19064   contig_175      OU443161.1      94.567  12480   433     76      51483   63841   3254983 3267338 0.0     Bombus sylvestris       Eukaryota       Bombus sylvestris genome assembly, chromosome: 21
contig_177      85660   27481   contig_177      HG995195.1      94.235  18196   714     156     8176    26225   3898724 3880718 0.0     Bombus hortorum Eukaryota       Bombus hortorum genome assembly, chromosome: 8
contig_189      207624  24701   contig_189      HG995140.1      94.622  16085   638     81      16613   32598   4823929 4839885 0.0     Bombus campestris       Eukaryota       Bombus campestris genome assembly, chromosome: 15
contig_42       30191   16907   contig_42       OU427026.1      93.961  11260   550     45      2670    13889   7110579 7121748 0.0     Bombus hypnorum Eukaryota       Bombus hypnorum genome assembly, chromosome: 7

Longest contig is 58 kb.
```

