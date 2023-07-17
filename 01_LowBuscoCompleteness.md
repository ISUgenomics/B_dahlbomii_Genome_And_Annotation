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

### run compleasm
```

```


