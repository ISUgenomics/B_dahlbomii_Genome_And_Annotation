
# Create a functional annotation for B. dahlbomii

##### Blast to NR and Swissprot
```
/work/gif3/masonbrink/Toth/01_dahlbomii

cp -rf ../../Baum/01_ReannotateAllSCNGenomes/02_OP50/01_Mikado/02_SwissProt/runMegablastSP.sh .
cp -rf ../../Baum/01_ReannotateAllSCNGenomes/02_OP50/01_Mikado/03_NR/runMegablastNR.sh .

ln -s /work/gif/remkv6/Toth/12_Bombus_dahlbomii/13_ReAnnotate/braker/Dahlbomii_Proteins.fasta
sed -i '/^[^>]/s/\.//' Dahlbomii_Proteins.fasta
cp -rf  /work/gif/remkv6/Toth/12_Bombus_dahlbomii/13_ReAnnotate/braker/SortedB_dahlbomiiBraker.gff3.gz .
gunzip SortedB_dahlbomiiBraker.gff3.gz

echo "sh runMegablastNR.sh Dahlbomii_Proteins.fasta" >NRBlast.sh
echo "sh runMegablastSP.sh Dahlbomii_Proteins.fasta" >SPBlast.sh

#compile best hit annotations
cat Dahlbomii_Proteins.vs.sp.diamond.out |grep -v -i -e "hypothetical" -e "unnamed" -e "uncharacterized" |sort -k1,1V -u >DefinedSPGenes.tab
cat DefinedSPGenes.tab <(sort -u -k1,1V Dahlbomii_Proteins.vs.sp.diamond.out ) |sort -k1,1V -u |cut -f 1,13 | sed  's/ RecName: Full=/ /g' |sed  's/;.*\[/\ [/g' >AllSPAnnotations.tab


#compile best hit annotations
cat Dahlbomii_Proteins.vs.nr.diamond.out |grep -v -i -e "hypothetical" -e "unnamed" -e "uncharacterized" |sort -k1,1V -u >DefinedNRGenes.tab
cat DefinedNRGenes.tab <(sort -u -k1,1V Dahlbomii_Proteins.vs.nr.diamond.out ) |sort -k1,1V -u |cut -f 1,13 | sed  's/ RecName: Full=/ /g' |sed  's/;.*\[/\ [/g' >AllNRAnnotations.tab

sed -i 's/\t/\tNRBLAST_/1' AllNRAnnotations.tab
sed -i 's/\t/\tSPBLAST_/1' AllSPAnnotations.tab
```

##### Combine Annotations
```
 awk -F"\t" '{arr[$1]=arr[$1] ";" $2}END{for(i in arr)print i,arr[i]}' AllNRAnnotations.tab AllSPAnnotations.tab |sed 's/;//1' |sed 's/;/#/g' >CombinedNRandSPAnnotations.tab

sed -i 's/ /\t/1' CombinedNRandSPAnnotations.tab
```



### Generate statistics on genes and annotations
```
/work/gif/remkv6/Baum/16_NewTN10Genome/22_FunctionalAnnotation
/work/gif/remkv6/Olsen/Bison/38_FunctionalAnnotation
#how many mRNAs and genes are there in the genome?

Count:  15,767
Count:  19,737

#total number of annotations obtained.
wc CombinedNRandSPAnnotations.tab
  15865  205986 2293875 CombinedNRandSPAnnotations.tab


#what percentages of the proteins have a functional annotation?
15865/19,737 = 80.4%

#how many transposon associated transcripts?
grep -i -e "transposase" -e "transcriptase" -e "retroelement" -e "helitron" -e "reverse" -e "helitron" -e "transposon"  CombinedNRandSPAnnotations.tab |wc
     38     420    4663
```


### Attach functional annotations to mRNA in gff3s
```
awk -F"\t" '{arr[$1]=arr[$1] ";" $2}END{for(i in arr)print i,arr[i]}'  <(less SortedB_dahlbomiiBraker.gff3 |awk 'substr($1,1,1)!="#"' |awk '$3=="mRNA"' |sed 's/;/\t/1'  |sed 's/ID=/ID=\t/g' |sed 's/;/;\t/1' |awk '{print $10"\t"$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' |sed 's/ //9' |sed 's/ /;/9'  ) CombinedNRandSPAnnotations.tab |sed 's/ ;H/\tH/g' |cut -f 2- |sed 's/;;/;Note=/g' |sed 's/ /\t/1'|sed 's/ /\t/1'|sed 's/ /\t/1'|sed 's/ /\t/1'|sed 's/ /\t/1'|sed 's/ /\t/1'|sed 's/ /\t/1'|sed 's/ /\t/1' >mRNAsAnnotated.gff3


less mRNAsAnnotated.gff3 |cut -f 9 |sed 's/ID=//g' |sed 's/Parent=//g' |sed 's/;/\t/1' |sed 's/;/\t/1'|sort -k3,3Vr |sort -k2,2 -u |sed '/^$/d' |cut -f 2-  >Gene2Annotation.list

awk -F"\t" '{arr[$1]=arr[$1] ";" $2}END{for(i in arr)print i,arr[i]}'  <(less SortedB_dahlbomiiBraker.gff3 |awk 'substr($1,1,1)!="#"' |awk '$3=="gene"'  |sed 's/;/\t/1' |sed 's/ID=/ID=\t/g'| cut -f 1-10 |awk '{print $10"\t"$1,$2,$3,$4,$5,$6,$7,$8,$9,$10";"}' |sed 's/ //9'  )  Gene2Annotation.list  |sed 's/ ;/\t/1' |cut -f 2- |sed 's/ /\t/1'|sed 's/ /\t/1'|sed 's/ /\t/1'|sed 's/ /\t/1'|sed 's/ /\t/1'|sed 's/ /\t/1'|sed 's/ /\t/1'|sed 's/ /\t/1' |sed 's/;;/;/g' >genesAnnotated.gff3

awk '$3!="gene" && $3!="mRNA" && substr($1,1,1)!="#"' SortedB_dahlbomiiBraker.gff3 >OtherFeatures.gff3
cat OtherFeatures.gff3 genesAnnotated.gff3 mRNAsAnnotated.gff3  >UnorderedMerge.gff3
gff3sort.pl --precise --chr_order natural UnorderedMerge.gff3  >SortedFunctionalAnnotationBdahlbomii.gff3
ml tabix
bgzip SortedFunctionalAnnotationBdahlbomii.gff3
tabix -p gff SortedFunctionalAnnotationBdahlbomii.gff3.gz

```
