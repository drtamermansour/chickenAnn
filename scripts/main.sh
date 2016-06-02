cd chickenAnn
chick_ann=$(pwd)
mkdir -p $chick_ann/{scripts,resources,output,annotation}

## read the user configurations
#source $chick_ann/user_config.txt
#cat $chick_ann/user_config.txt

## create a config file to contain all the pathes to be used by all pipelines
> $chick_ann/config.txt
echo "script_path=$chick_ann/scripts" >> $chick_ann/config.txt
echo "resources=$chick_ann/resources" >> $chick_ann/config.txt
echo "output=$chick_ann/output" >> $chick_ann/config.txt
echo "ann=$chick_ann/annotation" >> $chick_ann/config.txt
source $chick_ann/config.txt

cd $resources
wget ftp://ftp.ncbi.nih.gov/genomes/Gallus_gallus/comparisons/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt.gz
gunzip GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt.gz

cd $output
## classes
tail -n+2 $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '{A[$1]++}END{for(i in A)print i,A[i]}' > gene_category.freq.report
tail -n+2 $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '{A[$3]++}END{for(i in A)print i,A[i]}' > current_gene_biotype.freq.report

## protein coding genes 
cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '{if($3 == "protein_coding") print $2}' | sort | uniq | wc -l ## 19060
cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '{if($9 == "protein_coding") print $8}' | sort | uniq | wc -l ## 17139

## protein coding transcripts
cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '{if($3 == "protein_coding") print $15}' | sort | uniq | wc -l ## 48419
cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '{if($9 == "protein_coding") print $18}' | sort | uniq | wc -l ## 33105

## proteins
cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '{if($3 == "protein_coding") print $16}' | sort | uniq | wc -l ## 46334
cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '{if($9 == "protein_coding") print $19}' | sort | uniq | wc -l ## 32118

## identification of novel protein coding genes
cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '(($3 == "protein_coding") && ($9 == "protein_coding") && ($2 == $8))' > noChange.txt ## genes that did not change
cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '{if(($3 == "protein_coding") && ($9 == "protein_coding") && ($2 == $8)) print $2}' | sort | uniq | wc -l ## 15260 genes that did not change
cat noChange.txt | awk -F $'\t' '{A[$1]++}END{for(i in A)print i,A[i]}' > noChange.1.report  ## gene categories 

#cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '(($3 != "protein_coding") && ($9 == "protein_coding"))' > ptnChange_a.txt ## lost ptn coding genes
#cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '{if(($3 != "protein_coding") && ($9 == "protein_coding")) print $8}' | sort | uniq | wc -l ## 1727
#cat ptnChange_a.txt | awk -F $'\t' '{A[$1]++}END{for(i in A)print i,A[i]}' > ptnChange_a.1.report  ## gene categories 

cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '(($3 == "protein_coding") && ($9 != "protein_coding"))' > ptnChange_b.txt ## new ptncoding genes  
cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '{if(($3 == "protein_coding") && ($9 != "protein_coding")) print $2}' | sort | uniq | wc -l ## 3716
cat ptnChange_b.txt | awk -F $'\t' '{A[$1]++}END{for(i in A)print i,A[i]}' > ptnChange_b.1.report  ## gene categories 
cat ptnChange_b.txt | awk -F $'\t' '(($1 == "Current-novel") || ($1 == "Current-unmapped") || ($1 == "Current-other"))' > ptnChange_b_sig.txt
cat ptnChange_b.txt | awk -F $'\t' '{if(($1 == "Current-novel") || ($1 == "Current-unmapped") || ($1 == "Current-other")) print $2}' | sort | uniq | wc -l ## 3459
cat ptnChange_b.txt | awk -F $'\t' '(($1 != "Current-novel") && ($1 != "Current-unmapped") && ($1 != "Current-other"))' >  ptnChange_b_insig.txt
cat ptnChange_b.txt | awk -F $'\t' '{if(($1 != "Current-novel") && ($1 != "Current-unmapped") && ($1 != "Current-other")) print $2}' | sort | uniq | wc -l ## 257 (has 44 repeated IDs in noChange.xtx)

cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '(($3 == "protein_coding") && ($9 == "protein_coding") && ($2 != $8))' > ptnChange_c.txt ## genes that ware and are still ptncoding genes  
cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '{if(($3 == "protein_coding") && ($9 == "protein_coding") && ($2 != $8)) print $2}' | sort | uniq | wc -l ## 285 (has repeated IDs in noChange.txt, ptnChange_b_sig.txt, and ptnChange_b_insig.txt)
cat ptnChange_c.txt | awk -F $'\t' '{A[$1]++}END{for(i in A)print i,A[i]}' > ptnChange_c.1.report  ## gene categories 

#cat noChange.txt ptnChange_b_insig.txt | awk -F $'\t' '{print $2}' | sort | uniq | wc -l  ## 15473 (15517) => 44
#cat noChange.txt ptnChange_c.txt | awk -F $'\t' '{print $2}' | sort | uniq | wc -l  ## 15398 (15545) => 147
#cat ptnChange_b_sig.txt ptnChange_c.txt | awk -F $'\t' '{print $2}' | sort | uniq | wc -l  ## 3743 (3744) => 1
#cat ptnChange_b_insig.txt ptnChange_c.txt | awk -F $'\t' '{print $2}' | sort | uniq | wc -l  ## 528 (542) => 14
#cat ptnChange_b.txt ptnChange_c.txt | awk -F $'\t' '{print $2}' | sort | uniq | wc -l  ## 3986


## define the lists of improtant new genes 
head -n1 $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt > novels.txt
cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk '$1 == "Current-novel"' >> novels.txt
tail -n+2 novels.txt | awk -F $'\t' '{A[$3]++}END{for(i in A)print i,A[i]}' > novels.freq.report  ## classify by "current gene biotype"
tail -n+2 novels.txt | awk -F $'\t' '$3 == "protein_coding"' > novels.protein_coding.txt
cat novels.protein_coding.txt | awk -F $'\t' '{print $2}' | sort | uniq > novels.protein_coding.geneIDs  ## 980
cat novels.protein_coding.txt | awk -F $'\t' '{print $15}' | sort | uniq > novels.protein_coding.transIDs  ## 1436
cat novels.protein_coding.txt | awk -F $'\t' '{if($16!="NA") print $16}' | sort | uniq > novels.protein_coding.ptnIDs  ## 1357 ## 79 transcript has an "NA" protein
cat novels.protein_coding.txt | awk -F $'\t' 'BEGIN{OFS="\t";} {print $2,$15,$16}' > novels.protein_coding.map

head -n1 $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt > unmapped.txt
cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk '$1 == "Current-unmapped"' >> unmapped.txt
tail -n+2 unmapped.txt | awk -F $'\t' '{A[$3]++}END{for(i in A)print i,A[i]}' > unmapped.freq.report  ## classify by "current gene biotype"
tail -n+2 unmapped.txt | awk -F $'\t' '$3 == "protein_coding"' > unmapped.protein_coding.txt
cat unmapped.protein_coding.txt | awk -F $'\t' '{print $2}' | sort | uniq > unmapped.protein_coding.geneIDs  ## 2129
cat unmapped.protein_coding.txt | awk -F $'\t' '{print $15}' | sort | uniq > unmapped.protein_coding.transIDs  ## 2675
cat unmapped.protein_coding.txt | awk -F $'\t' '{if($16!="NA") print $16}' | sort | uniq > unmapped.protein_coding.ptnIDs  ## 2561
cat unmapped.protein_coding.txt | awk -F $'\t' 'BEGIN{OFS="\t";} {print $2,$15,$16}' > unmapped.protein_coding.map ## 2675

head -n1 $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt > other.txt
cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk '$1 == "Current-other"' >> other.txt
tail -n+2 other.txt | awk -F $'\t' '{A[$3]++}END{for(i in A)print i,A[i]}' > other.freq.report  ## classify by "current gene biotype"
tail -n+2 other.txt | awk -F $'\t' '$3 == "protein_coding"' > other.protein_coding.txt
cat other.protein_coding.txt | awk -F $'\t' '{print $2}' | sort | uniq > other.protein_coding.geneIDs  ## 350
cat other.protein_coding.txt | awk -F $'\t' '{print $15}' | sort | uniq > other.protein_coding.transIDs  ## 522
cat other.protein_coding.txt | awk -F $'\t' '{if($16!="NA") print $16}' | sort | uniq > other.protein_coding.ptnIDs  ## 478
cat other.protein_coding.txt | awk -F $'\t' 'BEGIN{OFS="\t";} {print $2,$15,$16}' > other.protein_coding.map ## 522

#head -n1 $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt > changeCoding.txt ## 7100
#cat $resources/GCF_000002315.4_Gallus_gallus-5.0_compare_prev.txt | awk -F $'\t' '(($1 != "Current-novel") && ($3 == "protein_coding") && ($9 != "protein_coding"))' >> changeCoding.txt
#tail -n+2 changeCoding.txt | awk -F $'\t' '{A[$1]++}END{for(i in A)print i,A[i]}' > changeCoding_geneCat.freq.report 	     ## classify acc to "gene category"
#tail -n+2 changeCoding.txt | awk -F $'\t' '{A[$14]++}END{for(i in A)print i,A[i]}' > changeCoding_transCat.freq.report  ## classify acc to "transcript category"
#tail -n+2 changeCoding.txt | awk -F $'\t' '{A[$9]++}END{for(i in A)print i,A[i]}' > changeCoding_geneBiotype.freq.report  ## classify acc to "previous gene biotype"
#tail -n+2 changeCoding.txt | awk -F $'\t' '{print $15}' | sort | uniq > changeCoding.transIDs  ## 3689
#tail -n+2 changeCoding.txt | awk -F $'\t' '{print $2}' | sort | uniq > changeCoding.geneIDs  ## 2736

#tail -n+2 changeCoding.txt | awk -F $'\t' '$9 == "NA"' > changeCoding.novel.txt ## previous gene biotype= Current-other 1044 & Current-unmapped 5350
#cat changeCoding.novel.txt | awk -F $'\t' '{print $15}' | sort | uniq > changeCoding.novel.transIDs  ## 3197
#cat changeCoding.novel.txt | awk -F $'\t' '{print $2}' | sort | uniq > changeCoding.novel.geneIDs  ## 2479
#cat changeCoding.novel.txt | awk -F $'\t' '{if($1 == "Current-other") print $2}' | sort | uniq > changeCoding.novel.CurrentOther.geneIDs  ## 350
#cat changeCoding.novel.txt | awk -F $'\t' '{if($1 == "Current-unmapped") print $2}' | sort | uniq > changeCoding.novel.CurrentUnmapped.geneIDs  ## 2129

#tail -n+2 changeCoding.txt | awk -F $'\t' '$9 != "NA"' > changeCoding.previous.txt ## 705
#cat changeCoding.previous.txt | awk -F $'\t' '{print $15}' | sort | uniq > changeCoding.previous.transIDs  ## 492
#cat changeCoding.previous.txt | awk -F $'\t' '{print $2}' | sort | uniq > changeCoding.previous.geneIDs  ## 257

## find GO terms on https://biodbnet-abcc.ncifcrf.gov/db/dbAnnot.php
## organism (Taxon ID): 9031   && Outputs: GO terms

####################################
cd $resources
#wget ftp://ftp.ncbi.nih.gov/genomes/Gallus_gallus/RNA/rna.fa.gz
#gunzip rna.fa.gz
#module load QIIME/1.8.0
#grep -F -w -f $output/novels.protein_coding.transIDs rna.fa | sed 's/>//' > novels.protein_coding.transIDs2
#filter_fasta.py --input_fasta_fp rna.fa --output_fasta_fp novels.protein_coding.fa --seq_id_fp novels.protein_coding.transIDs2

wget ftp://ftp.ncbi.nih.gov/genomes/Gallus_gallus/protein/protein.fa.gz
gunzip protein.fa.gz
module load QIIME/1.8.0

grep -F -w -f $output/novels.protein_coding.ptnIDs protein.fa | sed 's/>//' > novels.protein_coding.ptnIDs2
filter_fasta.py --input_fasta_fp protein.fa --output_fasta_fp novels.protein_coding.pep --seq_id_fp novels.protein_coding.ptnIDs2
grep "PREDICTED" novels.protein_coding.pep > novels.protein_coding.pep.predicted
grep "^>" novels.protein_coding.pep | grep -v "PREDICTED" > novels.protein_coding.pep.confirmed
cat novels.protein_coding.pep.predicted | awk -F '[|:]' 'BEGIN {OFS="\t"} {print $4,$6}' | sed 's/ \[Gallus gallus\]//g' > novels.protein_coding.pep.map
#cat novels.protein_coding.pep.predicted | awk -F '[|:]' '{print $4}' > novels.protein_coding.pep.map.ids
#cat novels.protein_coding.pep.predicted | awk -F '[|:]' '{print $6}' | sed 's/ \[Gallus gallus\]//g' | sed 's/^[[:space:]]//g' | sed 's/[[:space:]]$//g' > novels.protein_coding.pep.map.desc
#paste novels.protein_coding.pep.map.ids novels.protein_coding.pep.map.desc > novels.protein_coding.pep.map
cat novels.protein_coding.pep.confirmed | awk -F '|' 'BEGIN {OFS="\t"} {print $4,$5}' | sed 's/ \[Gallus gallus\]//g' >> novels.protein_coding.pep.map

grep -F -w -f $output/unmapped.protein_coding.ptnIDs protein.fa | sed 's/>//' > unmapped.protein_coding.ptnIDs2
filter_fasta.py --input_fasta_fp protein.fa --output_fasta_fp unmapped.protein_coding.pep --seq_id_fp unmapped.protein_coding.ptnIDs2
grep "PREDICTED" unmapped.protein_coding.pep > unmapped.protein_coding.pep.predicted
grep "^>" unmapped.protein_coding.pep | grep -v "PREDICTED" > unmapped.protein_coding.pep.confirmed
cat unmapped.protein_coding.pep.predicted | awk -F '[|:]' 'BEGIN {OFS="\t"} {print $4,$6}' | sed 's/ \[Gallus gallus\]//g' > unmapped.protein_coding.pep.map
cat unmapped.protein_coding.pep.confirmed | awk -F '|' 'BEGIN {OFS="\t"} {print $4,$5}' | sed 's/ \[Gallus gallus\]//g'>> unmapped.protein_coding.pep.map

grep -F -w -f $output/other.protein_coding.ptnIDs protein.fa | sed 's/>//' > other.protein_coding.ptnIDs2
filter_fasta.py --input_fasta_fp protein.fa --output_fasta_fp other.protein_coding.pep --seq_id_fp other.protein_coding.ptnIDs2
grep "PREDICTED" other.protein_coding.pep > other.protein_coding.pep.predicted
grep "^>" other.protein_coding.pep | grep -v "PREDICTED" > other.protein_coding.pep.confirmed
cat other.protein_coding.pep.predicted | awk -F '[|:]' 'BEGIN {OFS="\t"} {print $4,$6}' | sed 's/ \[Gallus gallus\]//g' > other.protein_coding.pep.map
cat other.protein_coding.pep.confirmed | awk -F '|' 'BEGIN {OFS="\t"} {print $4,$5}' | sed 's/ \[Gallus gallus\]//g' >> other.protein_coding.pep.map

## prepare the swiss-prot DB
module load BLAST+/2.2.29
mkdir ${resources}/uniprot
cd ${resources}/uniprot
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -dbtype prot
######################################
## Functional annotation
mkdir $output/iprscan_out
cd $output/iprscan_out
module load Java/1.8.0_31
module load InterProScan/5.17-56.0
cp $resources/novels.protein_coding.pep novel.fa
iprscan --goterms --pathways --seqtype p --input novel.fa > novel.iprscan.log
cp $resources/unmapped.protein_coding.pep unmapped.fa 
iprscan --goterms --pathways --seqtype p --input unmapped.fa > unmapped.iprscan.log
cp $resources/other.protein_coding.pep other.fa
iprscan --goterms --pathways --seqtype p --input other.fa > other.iprscan.log

## keep copy of these files for record 
#cp  /opt/software/signalp/4.1--Binary/signalp-4.1.license.txt .
#cp /mnt/research/common-data/Bio/iprscan/interproscan.properties .
for f in *.tsv;do
  echo $f;
  echo "no of annotated genes"
  cat $f | awk -F '\t' '{print $1}' | sort | uniq | wc -l
  echo "frequencies of unique models in memeber databases: check" $f.freq
  cat $f | awk -F $'\t' '{A[$5]++;if(A[$5]==1) B[$4]++}END{for(i in B)print i,B[i]}' | sort -k2,2nr > $f.freq
  echo "no of unique interpro models"
  cat $f | awk -F '\t' '{print $12}' | sort | uniq | wc -l
  echo "no of proteins with every interpro model: check" $f.interpro.freq
  awk -F $'\t' 'BEGIN{OFS="\t";} {print $1,$12,$13;}' $f | sort | uniq | awk -F $'\t' '{A[$2"\t"$3]++}END{for(i in A)print A[i]"\t"i}' | sort -k1,1nr > $f.interpro.freq
  echo "no of genes with GO annotation"
  cat $f | awk -F '\t' '{if($14!="") print $1}' | sort | uniq | wc -l
  echo "no of unique GO annotation terms" 
  cat $f | awk -F '\t' '{if($14!="") print $14}' | tr '|' '\n' | sort | uniq | wc -l
  cat $f | awk -F '\t' '{if($14!="") print $14}' | tr '|' '\n' | sort | uniq > $f.interpro.models
done
## frequencies of unique models in memeber databases for all genes
cat *.tsv | awk -F $'\t' '{A[$5]++;if(A[$5]==1) B[$4]++}END{for(i in B)print i,B[i]}' | sort -k2,2nr > all.freq
## no of unique interpro models for all genes
cat *.tsv | awk -F '\t' '{print $12}' | sort | uniq | wc -l
## no of proteins with every interpro model for all genes
cat *.tsv | awk -F $'\t' 'BEGIN{OFS="\t";} {print $1,$12,$13;}' | sort | uniq | awk -F $'\t' '{A[$2"\t"$3]++}END{for(i in A)print A[i]"\t"i}' | sort -k1,1nr > all.interpro.freq
## no of unique GO annotation terms for all genes
cat *.tsv | awk -F '\t' '{if($14!="") print $14}' | tr '|' '\n' | sort | uniq | wc -l
#########################################
## Run blast aganist swiss-prot DB
mkdir $output/blastp_out
cd $output/blastp_out
for f in $resources/*.protein_coding.pep; do
  newf=$(basename $f .protein_coding.pep)
  cp $f $newf.fa;
  awk '{print $1}' $newf.fa > "$newf"_simple.fa
  sed -i 's/|$//g' "$newf"_simple.fa
  sed -i 's/>.*|/>/g' "$newf"_simple.fa
  perl ${script_path}/splitFasta.pl "$newf"_simple.fa 4  ## 500 for uniprot_uniref90
done

module load BLAST+/2.2.29
for f in subset*_*_simple.fa; do
#  qsub -v input=$f,DB=$resources/uniprot/uniprot_sprot.fasta,label="uniprot_sprot" ${script_path}/blastp.sh;
  blastp -query "$f" -db "$resources/uniprot/uniprot_sprot.fasta" -num_threads 4 -max_target_seqs 20 -outfmt 5 -seg yes -evalue 1e-3 > "$f".uniprot_sprot.blastp.xml
done
#cat subset*_trans_truncated.fa.xml > ../uniprot_sprot.blastx.outfmt6             ## 5614133

cd $output
module load blast2go  ## version 2.5
cp /mnt/research/common-data/Bio/blast2go/b2gPipe.properties .
## java es.blast2go.prog.B2GAnnotPipe -prop b2gPipe.properties <args>
## Maximum usage example including the full functionality of the pipeline:
## java -Xmx1000m -cp *:ext/*: es.blast2go.prog.B2GAnnotPipe -in 10_BlastResults_2011.xml -out results/myproject -prop b2gPipe.properties -annot -dat -img -ips ipsr -annex -goslim -wiki html_template.html -v	
## Minimum uasge example to create only the annotation:
## java -Xmx500m -cp *:ext/*: es.blast2go.prog.B2GAnnotPipe -in 10_BlastResults_2011.xml -out results/myproject -prop b2gPipe.properties -annot
java es.blast2go.prog.B2GAnnotPipe -prop b2gPipe.properties -in novel.fa.uniprot_sprot.blastp.xml -annot -dat -img -annex -goslim -v > blast2go.log 

java -cp *:ext/*: es.blast2go.prog.B2GAnnotPipe -in 10_BlastResults_2011.xml -out results/myproject -prop b2gPipe.properties -v -annot -dat -img -ips ipsr -annex -goslim -wiki html_template.html

##########################################
mkdir $output/ncbiSearchByName
cd $output/ncbiSearchByName
for pepMap in $resources/*.protein_coding.pep.map;do
 geneMap=$output/$(basename $pepMap .pep.map).map
 newMap=$(basename $pepMap).comp
 Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],header=F,row.names=NULL); head(data1); data2=read.table(args[2],header=F,row.names=NULL,sep="\t",quote=""); head(data2); dataMerge=merge(data1,data2,by.x="V3",by.y="V1",all.y=T);write.table(dataMerge[,c(2,3,1,4)],args[3], sep="\t", quote=F, row.names=F, col.names=F);' $geneMap $pepMap $newMap
done

for map in *.protein_coding.pep.map.comp;do  ## 1357
 grep -v "LOW QUALITY PROTEIN" $map > $map.simple  ## 1312
 grep -v "uncharacterized protein LOC" $map.simple > $map.simple2 ## 1056
 sed -i 's/ isoform .*$//g' $map.simple2
 sed -i 's/, partial$//g' $map.simple2
 sed -i 's/-like$//g' $map.simple2
done

for map in *.simple2;do cat $map | awk -F '\t' '{print $4}' | sort | uniq > $map.id; done  ## 490

#while read gene desc;do
#  echo $desc;
#  esearch -db gene -query "$desc [KYWD] AND human [ORGN]" | efetch -format native -mode xml | xtract -pattern Entrezgene -element Gene-track_geneid Org-ref_taxname Gene-ref_locus Gene-ref_desc -group Gene-ref_db -match "Dbtag_db:HGNC" -first Object-id_str -group Gene-ref_syn -sep "|" -element Gene-ref_syn_E -group Gene-commentary -match "Gene-commentary_heading:GeneOntology" -block Gene-commentary_comment -subset Gene-commentary -match "Dbtag_db:GO" -sep "|" -element Object-id_id Other-source_anchor > tempSearch;
#  if [ $(cat tempSearch | wc -l) -gt 0 ];then line=$(head -n1 tempSearch); echo -e "$gene\t$desc\t$line";else echo -e "$gene\t$desc";fi
#  echo "end of loop"
#done < $resources/temp

names=()
while read desc;do echo $desc; names+=("$desc");done < novels.protein_coding.pep.map.comp.simple2.id
for desc in "${names[@]:1:10}"; do 
  esearch -db gene -query "$desc AND human [ORGN]" -sort "Relevance" | efetch -format native -mode xml | xtract -pattern Entrezgene -element Gene-track_geneid Org-ref_taxname Gene-ref_locus Gene-ref_desc -group Gene-ref_db -match "Dbtag_db:HGNC" -first Object-id_str -group Gene-ref_syn -sep "|" -element Gene-ref_syn_E -group Gene-commentary -match "Gene-commentary_heading:GeneOntology" -block Gene-commentary_comment -subset Gene-commentary -match "Dbtag_db:GO" -sep "|" -element Object-id_id Other-source_anchor > tempSearch;
  if [ $(cat tempSearch | wc -l) -gt 0 ];then 
    cat tempSearch | awk -F '\t' '{print $4}' > tempSearch.IDs
    while read term;do python $script_path/distance2.py "$desc" "$term"; done < tempSearch.IDs > tempSearch.scores
    NUM=$(cat tempSearch.scores | sort -n | head -n1 | grep -n -w -f - tempSearch.scores | awk -F ':' '{print $1}')
    line=$(sed "${NUM}q;d" tempSearch); 
    echo -e "$desc\t$line";
  else echo -e "$desc";fi
done > ncbi.report.1.10

#  desc2=$(echo $desc | sed 's/^[[:space:]]//g' | sed 's/[[:space:]]$//g' | sed 's/ /[ALL] /g')
#  esearch -db gene -query "($desc2[ALL]) AND human [ORGN]" -sort "Relevance" | efetch -format native -mode xml | xtract -pattern Entrezgene -element Gene-track_geneid Org-ref_taxname Gene-ref_locus Gene-ref_desc -group Gene-ref_db -match "Dbtag_db:HGNC" -first Object-id_str -group Gene-ref_syn -sep "|" -element Gene-ref_syn_E -group Gene-commentary -match "Gene-commentary_heading:GeneOntology" -block Gene-commentary_comment -subset Gene-commentary -match "Dbtag_db:GO" -sep "|" -element Object-id_id Other-source_anchor > tempSearch;
#    cat tempSearch | awk -F '\t' '{print $4}' | xargs -I '{}' bash -c 'python $0/distance2.py "$desc" "{}"' "$script_path" > tempSearch.scores


## Blast2Go GUI
# 1. load blastp files *
# 2. Blast tool bar: Run blast description annotator 
# 3. load interproscan results
# 4. interproscan tool bar: Run Merge Interproscan GOs to Annotation
# 5. Run mapping
# 6. Run Annotation *
# 7. Annotation tool bar: Run Annex
# 8. Menu (Analysis) -> Enzyme code and KEGG -> Run GO-enzyme code ampping
# 9. Menu (Analysis) -> Enzyme Code and KEGG -> Load Pathway-Maps from KEGG *
# 10. save all maps
# 11. graphs
# 12. Menu (Analysis) -> GO-slim -> Run GO-slim
# 13. graphs

## GO Annotations: C-cellular component, F-molecular Function, P-biological process
cd GO_annotation
for dir in Current-*;do
  echo $dir
  Genespring=$(ls $dir/*.Genespring)
  GOStat=$(ls $dir/*.GOStat) 
  tail -n+2 $Genespring | sed 's/.$//' > $Genespring.2
  tail -n+2 $GOStat | sed 's/.$//' | sed 's/,/;/g' > $GOStat.2
  echo -e "Sequence Name\tGO Number\tBiological Process\tCellular Component\t Molecular Function" > $dir/GO.summary 
  paste $GOStat.2 $Genespring.2 | cut -f 1,2,4,5,6 >> $dir/GO.summary 
  echo "no of annotated genes"
  cat $dir/GO.summary | awk -F '\t' '{print $1}' | sort | uniq | wc -l;
  echo "no of unique GO annotation terms" 
  tail -n+2 $dir/GO.summary | awk -F '\t' '{print $2}' | tr ';' '\n' | sort | uniq | wc -l
  tail -n+2 $dir/GO.summary | awk -F '\t' '{print $2}' | tr ';' '\n' | sort | uniq > $dir/GO.models
  echo "Frequnecy of GO terms"
  tail -n+2 $dir/GO.summary | awk -F '\t' '{print $2}' | tr ';' '\n' | awk  -F '\t' '{A[$1]++}END{for(i in A)print A[i]"\t"i}' | sort -k1,1nr > $dir/GO.freq
  tail -n+2 $dir/GO.summary | awk -F '\t' '{if($3!="") print $3}' | tr ';' '\n' | awk  -F '\t' '{A[$1]++}END{for(i in A)print A[i]"\t"i}' | sort -k1,1nr > $dir/GO.BP.freq
  tail -n+2 $dir/GO.summary | awk -F '\t' '{if($4!="") print $4}' | tr ';' '\n' | awk  -F '\t' '{A[$1]++}END{for(i in A)print A[i]"\t"i}' | sort -k1,1nr > $dir/GO.CC.freq
  tail -n+2 $dir/GO.summary | awk -F '\t' '{if($5!="") print $5}' | tr ';' '\n' | awk  -F '\t' '{A[$1]++}END{for(i in A)print A[i]"\t"i}' | sort -k1,1nr > $dir/GO.MF.freq
done
## no of annotated genes for all genes 
cat */GO.summary | awk -F '\t' '{print $1}' | sort | uniq | wc -l
## no of unique GO annotation terms for all genes
cat */GO.models | sort | uniq | wc -l
cat */GO.models | sort | uniq > GO.models
## froquency of GO annotation terms for all genes
cat */GO.models | sort | uniq | wc -l
cat */GO.models | sort | uniq > GO.models

cd GOSlim_annotation
for dir in Current-*;do
  echo $dir
  Genespring=$(ls $dir/*.Genespring)
  cat $Genespring | sed 's/.$//' > $Genespring.2
  echo "Frequnecy of GO terms"
  tail -n+2 $Genespring.2 | awk -F '\t' '{if($2!="") print $2}' | tr ';' '\n' | awk  -F '\t' '{A[$1]++}END{for(i in A)print A[i]"\t"i}' | sort -k1,1nr > $dir/GO.BP.freq
  tail -n+2 $Genespring.2 | awk -F '\t' '{if($3!="") print $3}' | tr ';' '\n' | awk  -F '\t' '{A[$1]++}END{for(i in A)print A[i]"\t"i}' | sort -k1,1nr > $dir/GO.CC.freq
  tail -n+2 $Genespring.2 | awk -F '\t' '{if($4!="") print $4}' | tr ';' '\n' | awk  -F '\t' '{A[$1]++}END{for(i in A)print A[i]"\t"i}' | sort -k1,1nr > $dir/GO.MF.freq
done
cat */*.Genespring.2 | grep -v "^Sequence" > all.GenespringNoheader
cat all.GenespringNoheader | awk -F '\t' '{if($2!="") print $2}' | tr ';' '\n' | awk  -F '\t' '{A[$1]++}END{for(i in A)print A[i]"\t"i}' | sort -k1,1nr > GO.BP.freq
cat all.GenespringNoheader | awk -F '\t' '{if($3!="") print $3}' | tr ';' '\n' | awk  -F '\t' '{A[$1]++}END{for(i in A)print A[i]"\t"i}' | sort -k1,1nr > GO.CC.freq
cat all.GenespringNoheader | awk -F '\t' '{if($4!="") print $4}' | tr ';' '\n' | awk  -F '\t' '{A[$1]++}END{for(i in A)print A[i]"\t"i}' | sort -k1,1nr > GO.MF.freq

