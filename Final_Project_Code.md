#Commands used to process and generate the files

## Downloading Rice protein sequences

    wget "http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.pep"

## Downloading Arabidopsis protein sequences

    wget "ftp://ftp.arabidopsis.org/home/tair/Proteins/TAIR10_protein_lists/TAIR10_pep_20101214"

## Downloading Bacteria protein sequences

`wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria"`

## Arabidopsis Domain AND GO information:
`wget "ftp://ftp.arabidopsis.org/home/tair/Proteins/Domains/TAIR10_all.domains"
wget "ftp://ftp.arabidopsis.org/home/tair/Ontologies/Gene_Ontology/ATH_GO_GOSLIM.txt"``
`
## Rice Domain information:

    `http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.pfam
`
## Processing files

## Making blast database with bacteria protein sequences

    makeblastdb –in prot2003-2014.fa –dbtype prot

## To make the processing faster, rice and protein sequences were divided into subset of 1000 sequences

`awk -v size=1000 -v pre=prefix -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' all.pep 
awk -v size=1000 -v pre=prefix -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' TAIR10_pep_20101214`

## Adding suffix '.txt' to all files to make processing easier

    
    for i in *; do mv $i ${i}.txt; done

## Running blast against protein sequences with arabidopsis and rice sequences as queries using ksh shell
    bash blast_batch.sh 

## 'bash blast_batch.sh' script:
#!/usr/bin/env ksh93

`set -o errexit
set -o nounset
INDIR=rice_seq_1000_split
OUTDIR=rice_seq_1000_out
JOBMAX=100`
for path in $INDIR/*.txt; do
    file=`basename $path .txt`;
    blastp -db prot2003-2014.fa -query $path -outfmt 6 -out $OUTDIR/$file-blast_out.txt &
done
wait

(simlarly done for arabidopsis)

## Processing arabidopsis BLAST output files
## Concatenating all output files

    cat *.txt > arab_prtn.txt

## Filtering blast hit file with evalue (1 X 10^-6)

    cat arab_prtn.txt | awk '{if($11<1e-6){print}}' > arab_prtn_filtered.txt

## Inspecting unique protein Ids
    awk '{print $1}' arab_prtn_filtered.txt | sort | uniq > arab_prtn_filtered_uniq.txt

## Extracting uniq protein Ids

    awk '{print $1}' arab_prtn_filtered_uniq.txt > arab_prtn_filtered_uniq_ID.txt

## Extracting chromosome information from protein Ids

    awk '{print $1, substr($1,3,1)}' arab_prtn_filtered_uniq_ID.txt > old_arab_prtn_chr_name.txt

## Subtracting old proteins IDs (processed blast output) from original prtein file to get new proteins IDs   

    grep -Fvf old_arab_prtn_chr_name_sorted.txt all_prtn_ID_sorted.txt > new_prtn_arab_ID.txt

## Processing old and new protein files to make plots in R

## Processing rice BLAST output files
## Concatenating all output files
    cat *.txt > rice_prtn.txt

## Filtering blast hit file with evalue (1 X 10^-6)
    cat rice_prtn.txt | awk '{if($11<1e-6){print}}' > rice_prtn_filtered.txt

## Inspecting unique protein Ids

    awk '{print $1}' rice_prtn_filtered.txt | sort | uniq > rice_prtn_filtered_uniq.txt

# Extracting uniq protein Ids
awk '{print $1}' rice_prtn_filtered_uniq.txt > rice_prtn_filtered_uniq_ID.txt

# Extracting chromosome information from protein Ids
awk '{print $1, substr($1,3,2)}' rice_prtn_filtered_uniq_ID.txt > old_rice_prtn_chr_name.txt

## Subtracting old proteins IDs (processed blast output) from original prtein file to get new proteins IDs
grep -Fvf old_rice_prtn_chr_name_sorted.txt all_prtn_ID_sorted.txt > new_prtn_rice_ID.txt

## Processing old and new protein files to make plots in R

##Processing Arabidopsis files
`awk 'BEGIN { OFS = " " } {print $0, "old"}' total_filtered_chr_name.txt > old_proteins_arabidopsis.txt
head old_proteins_arabidopsis.txt
tail old_proteins_arabidopsis.txt
head new_prtn_ID_chr_name.txt
awk 'BEGIN { OFS = " " } {print $0, "new"}' new_prtn_ID_chr_name.txt > new_proteins_arabidopsis.txt
head new_proteins_arabidopsis.txt
tail new_proteins_arabidopsis.txt
cat old_proteins_arabidopsis.txt new_proteins_arabidopsis.txt > total_proteins_arabidopsis.txt
head tail total_proteins_arabidopsis.txt`

## Inserting header for total_proteins_arabidopsis.txt - 
`vi total_proteins_arabidopsis.txt
I - Insert
Gene Chromosome Protein_Type
Escape
:wq!``

##Processing Rice files
`
awk 'BEGIN { OFS = " " } {print $0, "old"}' old_rice_prtn_chr_name.txt > old_proteins_rice.txt
head old_proteins_rice.txt
tail old_proteins_rice.txt
awk 'BEGIN { OFS = " " } {print $0, "new"}' new_prtn_rice_ID.txt > new_proteins_rice.txt
head new_proteins_rice.txt
tail new_proteins_rice.txt
cat old_proteins_rice.txt new_proteins_rice.txt > total_proteins_rice.txt
head tail total_proteins_rice.txt`
## Inserting header for total_proteins_rice.txt - 
vi total_proteins_rice.txt
I ## Insert
Gene Chromosome Protein_Type
Escape
:wq!



## Processing Domain and Go files for arabidopsis to make plots in R
`cut -f3- -d $'\t' ATH_GO_GOSLIM.txt ATH_GO_GOSLIM_V2.txt
sort -k1,1 ATH_GO_GOSLIM_V2.txt > ATH_GO_GOSLIM_V2_sorted.txt
join -t $'\t' -1 1 -2 1 ATH_GO_GOSLIM_V2_sorted.txt old_arab_prtn_sorted.txt > old_arab_GO.txt
sort -k1,1 TAIR10_all.domains > TAIR10_all.domains_sorted.txt
join -t $'\t' -1 1 -2 1 TAIR10_all.domains_sorted.txt old_arab_prtn_sorted.txt > old_arab_domain.txt `

## Processing Domain and Go files for rice to make plots in R

`sort -k1,1 all.GOSlim_assignment > all.GOSlim_assignment_sorted.txt
sort -k1,1 all.interpro > all.interpro_sorted.txt
sort -k1,1 all.pfam > all.pfam_sorted.txt
sort -k1,1 old_rice_prtn.txt > old_rice_prtn_sorted.txt
join -t $'\t' -1 1 -2 1  all.GOSlim_assignment_sorted.txt old_rice_prtn_sorted.txt > old_rice_GO_V1.txt
join -t $'\t' -1 1 -2 1  all.interpro_sorted.txt old_rice_prtn_sorted.txt > old_rice_interpro_V1.txt
join -t $'\t' -1 1 -2 1  all.pfam_sorted.txt old_rice_prtn_sorted.txt > old_rice_pfam_V1.txt`






