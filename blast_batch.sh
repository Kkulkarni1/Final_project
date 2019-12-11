#!/usr/bin/env ksh93

set -o errexit
set -o nounset

INDIR=arab_1000_split
OUTDIR=arab_1000_out
JOBMAX=100

for path in $INDIR/*.txt; do
    file=`basename $path .txt`;
    blastp -db prot2003-2014.fa -query $path -outfmt 6 -out $OUTDIR/$file-blast_out.txt & 
done
wait 

