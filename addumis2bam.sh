#!/bin/sh
# -*- mode: sh; -*-

## This scripts finds the UMIs in the fastq files and adds them to the bam file,
## and calling the output *.withumis.bam.
## UMI's are the first 8 characters in both R1 and R2. They are added as optional tags
## at the end of each bam line using the tags (e.g.)  u1:Z:TGTGTGTG  and u2:Z:CTGAGTGT
## (both the R1 and R2 get the same u1 and u2 tags).
## For speed I'm using the venerable Unix join and paste, which requires that 
## reads are properly sorted in strict ASCII order. samtools view -n does *not* do this; 
## use sambamba view or picard tools, or even Unix sort (make sure you set LC_ALL=C for that)
##
## bam files are expected in $bamdir, fastq files in $fastqdir
## and things are written to $outdir, all defined below
##
## When submitting to the queue, be sure to request tmpspace (say
## twice the size of sam file, 10x that of the bam file). It is needed for the
## intermediates files and for sorting. The script did ~ 800k reads per
## minute.

set -x
set -u
set -e

fastqdir=$2
bamdir=$3
outdir=$4
## outdir=.

name="$1"

inbam=$bamdir/$name.namesorted.bam
fq1=$(ls $fastqdir/$name*R1*fastq.gz)
fq2=$(ls $fastqdir/$name*R2*fastq.gz)

ls -l $inbam $fq1 $fq2

if [ -f $inbam -a -f $fq1 -a $fq2 ]; then
: # aok
else
    echo "Not all of files: $inbam $fq1 $fq2 found, exiting" >&2
    exit 2
fi

umis1=$TMPDIR/$name-R1.umis
umis2=$TMPDIR/$name-R2.umis
umisboth=$TMPDIR/$name-both.umis
sam=$TMPDIR/$name.namesorted.sam

zcat $fq1 | paste -d '	' - - - - | awk '{print substr($1,2) "\t" "u1:Z:" substr($3, 1, 8)}'  > $umis1
zcat $fq2 | paste -d '	' - - - - | awk '{print                   "u2:Z:" substr($3, 1, 8)}'  > $umis2
paste -d '	' $umis1 $umis2  | env LC_ALL=C  nice sort -T $TMPDIR > $umisboth

samtools view $inbam > $sam

( samtools view -H $inbam
  echo -e "@PG\tID:addumis2bam.sh\tPG:addumis2bam.sh\tVN:0.0\tCL:$0 $@"
  join  -t '	' $sam $umisboth
)  | sambamba view -S -f bam /dev/stdin > $outdir/$name.withumis.bam

rm $umis1 $umis2 $umisboth $sam


