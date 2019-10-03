#!/bin/bash

	# BASH Script for Separating VCF into Multiallelic

LIST=$1 # LIST of Initial Run PLATYPUS Calling
OUTDIR=$2

while read FILE; do 
            NAME=`basename $FILE`
		if [ -f $OUTDIR/2_filter_individual/"${NAME%.*}".multiallelic.vcf ]; then 
	   		echo "$NAME  exists\n"
         else
        	echo "$NAME"
        	
	# 1. Break any Multi-allelic Variants into SNVs and merge it into the existing VCF
				vcfbreakmulti $FILE  >$OUTDIR/2_filter_individual/"${NAME%.*}".multiallelic.vcf
				bgzip $FILE
		fi
done<$LIST
