#!/bin/bash

REFERENCE=$1
BAMLIST=$2
REGIONS=$3
VCF=$4
OUTDIR=$5

NAME=`basename $VCF`

 	Platypus.py callVariants \
	       	--logFileName="${NAME%.*}"_default.log \
    		--refFile=$REFERENCE \
         	--bamFiles=$BAMLIST \
            	--minReads=1 \
		--minFlank=0 \
		--minBaseQual=30 \
      		--filteredReadsFrac=.8 \
		--badReadsThreshold=30 \
	        --badReadsWindow=15 \
            	--minPosterior=0 \
	        --getVariantsFromBAMs=0 \
	        --bufferSize=500 \
		--regions=$REGIONS \
	        --source=$VCF \
            	--nCPU=1 \
            	-o $OUTDIR/"${NAME%.*}".genotyped.vcf
                    
