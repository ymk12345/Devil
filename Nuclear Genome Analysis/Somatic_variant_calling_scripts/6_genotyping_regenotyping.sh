#!/bin/bash

REFERENCE=$1
BAMLIST=$2
REGIONS=$3
VCF=$4
OUTDIR=$5 # This is the directory of 6_genotyping

	NAME=`basename $VAR`
	GTP=1

 	Platypus.py callVariants \
	    --logFileName="${NAME%.*}"_genotyped_${GTP}_default.log \
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
        -o $OUTDIR/"${NAME%.*}"_genotyped_${GTP}.vcf

	bgzip -d $VAR
	cut -f1-5 $VAR > $OUTDIR/"${NAME%.*}"_merged_pos.txt
	vcfbreakmulti $OUTDIR/"${NAME%.*}"_genotyped_${GTP}.vcf >$OUTDIR/"${NAME%.*}"_default.multiallelic.vcf


	# start re-iterative genotyping
	DONE=false

	while [ "$DONE" != true ]; do
	
	rm $OUTDIR/"${NAME%.*}"_genotyped.vcf
	rm $OUTDIR/"${NAME%.*}"_coords.vcf
	
	tail -n +49 $OUTDIR/"${NAME%.*}"_default.multiallelic.vcf | cut -f1-5 > $OUTDIR/"${NAME%.*}"_genotyped.vcf
	grep -vxFf $OUTDIR/"${NAME%.*}"_genotyped.vcf $OUTDIR/"${NAME%.*}"_merged_pos.txt > $OUTDIR/"${NAME%.*}"_coords.vcf

	vcf-sort $OUTDIR/"${NAME%.*}"_coords.vcf
	
    # Create a new regions file containing only the bases of the missing variants
	    # The size of the region is 200bp for the first genotyping
		GTP=$(( GTP + 1 ))

	if[ "$GTP" == "2" ]; then
		awk '{if (length($4) >= length($5)) { print $1 ":" $2-50 "-" $2+length($4)+50 } else { print $1 ":" $2-50 "-" $2+length($5)+50 }}' $OUTDIR/"${NAME%.*}"_coords.vcf > $OUTDIR/"${NAME%.*}"_regions_${GTP}.txt
	
	else
		awk '{if (length($4) >= length($5)) { print $1 ":" $2-1 "-" $2+length($4)+1 } else { print $1 ":" $2-1 "-" $2+length($5)+1 }}' $OUTDIR/"${NAME%.*}"_coords.vcf > $OUTDIR/"${NAME%.*}"_regions_${GTP}.txt
	
	fi
	
	Somatypus_MergeRegions.py $OUTDIR/"${NAME%.*}"_regions_${GTP}.txt
	
	
	if [ -s $OUTDIR/"${NAME%.*}"_regions_${GTP}.txt ]; then 
	echo -e "\nGenotyping missing calls"
	
	# Prepare VCF of missing variants
		bgzip $OUTDIR/"${NAME%.*}"_coords.vcf
		tabix -p vcf $OUTDIR/"${NAME%.*}"_coords.vcf.gz
	
	# Run Platypus
	 	Platypus.py callVariants \
	    	--logFileName="${NAME%.*}"_genotyped_${GTP}_default.log \
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
			--regions=$OUTDIR/"${NAME%.*}"_regions_${GTP}_merged.txt \
	    	--source=$OUTDIR/"${NAME%.*}"_coords.vcf.gz \
        	--nCPU=1 \
        	-o $OUTDIR/"${NAME%.*}"_genotyped_${GTP}.vcf
	
	
	vcfbreakmulti $OUTDIR/"${NAME%.*}"_genotyped_${GTP}.vcf >$OUTDIR/"${NAME%.*}"_default.multiallelic_${GTP}.vcf
	
	tail -n +49 $OUTDIR/"${NAME%.*}"_default.multiallelic_${GTP}.vcf >> $OUTDIR/"${NAME%.*}"_default.multiallelic.vcf
	
	else
		echo -e "\nNo missing calls"
		
		DONE=true

	done

	



