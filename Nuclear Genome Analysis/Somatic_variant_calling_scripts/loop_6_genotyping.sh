#!/bin/bash

REFERENCE=$1
BAMLIST=$2
VCFREGIONS=$3
#VCFLIST=$4
OUTDIR=$4 # This is the directory of 6_genotyping


while read -r VCF REGIONS; do
	NAME=`basename $VCF`
	GTP=1
	echo "$NAME"
#if [ -s $OUTDIR/"${NAME%.*}"_genotyped_${GTP}.vcf ]; then
        cut -f1-5 $VCF > $OUTDIR/"${NAME%.*}"_merged_pos.txt
	bgzip $VCF
	tabix -p vcf $VCF.gz


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
		--bufferSize=50000 \
		--regions=$REGIONS \
		--source=$VCF.gz \
        --nCPU=1 \
        -o $OUTDIR/"${NAME%.*}"_genotyped_${GTP}.vcf

	vcfbreakmulti $OUTDIR/"${NAME%.*}"_genotyped_${GTP}.vcf >$OUTDIR/"${NAME%.*}"_default.multiallelic.vcf

#	fi
	# start re-iterative genotyping
	DONE="false"

while [[ "$DONE" != "true" && $GTP < "10" ]]; do
	
	#rm $OUTDIR/"${NAME%.*}"_genotyped.vcf
	#rm $OUTDIR/"${NAME%.*}"_coords.vcf*

    # Create a new regions file containing only the bases of the missing variants
       # The size of the region is 200bp for the first genotyping
        GTP=$((GTP +1))


	tail -n +49 $OUTDIR/"${NAME%.*}"_default.multiallelic.vcf | cut -f1-5 > $OUTDIR/"${NAME%.*}"_genotyped.vcf
	grep -vxFf $OUTDIR/"${NAME%.*}"_genotyped.vcf $OUTDIR/"${NAME%.*}"_merged_pos.txt > $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf

	vcf-sort $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf >$OUTDIR/"${NAME%.*}"_coords_${GTP}_sorted.vcf
	rm $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf
	mv $OUTDIR/"${NAME%.*}"_coords_${GTP}_sorted.vcf $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf
	cp $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf $OUTDIR/"${NAME%.*}"_coords_${GTP}_backup.vcf

	if [ $GTP=="2" ]; then
		awk '{if (length($4) >= length($5)) { print $1 ":" $2-50 "-" $2+length($4)+50 } else { print $1 ":" $2-50 "-" $2+length($5)+50 }}' $OUTDIR/"${NAME%.*}"_coords_${GT
P}.vcf > $OUTDIR/"${NAME%.*}"_regions_${GTP}.txt
	
	else
		awk '{if (length($4) >= length($5)) { print $1 ":" $2-1 "-" $2+length($4)+1 } else { print $1 ":" $2-1 "-" $2+length($5)+1 }}' $OUTDIR/"${NAME%.*}"_coords_${GTP}.v
cf > $OUTDIR/"${NAME%.*}"_regions_${GTP}.txt
	
	fi
	
	Somatypus_MergeRegions.py $OUTDIR/"${NAME%.*}"_regions_${GTP}.txt
	
	#line='wc -l < $OUTDIR/"${NAME%.*}"_regions_${GTP}.txt' 
	
	if [[ $(wc -l <$OUTDIR/"${NAME%.*}"_regions_${GTP}.txt) -ge 1 ]]; then 
	echo -e "\nGenotyping missing calls"
	
	# Prepare VCF of missing variants
		bgzip $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf
		tabix -p vcf $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf.gz
	
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
	    	--bufferSize=50000 \
			--regions=$OUTDIR/"${NAME%.*}"_regions_${GTP}_merged.txt \
	    	--source=$OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf.gz \
        	--nCPU=1 \
        	-o $OUTDIR/"${NAME%.*}"_genotyped_${GTP}.vcf
	
	
	vcfbreakmulti $OUTDIR/"${NAME%.*}"_genotyped_${GTP}.vcf >$OUTDIR/"${NAME%.*}"_default.multiallelic_${GTP}.vcf
	
	tail -n +49 $OUTDIR/"${NAME%.*}"_default.multiallelic_${GTP}.vcf >> $OUTDIR/"${NAME%.*}"_default.multiallelic.vcf
	
	   if [ $GTP == "10" ]; then
	        GTP=$((GTP +1))


		 tail -n +49 $OUTDIR/"${NAME%.*}"_default.multiallelic.vcf | cut -f1-5 > $OUTDIR/"${NAME%.*}"_genotyped.vcf
        	 grep -vxFf $OUTDIR/"${NAME%.*}"_genotyped.vcf $OUTDIR/"${NAME%.*}"_merged_pos.txt > $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf

	         vcf-sort $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf >$OUTDIR/"${NAME%.*}"_coords_${GTP}_sorted.vcf
        	 rm $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf
        	 mv $OUTDIR/"${NAME%.*}"_coords_${GTP}_sorted.vcf $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf
        	 cp $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf $OUTDIR/"${NAME%.*}"_coords_${GTP}_backup.vcf
   		 awk '{if (length($4) >= length($5)) { print $1 ":" $2 "-" $2+length($4)-1 } else { print $1 ":" $2 "-" $2+length($5)-1 }}' $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf > $OUTDIR/"${NAME%.*}"_regions_${GTP}.txt 

	  

		# Prepare VCF of missing variants
                 bgzip $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf
                 tabix -p vcf $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf.gz
        
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
                		 --bufferSize=50000 \
                		 --regions=$OUTDIR/"${NAME%.*}"_regions_${GTP}.txt \
                		 --source=$OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf.gz \
                		 --nCPU=1 \
                		 -o $OUTDIR/"${NAME%.*}"_genotyped_${GTP}.vcf
        
        
        	vcfbreakmulti $OUTDIR/"${NAME%.*}"_genotyped_${GTP}.vcf >$OUTDIR/"${NAME%.*}"_default.multiallelic_${GTP}.vcf
        
        	tail -n +49 $OUTDIR/"${NAME%.*}"_default.multiallelic_${GTP}.vcf >> $OUTDIR/"${NAME%.*}"_default.multiallelic.vcf
        
	   fi
	
	else

		echo -e "\nNo missing calls"
		
		DONE="true"
	fi
	done

		cp $OUTDIR/"${NAME%.*}"_default.multiallelic.vcf $OUTDIR/genotyped_"${NAME%.*}".multiallelic.vcf
		cp $VCF.gz $OUTDIR/prior_genotyping_"${NAME%.*}".multiallelic.vcf
 		
		tail -n +49 $OUTDIR/"${NAME%.*}"_default.multiallelic.vcf | cut -f1-5 > $OUTDIR/"${NAME%.*}"_genotyped.vcf
                grep -vxFf $OUTDIR/"${NAME%.*}"_genotyped.vcf $OUTDIR/"${NAME%.*}"_merged_pos.txt > $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf

                vcf-sort $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf >$OUTDIR/"${NAME%.*}"_coords_${GTP}_sorted.vcf
                rm $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf
                mv $OUTDIR/"${NAME%.*}"_coords_${GTP}_sorted.vcf $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf
                cp $OUTDIR/"${NAME%.*}"_coords_${GTP}.vcf $OUTDIR/"${NAME%.*}"_coords_${GTP}_backup.vcf
         
		cp $OUTDIR/"${NAME%.*}"_coords_${GTP}_backup.vcf $OUTDIR/variants_remaining_"${NAME%.*}".vcf
		#rm $OUTDIR/"${NAME%.*}"*

done<$VCFREGIONS



