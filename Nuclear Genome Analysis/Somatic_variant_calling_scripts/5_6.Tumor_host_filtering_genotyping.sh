#!/bin/bash

GENOTYPEDVCF=$1
OUTDIR=$2


	# For Each Tumor, remove all variants identified in the host using R:
		mkdir -p $OUTDIR/5.Prep.Genotyping/

		while read FILE; do 
	
			# Remove any Variant observed in the Host, Quality Score < 10, and MQ, SC, badReads Flags
				Rscript 5.Variant_Filterings.R $FILE $OUTDIR/5.Prep.Genotyping/
	
		done<$GENOTYPEDVCF
	
	
	
	# Genotyping Across All Genotyped VCF Files
		mkdir -p $OUTDIR/6.Genotyping.tumors/
		mkdir -p $OUTDIR/6.Genotyping.hosts/
		
		cat $(ls -t $OUTDIR/5.Prep.Genotyping/*.vcf) > $OUTDIR/5.Prep.Genotyping/final.genotypinglist.vcf
		vcf-sort $OUTDIR/5.Prep.Genotyping/final.genotypinglist.vcf | uniq > $OUTDIR/5.Prep.Genotyping/final.genotypinglist.unique.sorted.vcf
		
		
		# Split Genotyping VCF File into 20 chunks
		a=(`wc -l $OUTDIR/5.Prep.Genotyping/final.genotypinglist.unique.sorted.vcf`) ; lines=`echo $a/20 | bc -l` ; split -l=$lines -d  $OUTDIR/5.Prep.Genotyping/final.genotypinglist.unique.sorted.vcf split.genotypinglist.unique.sorted


		
		# Preparation for Genotyping List
			bgzip $OUTDIR/5.Prep.Genotyping/final.genotypinglist.unique.sorted.vcf
			tabix -p vcf $OUTDIR/5.Prep.Genotyping/final.genotypinglist.unique.sorted.vcf.gz
		
		
		
		
		# Genotyping Across All Tumors:
			Platypus.py callVariants \
		        	--logFileName=$OUTDIR/"${NAME%.*}"_default.log \
    		                --refFile=$REFERENCE \
         		        --bamFiles=$OUTDIR/4_Tumor_Host_Filtering/bamlist.txt \
                        	--minReads=1 \
		                --minFlank=0 \
		                --minBaseQual=30 \
      			        --filteredReadsFrac=.8 \
		                --badReadsThreshold=30 \
	                        --badReadsWindow=15 \
                                --minPosterior=0 \
	                        --getVariantsFromBAMs=0 \
				--regions=$OUTDIR/4_Tumor_Host_Filtering/${NAME%.*}_varRegions_snvs.txt \
	                        --source=$OUTDIR/5.Prep.Genotyping/final.genotypinglist.sorted.vcf.gz \
                                --nCPU=8 \
                                -o $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped4.vcf

		
		# Genotyping Across All Hosts:
					
		
		
		
		
		
		
	
	
		

	
	

