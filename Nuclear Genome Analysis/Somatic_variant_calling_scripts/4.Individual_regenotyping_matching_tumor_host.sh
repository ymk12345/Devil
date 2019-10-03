#!/bin/bash

# BASH Script for 1. Matching Tumors and Host and 2. Re-iterative genotyping for SNV or INDEL Variant Calling


IDS=$1 # IDs where there is a matching Host for the Tumors (Use the Matching.ID.R Script)
BAMSDIR=$2 # BAM FILE Directory
REFERENCE=$3 # Reference
OUTDIR=$4 # The Directory where all the analyses is held


while read NAME; do
		
		# Remove Existing lists that we use for running Platypus on the matching tumor/host sets
				rm $OUTDIR/4_Tumor_Host_Filtering/bamlist.txt
				rm $OUTDIR/4_Tumor_Host_Filtering/temp.vcf.list
			
		
		# Identify the filtered VCFs for matching tumor(s) and Host
				ls -1 $OUTDIR/3_SNV_filtering_individual/platypusVariants_${NAME%.*}[A-Z]*.indels.vcf > $OUTDIR/4_Tumor_Host_Filtering/temp.vcf.list
		
			
		# Identify the Appropriate BAM Files that will be used for Platypus
				ls -1 $BAMSDIR/${NAME%.*}[A-Z]*.bam > $OUTDIR/4_Tumor_Host_Filtering/bamlist.txt


		# Merge the Tumor(s) and Host Filtered SNV VCFs together
				while read FILE; do 
			
					egrep -v "^#" $FILE |cut -f1-5 >> $OUTDIR/4_Tumor_Host_Filtering/temp.vcf
				
				done<$OUTDIR/4_Tumor_Host_Filtering/temp.vcf.list


	# Prepare the VCF File for Genotyping
		
				# Sort Merged VCF and Extract only the Unique Values		
						vcf-sort $OUTDIR/4_Tumor_Host_Filtering/temp.vcf | uniq > $OUTDIR/4_Tumor_Host_Filtering/${NAME%.*}_to.genotype.vcf
						rm $OUTDIR/4_Tumor_Host_Filtering/temp.vcf
			
				# Save a copy
						cp $OUTDIR/4_Tumor_Host_Filtering/${NAME%.*}_to.genotype.vcf $OUTDIR/4_Tumor_Host_Filtering/${NAME%.*}_to.genotype.backup.vcf

				# Compress and Index VCF File 
                        bgzip $OUTDIR/4_Tumor_Host_Filtering/${NAME%.*}_to.genotype.vcf
                        tabix -p vcf $OUTDIR/4_Tumor_Host_Filtering/${NAME%.*}_to.genotype.vcf.gz


	# First Genotyping Step:
				# Call Using Platypus (don't use regions here, actually takes longer)		
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
	   			 				--source=$OUTDIR/4_Tumor_Host_Filtering/${NAME%.*}_to.genotype.vcf.gz \
         			   			--nCPU=8 \
            					-o $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped.vcf


		## As certain values are multi-allelic, split the VCF into one where multi-allelic is in each row
			vcfbreakmulti $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped.vcf > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped.multiallelic.vcf



	# Identify Missing Calls
		 			# Create VCF of Missing Variants
							tail -n +49 $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped.multiallelic.vcf | cut -f1-5 > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_geno_pos.txt
							cut -f1-5 $OUTDIR/4_Tumor_Host_Filtering/${NAME%.*}_to.genotype.backup.vcf > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_merged_pos.txt
		        
		        			grep -vxFf $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_geno_pos.txt $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_merged_pos.txt > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.vcf
							rm $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_merged_pos.txt
			



# Run Second Genotype if all variants are not called in the first genotype:
		if [ -s $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.vcf ]; then 
			
					# Create Coordinates Positions to call Variants using 200bp Window and Merge
							awk '{if (length($4) >= length($5)) { print $1 ":" $2-200 "-" $2+length($4)+200 } else { print $1 ":" $2-200 "-" $2+length($5)+200 }}' $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.vcf > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_varRegions_snvs.txt
							Somatypus_MergeRegions.py $OUTDIR/4_Tumor_Host_Filtering/${NAME%.*}_varRegions_snvs.txt

					# Compress and Index VCF File of Missing Variants 
							vcf-sort $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.vcf > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.sorted.vcf
							bgzip $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.sorted.vcf
							tabix -p vcf $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.sorted.vcf.gz




	# Second Genotyping Step:
					# Call Using Platypus: Use Regions here	(200bp window)

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
                                --regions=$OUTDIR/4_Tumor_Host_Filtering/${NAME%.*}_varRegions_snvs_merged.txt \
                                --source=$OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.sorted.vcf.gz \
                                --nCPU=8 \
                                -o $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped2.vcf
			
			
			## As certain values are multi-allelic, split the VCF into one where multi-allelic is in each row
			 	vcfbreakmulti $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped2.vcf > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped2.multiallelic.vcf


			## Merge the output of the second genotype to the first genotyped vcf file (for both unsplit and split)
				egrep -v "^#" $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped2.vcf >>$OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped.vcf
				egrep -v "^#" $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped2.multiallelic.vcf >>$OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped.multiallelic.vcf




	# Identify Missing Calls from Second Genotyping 
		 			# Create VCF of Missing Variants
		 					rm $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_geno_pos.txt
	 						tail -n +49 $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped2.multiallelic.vcf | cut -f1-5 > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_geno_pos.txt
                        	grep -vxFf $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_geno_pos.txt $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.vcf > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords2.vcf
 
 
 
 # Run Third Genotype if all variants are not called in the second genotype:
 		if [ -s $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords2.vcf ]; then 
 
 
 					# Create Coordinates Positions to call Variants using 50bp Window and Merge
 							rm $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_varRegions_snvs.txt
 							rm $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_varRegions_snvs_merged.txt
 							
							awk '{if (length($4) >= length($5)) { print $1 ":" $2-50 "-" $2+length($4)+50 } else { print $1 ":" $2-50 "-" $2+length($5)+50 }}' $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords2.vcf > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_varRegions_snvs.txt
							Somatypus_MergeRegions.py $OUTDIR/4_Tumor_Host_Filtering/${NAME%.*}_varRegions_snvs.txt

 
					# Compress and Index VCF File of Missing Variants                         
                        		vcf-sort $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords2.vcf > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.sorted2.vcf
                        		bgzip $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.sorted2.vcf
                        		tabix -p vcf $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.sorted2.vcf.gz




		# Third Genotyping Step:
						# Call Using Platypus: Use Regions here	(50bp window)

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
                	                --source=$OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.sorted2.vcf.gz \
                    	            --nCPU=8 \
                        	        -o $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped3.vcf


				## As certain values are multi-allelic, split the VCF into one where multi-allelic is in each row
        	            vcfbreakmulti $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped3.vcf > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped3.multiallelic.vcf


				## Merge the output of Third Genotyping to the First genotyped vcf file
    	                    egrep -v "^#" $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped3.vcf >>$OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped.vcf
        	                egrep -v "^#" $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped3.multiallelic.vcf >>$OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped.multiallelic.vcf




		# Identify Missing Calls from Third Genotyping 
			 			# Create VCF of Missing Variants
			 					rm $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_geno_pos.txt
	 							tail -n +49 $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped3.multiallelic.vcf | cut -f1-5 > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_geno_pos.txt
                	        	grep -vxFf $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_geno_pos.txt $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.vcf > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords3.vcf
 
 
	  	# Run Fourth Genotype if all variants are not called in the third genotype:
 				if [ -s $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords3.vcf ]; then 
 
 
 
 							# Create Coordinates Positions to call Variants using 1bp Window and Merge
 									rm $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_varRegions_snvs.txt
 									rm $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_varRegions_snvs_merged.txt
 							
									awk '{if (length($4) >= length($5)) { print $1 ":" $2-1 "-" $2+length($4)+1 } else { print $1 ":" $2-1 "-" $2+length($5)+1 }}' $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords3.vcf > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_varRegions_snvs.txt
									Somatypus_MergeRegions.py $OUTDIR/4_Tumor_Host_Filtering/${NAME%.*}_varRegions_snvs.txt

 
							# Compress and Index VCF File of Missing Variants                         
            		            		vcf-sort $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords3.vcf > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.sorted3.vcf
                    		    		bgzip $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.sorted3.vcf
                        				tabix -p vcf $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.sorted3.vcf.gz



			# Fourth Genotyping Step:
							# Call Using Platypus: Use Regions here	(50bp window)

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
		                                --source=$OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_coords.sorted3.vcf.gz \
		                                --nCPU=8 \
 		                                -o $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped4.vcf


					## As certain values are multi-allelic, split the VCF into one where multi-allelic is in each row
        		            vcfbreakmulti $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped4.vcf > $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped4.multiallelic.vcf


					## Merge the output of Third Genotyping to the First genotyped vcf file
		                        egrep -v "^#" $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped4.vcf >>$OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped.vcf
        		                egrep -v "^#" $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped4.multiallelic.vcf >>$OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped.multiallelic.vcf

				fi
			fi
		fi

			## Rename FINAL Output
						mv $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped.vcf $OUTDIR/4_Tumor_Host_Filtering/final_regenotyped_"${NAME%.*}"_genotyped.vcf
						mv $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"_genotyped.multiallelic.vcf $OUTDIR/4_Tumor_Host_Filtering/final_regenotyped_"${NAME%.*}"_genotyped.multiallelic.vcf



			# Decided NOT to do past 4th Genotyping as if it didn't pick it up by this step, then these are likely not worth pursuing
						rm $OUTDIR/4_Tumor_Host_Filtering/"${NAME%.*}"*

	done<$IDS
