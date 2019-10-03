#!/bin/bash

	# BASH Script for FILTERING SNVs and Indels

VCFLIST=$1
OUTDIR=$2

while read FILE; do 
            NAME=`basename $FILE`
	 	echo "$NAME"

# CHECK IF Filtered SNV FILE ALREADY EXISTS
	if [ ! -f $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".snvs.filtered.vcf ]; then 
			
			# Unzip File
				bgzip -d $FILE
			
				egrep "^#" $FILE >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".metaheader.vcf		
			
			# 2. Extract SNVS
					awk '{if(length($5)==1 && length($4)==1) print $0;}' $FILE  >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_default.snvs.temp.vcf

				# Extract Indels
					awk '{if(length($5)!=length($4)) print $0;}' $FILE >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".indels.vcf

	
	
			# 3. Filter out the badReads, Mapping Quality, and Sequence Complexity
					awk '!(($7 ~ /badReads/) || ($7 ~ /MQ/) || ($7 ~ /SC/))'  $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_default.snvs.temp.vcf > $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_default.snvs.temp.flags.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_default.snvs.temp.vcf
	
	
	
			# 4. Remove Chromosome Unknown and Mitochondria
					egrep -v "ChrU|MT" $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_default.snvs.temp.flags.vcf > $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_default.snvs.temp.flags.Chromosomes.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_default.snvs.temp.flags.vcf


	
			# 5. Merge Header and SNV VCF
					cat $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".metaheader.vcf $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_default.snvs.temp.flags.Chromosomes.vcf >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".snvs.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_default.snvs.temp.flags.Chromosomes.vcf

	
	
			# 6. Remove SNVs at the Ends of Each Supercontig
		 			vcfintersect -b $OUTDIR/3_SNV_filtering_individual/supercontigs.bed $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".snvs.vcf >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.snvs.filtered.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".snvs.vcf



			# 7. Remove SNVs at the Ends of Each Contigs
		 			vcfintersect -v  -b $OUTDIR/3_SNV_filtering_individual/contigGaps.bed $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.snvs.filtered.vcf >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.snvs.filtered1.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.snvs.filtered.vcf
	


			# 8. Remove Simple Repeats
		 			vcfintersect -v  -b $OUTDIR/3_SNV_filtering_individual/simpleRepeats.bed $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.snvs.filtered1.vcf >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.snvs.filtered2.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.snvs.filtered1.vcf



			# 9. Remove variants that are homozygous in 91H
					vcfintersect -v  -b $OUTDIR/3_SNV_filtering_individual/variants.91H.version2.bed $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.snvs.filtered2.vcf >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.snvs.filtered3.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.snvs.filtered2.vcf



			# 10. VCF-filter MMLQ
					vcffilter -f "MMLQ > 30" $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.snvs.filtered3.vcf >  $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".snvs.filtered.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.snvs.filtered3.vcf					
					
					
					
			# 11. Remove variants near indels
					awk '{if(length($5)!=1 || length($4)!=1) print $0;}' $FILE  >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".mnvs.temp.vcf
					awk '{if(length($5)==length($4)) print $0;}' $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".mnvs.temp.vcf  >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".mnvs.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".mnvs.temp.vcf

					egrep "^#" -v $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".indels.vcf | awk '{if ($2-21 < 0) start=0; else start=$2-21; if (length($4)>length($5)) ftprint=length($4)-1; else ftprint=0; printf("%s\t%s\t%s\n",$1,start,$2+21)}' > $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_indel.20.bed

					vcfintersect -v  -b $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_indel.20.bed $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".snvs.filtered.vcf >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".snvs.filtered.indels.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".metaheader.vcf


			# 12. Compress Files and Move them to Appropriate Folders
					bgzip $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".snvs.filtered.vcf
					bgzip $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".mnvs.vcf
					bgzip $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".indels.vcf

					bgzip $FILE

					mv $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".snvs.filtered.vcf.gz $OUTDIR/3_SNV_filtering_individual/Archive/
					mv $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".mnvs.vcf.gz $OUTDIR/3_SNV_filtering_individual/MNVs/
					mv $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".indels.vcf.gz $OUTDIR/3_SNV_filtering_individual/Indels/
					mv $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_indel.20.bed $OUTDIR/3_SNV_filtering_individual/Indels/
	else
			echo "$Name SNP exists\n"
	fi
		
		
		
# CHECK IF Filtered Indel FILE ALREADY EXISTS
	if [ ! -f $OUTDIR/3_SNV_filtering_individual/Indels/"${NAME%.*}".indels.filtered.vcf.gz ]; then 
			
				# Unzip File
				bgzip -d $FILE
			
				egrep "^#" $FILE >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".metaheader.vcf		
			
	
				# Extract Indels
					awk '{if(length($5)!=length($4)) print $0;}' $FILE >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_indels.vcf

			
			
			# 3. Filter out the badReads, Mapping Quality, and Sequence Complexity
					awk '!(($7 ~ /badReads/) || ($7 ~ /MQ/) || ($7 ~ /SC/))'  $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_indels.vcf > $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_indels.flags.vcf	
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_indels.vcf
	
	
	
			# 4. Remove Chromosome Unknown and Mitochondria
					egrep -v "ChrU|MT" $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_indels.flags.vcf > $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_indels.flags.Chromosomes.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_indels.flags.vcf


	
			# 5. Merge Header and SNV VCF
					cat $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".metaheader.vcf $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_indels.flags.Chromosomes.vcf >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".indels.filtering.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_indels.flags.Chromosomes.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".metaheader.vcf



			# 6. Remove SNVs at the Ends of Each Supercontig
					vcfintersect -b $OUTDIR/3_SNV_filtering_individual/supercontigs.bed $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".indels.filtering.vcf >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.indels.filtering.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".indels.filtering.vcf



			# 7. Remove SNVs at the Ends of Each Contigs
					vcfintersect -v  -b $OUTDIR/3_SNV_filtering_individual/contigGaps.bed $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.indels.filtering.vcf >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.indels.filtering1.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.indels.filtering.vcf
		
	

			# 8. Remove Simple Repeats
		 			vcfintersect -v  -b $OUTDIR/3_SNV_filtering_individual/simpleRepeats.bed $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.indels.filtering1.vcf >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.indels.filtering2.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.indels.filtering1.vcf



			# 9. Remove variants that are homozygous in 91H
					vcfintersect -v  -b $OUTDIR/3_SNV_filtering_individual/variants.91H.version2.bed $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.indels.filtering2.vcf >$OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.indels.filtering3.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.indels.filtering2.vcf



			# 10. VCF-filter MMLQ
					vcffilter -f "MMLQ > 30" $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.indels.filtering3.vcf >  $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".indels.filtered.vcf
					rm $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}"_multiallelic.indels.filtering3.vcf
					
					
			# Compress and Move Filtered Indels VCF
					bgzip $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".indels.filtered.vcf
					mv $OUTDIR/3_SNV_filtering_individual/"${NAME%.*}".indels.filtered.vcf.gz $OUTDIR/3_SNV_filtering_individual/Indels/

					

	else
			echo "$Name Indel exists\n"
	fi
		
					
		
done<$VCFLIST


