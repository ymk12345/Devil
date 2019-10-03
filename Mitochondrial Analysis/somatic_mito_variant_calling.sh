#!/bin/bash

## Platypus SNP calling pipeline adapted for MT genomes

# INPUT
# $1: path to folder containing the input BAM files (accompanied by BAI indices)
# $2: path to genome FASTA (must have FAI index)
# $3: path to file of regions in chr:start-end format
# $4: path to Platypus.py
# $5: path to output folder (will be created)


if [ "$#" -ne 5 ]; then
    echo
    echo "| platypus_snp_pipeline_mt.sh: Calls SNPs from a set of BAM files using Platypus."
    echo "|"
    echo "| IMPORTANT: This script should be run on the Sanger Farm, with these options:"
    echo "|            (Memory amount will depend on the number of samples)"
    echo "|            -q hugemem -n 8 -R 'span[hosts=1] select[mem>=40000] rusage[mem=40000]' -M40000"
    echo "|"
    echo "| Input: Path to folder containing the input BAM files (accompanied by BAI indices)."
    echo "|        Path to reference genome FASTA file (must be accompanied by FAI index)."
    echo "|        Path to text file of regions to use, in CHR:START-END format, one per line."
    echo "|        Path to your Platypus installation (Platypus.py file)."
    echo "|        Path to the output folder (it will be created if needed)."
    echo "|"
    echo "| Usage: platypus_snp_pipeline_mt.sh /path/to/bams_folder /path/to/genome.fna /path/to/regions.txt /path/to/Platypus.py /path/to/out_folder"
    echo
    exit
fi


BAMSDIR="$1"
REFERENCE="$2"
REGIONS="$3"
PLATYPUS="$4"
OUTDIR="$5"
SCRIPTS="/nfs/dog_n_devil/adrian/software/somatypus/src"

function runPlatypus {
# 1. RUN PLATYPUS (DEFAULT) INDIVIDUALLY ON EVERY SAMPLE
echo -e "\n(1) RUNNING PLATYPUS (DEFAULT) INDIVIDUALLY ON EVERY SAMPLE\n"

mkdir -p $OUTDIR/1_individual_calls
mkdir $OUTDIR/logs

for FILE in `ls $BAMSDIR/*.bam`; do 
    NAME=`basename $FILE`
    echo -e "\nCalling on $NAME\n"
    $PLATYPUS callVariants \
    --logFileName=$OUTDIR/logs/1_individual_default.log \
    --refFile=$REFERENCE \
    --bamFiles=$FILE \
    --regions=$REGIONS \
    --minPosterior=0 \
    --minReads=3 \
    --nCPU=8 \
    -o $OUTDIR/1_individual_calls/platypusVariants_"${NAME%.*}"_default.vcf
done
}

function runPlatypusMinFlank {
# 2. RUN PLATYPUS (MINFLANK=0) INDIVIDUALLY ON EVERY SAMPLE
echo -e "\n(2) RUNNING PLATYPUS (MINFLANK=0) INDIVIDUALLY ON EVERY SAMPLE\n"

for FILE in `ls $BAMSDIR/*.bam`; do 
    NAME=`basename $FILE`
    echo -e "\nCalling on $NAME\n"
    $PLATYPUS callVariants \
    --logFileName=$OUTDIR/logs/2_individual_minFlank.log \
    --refFile=$REFERENCE \
    --bamFiles=$FILE \
    --regions=$REGIONS \
    --minPosterior=0 \
    --minFlank=0 \
    --trimReadFlank=10 \
    --minReads=3 \
    --nCPU=8 \
    -o $OUTDIR/1_individual_calls/platypusVariants_"${NAME%.*}"_minFlank.vcf
done
}

function splitCalls {
# 3. SPLIT INDIVIDUAL CALLS
echo -e "\n(3) SPLITTING INDIVIDUAL CALLS\n($SCRIPTS/splitMAandMNPs.py)\n"

mkdir $OUTDIR/2_individual_split
for FILE in `ls $OUTDIR/1_individual_calls/platypusVariants_*`; do
    $SCRIPTS/splitMAandMNPs.py $FILE >> $OUTDIR/logs/3_split.log
    mv "${FILE%.*}".split.vcf $OUTDIR/2_individual_split/
done
}

function identifyIndelSNPS {
# 4. IDENTIFY SNPS CLOSE TO INDELS IN ANY SAMPLE
echo -e "\n(4) IDENTIFYING SNPS CLOSE TO INDELS IN ANY SAMPLE\n($SCRIPTS/platypusIndelFlag.py)\n"

mkdir $OUTDIR/4_merged
ls -1 $OUTDIR/2_individual_split/* > $OUTDIR/2_individual_split/list.txt
$SCRIPTS/platypusIndelFlag.py $OUTDIR/2_individual_split/list.txt $OUTDIR/4_merged/indel_flagged_SNPs.txt > $OUTDIR/logs/4_indel_flag.log
}

function filterIndividualCalls {
# 5. FILTER INDIVIDUAL CALLS
echo -e "\n(5) FILTERING INDIVIDUAL CALLS\n"

mkdir $OUTDIR/3_individual_filtered
for FILE in `ls $OUTDIR/2_individual_split/platypusVariants_*`; do
    NAME=`basename $FILE`
    awk '!(($7 ~ /badReads/) || ($7 ~ /MQ/) || ($7 ~ /strandBias/) || ($7 ~ /SC/) || ($7 ~ /QD/))' $FILE > $OUTDIR/3_individual_filtered/"${NAME%.*}".filtered.vcf
done
}

function mergeVCFs {
# 6. MERGE THE VCFS ACCORDING TO CHR,POS,REF,ALT VALUES
echo -e "\n(6) EXTRACTING SNPS AND INDELS FILTERED INDIVIDUAL VCF FILES\n($SCRIPTS/platypusSNPmerge.py, $SCRIPTS/platypusIndelMerge.py)\n"

# Create list of filtered individual VCFs, extract SNPs, extract indels, append indels to Allele 1 SNPs file
ls -1 $OUTDIR/3_individual_filtered/* > $OUTDIR/3_individual_filtered/list.txt
$SCRIPTS/platypusSNPmerge.py $OUTDIR/3_individual_filtered/list.txt $OUTDIR/4_merged/indel_flagged_SNPs.txt $OUTDIR/4_merged > $OUTDIR/logs/5.1_SNP_merge.log
$SCRIPTS/platypusIndelMerge.py $OUTDIR/3_individual_filtered/list.txt $OUTDIR/4_merged > $OUTDIR/logs/5.2_indel_merge.log
cat $OUTDIR/4_merged/MergedIndels.vcf >> $OUTDIR/4_merged/MergedSNPs_allele1.vcf
}

function genotypeMergedVariants {
# 7. GENOTYPE THE MERGED VARIANTS IN ALL SAMPLES
echo -e "\n(7) GENOTYPING MERGED VARIANTS IN ALL SAMPLES\n($SCRIPTS/platypusExtractExons.py)\n"

for FILE in `ls $OUTDIR/4_merged/MergedSNPs_allele*.vcf`; do
    vcf-sort $FILE > "${FILE%.*}".sorted.vcf
    bgzip -f "${FILE%.*}".sorted.vcf
    tabix -p vcf "${FILE%.*}".sorted.vcf.gz
done

# Extract only regions (exons) with variants
mkdir -p $OUTDIR/5_genotyped
$SCRIPTS/platypusExtractExons.py $REGIONS $OUTDIR/4_merged/MergedSNPs_allele1.vcf $OUTDIR/4_merged/MergedSNPs_allele2.vcf $OUTDIR/4_merged/MergedSNPs_allele3.vcf $OUTDIR/5_genotyped > $OUTDIR/logs/6_extract_exons.log

# Create list of BAM files for Platypus
ls -1 $BAMSDIR/*.bam > $OUTDIR/5_genotyped/bam_list.txt

for IND in `seq 1 3`; do
    echo -e "\nGenotyping allele $IND\n"
    $PLATYPUS callVariants \
    --logFileName=$OUTDIR/logs/$(( 6 + $IND ))_genotype_first_allele${IND}.log \
    --refFile=$REFERENCE \
    --regions=$OUTDIR/5_genotyped/regions_allele${IND}.txt \
    --bamFiles=$OUTDIR/5_genotyped/bam_list.txt \
    --minPosterior=0 \
    --nCPU=8 \
    --minReads=3 \
    --source=$OUTDIR/4_merged/MergedSNPs_allele${IND}.sorted.vcf.gz \
    --getVariantsFromBAMs=0 \
    -o $OUTDIR/5_genotyped/GenotypedSNPs_first_allele${IND}.vcf
done
}

function genotypeMissingCalls {
# 8. GENOTYPE THE MISSING CALLS USING CONSTRAINED REGIONS
echo -e "\n(8) GENOTYPING THE MISSING CALLS USING CONSTRAINED REGIONS\n"

for IND in `seq 1 3`; do
    # Extract missing calls by comparing merged and genotyped VCFs
    tail -n +49 $OUTDIR/5_genotyped/GenotypedSNPs_first_allele${IND}.vcf | cut -f1,2 > $OUTDIR/geno_pos.txt
    cut -f1,2 $OUTDIR/4_merged/MergedSNPs_allele${IND}.vcf > $OUTDIR/merged_pos.txt
    grep -vxFf $OUTDIR/geno_pos.txt $OUTDIR/merged_pos.txt > $OUTDIR/coords.txt
    awk -v OFS="" '{print $1 ":" $2 "-" $2}' $OUTDIR/coords.txt > $OUTDIR/5_genotyped/SNPregions_allele${IND}.txt
    rm $OUTDIR/geno_pos.txt $OUTDIR/merged_pos.txt $OUTDIR/coords.txt
    
    echo -e "\nGenotyping allele $IND\n"
    $PLATYPUS callVariants \
    --logFileName=$OUTDIR/logs/$(( 9 + $IND ))_genotype_snpRegions_allele${IND}.log \
    --refFile=$REFERENCE \
    --bamFiles=$OUTDIR/5_genotyped/bam_list.txt \
    --regions=$OUTDIR/5_genotyped/SNPregions_allele${IND}.txt \
    --minPosterior=0 \
    --nCPU=8 \
    --minReads=3 \
    --source=$OUTDIR/4_merged/MergedSNPs_allele${IND}.sorted.vcf.gz \
    --getVariantsFromBAMs=0 \
    -o $OUTDIR/5_genotyped/GenotypedSNPs_second_allele${IND}.vcf
done
}

function genotypeIndelExcludedSNPS {
# 9. GENOTYPE INDEL-EXCLUDED SNPS IN ALL SAMPLES
echo -e "\n(9) GENOTYPING INDEL-EXCLUDED SNPS IN ALL SAMPLES\n($SCRIPTS/platypusExtractExons.py)\n"

for FILE in `ls $OUTDIR/4_merged/IndelExcludedSNPs_allele*.vcf`; do
    vcf-sort $FILE > "${FILE%.*}".sorted.vcf
    bgzip -f "${FILE%.*}".sorted.vcf
    tabix -p vcf "${FILE%.*}".sorted.vcf.gz
done

# Rename previous region files to avoid overwriting
for IND in `seq 1 3`; do
    mv $OUTDIR/5_genotyped/regions_allele${IND}.txt $OUTDIR/5_genotyped/regions_allele${IND}_nonExcluded.txt
    mv $OUTDIR/5_genotyped/SNPregions_allele${IND}.txt $OUTDIR/5_genotyped/SNPregions_allele${IND}_nonExcluded.txt
done

# Extract only regions (exons) with SNPs
$SCRIPTS/platypusExtractExons.py $REGIONS $OUTDIR/4_merged/IndelExcludedSNPs_allele1.vcf $OUTDIR/4_merged/IndelExcludedSNPs_allele2.vcf $OUTDIR/4_merged/IndelExcludedSNPs_allele3.vcf $OUTDIR/5_genotyped > $OUTDIR/logs/13
_extract_exons_excluded.log

for IND in `seq 1 3`; do
    echo -e "\nGenotyping allele $IND\n"
    $PLATYPUS callVariants \
    --logFileName=$OUTDIR/logs/$(( 13 + $IND ))_genotype_first_allele${IND}_indelExcluded.log \
    --refFile=$REFERENCE \
    --regions=$OUTDIR/5_genotyped/regions_allele${IND}.txt \
    --bamFiles=$OUTDIR/5_genotyped/bam_list.txt \
    --minPosterior=0 \
    --nCPU=8 \
    --minReads=3 \
    --source=$OUTDIR/4_merged/IndelExcludedSNPs_allele${IND}.sorted.vcf.gz \
    --getVariantsFromBAMs=0 \
    -o $OUTDIR/5_genotyped/GenotypedSNPs_indelExcluded_first_allele${IND}.vcf
done
}

function genotypeMissingIndels {
# 10. GENOTYPE THE MISSING INDEL-EXCLUDED CALLS USING CONSTRAINED REGIONS
echo -e "\n(10) GENOTYPING MISSING INDEL-EXCLUDED CALLS USING CONSTRAINED REGIONS\n"

for IND in `seq 1 3`; do
    # Extract missing calls by comparing merged and genotyped VCFs
    tail -n +49 $OUTDIR/5_genotyped/GenotypedSNPs_indelExcluded_first_allele${IND}.vcf | cut -f1,2 > $OUTDIR/geno_pos.txt
    cut -f1,2 $OUTDIR/4_merged/IndelExcludedSNPs_allele${IND}.vcf > $OUTDIR/merged_pos.txt
    grep -vxFf $OUTDIR/geno_pos.txt $OUTDIR/merged_pos.txt > $OUTDIR/coords.txt
    awk -v OFS="" '{print $1 ":" $2 "-" $2}' $OUTDIR/coords.txt > $OUTDIR/5_genotyped/SNPregions_allele${IND}.txt
    rm $OUTDIR/geno_pos.txt $OUTDIR/merged_pos.txt $OUTDIR/coords.txt
    
    echo -e "\nGenotyping allele ${IND}\n"
    $PLATYPUS callVariants \
    --logFileName=$OUTDIR/logs/$(( 16 + $IND ))_genotype_snpRegions_allele${IND}_indelExcluded.log \
    --refFile=$REFERENCE \
    --bamFiles=$OUTDIR/5_genotyped/bam_list.txt \
    --regions=$OUTDIR/5_genotyped/SNPregions_allele${IND}.txt \
    --minPosterior=0 \
    --nCPU=8 \
    --minReads=3 \
    --source=$OUTDIR/4_merged/IndelExcludedSNPs_allele${IND}.sorted.vcf.gz \
    --getVariantsFromBAMs=0 \
    -o $OUTDIR/5_genotyped/GenotypedSNPs_indelExcluded_second_allele${IND}.vcf
done
}

function filterIndelExcluded {
# 11. MERGE AND FILTER INDEL-EXCLUDED SNPS
echo -e "\n(11) MERGING AND FILTERING INDEL-EXCLUDED SNPS\n($SCRIPTS/platypusIndelRescueFilter.py)\n"

cp $OUTDIR/5_genotyped/GenotypedSNPs_indelExcluded_first_allele1.vcf $OUTDIR/5_genotyped/GenotypedSNPs_indelExcluded_merged.vcf
tail -n +49 $OUTDIR/5_genotyped/GenotypedSNPs_indelExcluded_first_allele2.vcf >> $OUTDIR/5_genotyped/GenotypedSNPs_indelExcluded_merged.vcf
tail -n +49 $OUTDIR/5_genotyped/GenotypedSNPs_indelExcluded_first_allele3.vcf >> $OUTDIR/5_genotyped/GenotypedSNPs_indelExcluded_merged.vcf
tail -q -n +49 $OUTDIR/5_genotyped/GenotypedSNPs_indelExcluded_second_allele*.vcf >> $OUTDIR/5_genotyped/GenotypedSNPs_indelExcluded_merged.vcf

awk '!(($7 ~ /badReads/) || ($7 ~ /MQ/) || ($7 ~ /strandBias/) || ($7 ~ /SC/) || ($7 ~ /QD/))' $OUTDIR/5_genotyped/GenotypedSNPs_indelExcluded_merged.vcf > $OUTDIR/5_genotyped/GenotypedSNPs_indelExcluded_merged.filtered.
vcf
$SCRIPTS/platypusIndelRescueFilter.py $OUTDIR/5_genotyped/GenotypedSNPs_indelExcluded_merged.filtered.vcf > $OUTDIR/logs/20_indel_rescued_filter.log
}

function mergeGenotypedCalls {
# 12. MERGE AND FILTER GENOTYPED CALLS
echo -e "\n(12) MERGING AND FILTERING ALL GENOTYPED CALLS\n($SCRIPTS/platypusVAFFilter.py)\n"

cp $OUTDIR/5_genotyped/GenotypedSNPs_first_allele1.vcf $OUTDIR/5_genotyped/GenotypedVariants_merged.vcf
tail -n +49 $OUTDIR/5_genotyped/GenotypedSNPs_first_allele2.vcf >> $OUTDIR/5_genotyped/GenotypedVariants_merged.vcf
tail -n +49 $OUTDIR/5_genotyped/GenotypedSNPs_first_allele3.vcf >> $OUTDIR/5_genotyped/GenotypedVariants_merged.vcf
tail -q -n +49 $OUTDIR/5_genotyped/GenotypedSNPs_second_allele*.vcf >> $OUTDIR/5_genotyped/GenotypedVariants_merged.vcf
grep -v '#' $OUTDIR/5_genotyped/GenotypedSNPs_indelExcluded_merged.filtered.vcf |grep -v 'alleleBias' >> $OUTDIR/5_genotyped/GenotypedVariants_merged.vcf

awk '!(($7 ~ /badReads/) || ($7 ~ /MQ/) || ($7 ~ /strandBias/) || ($7 ~ /SC/) || ($7 ~ /QD/) || ($8 ~ /CGTACACGTA/))' $OUTDIR/5_genotyped/GenotypedVariants_merged.vcf > $OUTDIR/5_genotyped/GenotypedVariants_merged.filter
ed.vcf

$SCRIPTS/splitMAandMNPs.py $OUTDIR/5_genotyped/GenotypedVariants_merged.filtered.vcf 
vcf-sort $OUTDIR/5_genotyped/GenotypedVariants_merged.filtered.split.vcf > $OUTDIR/GenotypedVariants_final.vcf


echo -e "\nAll done!\nOutput in $OUTDIR/GenotypedVariants_final.vcf\n"
}

runPlatypus
runPlatypusMinFlank
splitCalls
identifyIndelSNPS
filterIndividualCalls
mergeVCFs
genotypeMergedVariants
genotypeMissingCalls
genotypeIndelExcludedSNPS
genotypeMissingIndels
filterIndelExcluded
mergeGenotypedCalls