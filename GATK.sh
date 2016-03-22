#!/bin/sh

#  GATK.sh
#  Pipeline for the GATK workflow for calling variants.
#  Uses HaplotypeCaller.
#
#  Created by Janu Verma on 10/24/14.
#  jv367@cornell.edu



##################################################################
#   SOFTWARES
# Provide locations of the softwares to be used. 

PICARD="java -Xmx50g -Djava.io.tmpdir=`pwd`/tmp -jar  /program/picard-tools-2.1.1/picard.jar"
SAMTOOLS="samtools"
GATK="java -jar /program/gatk/GenomeAnalysisTK.jar"


#################################################################
#   FILES 
# Aligned BAM file, reference FASTA file and VCF file of known variants. 

BAM="SRR1984571.sorted.bam"
REFERENCE="/ref/analysis/juntaehwan/ref/Cannuum/Pepper.v.1.55.total.chr.fa"
#KNOWN="known_sites.vcf"


#####################################################################
#   Check if the BAM file satisfies the requirements of the GATK

$SAMTOOLS view -H $BAM > alpha.txt
head_file="alpha.txt"
grep '^@SQ' $head_file > dot_com.txt
file="dot_com.txt"


echo "created $file"

#less $file
# check if the file is sorted.
T="sort -c -t':' -nk2 $file"
if [ "$T" ]; then
    echo "The file is sorted!"

else
    echo "The file is not sorted."
    echo "Sorting........."
    $PICARD SortSam INPUT=$BAM OUTPUT="${BAM%.bam}.sorted.bam" SORT_ORDER=coordinate
    mv "${BAM%.bam}.sorted.bam" $BAM

fi


#   check if the file contains RG information.
if grep -q '^@RG' $head_file; then
    echo "The file contains RG information."
else
    echo "The file does not contain RG information and the GATK will not work!"
    exit
fi



###################################################################
#   Prepping reference geonome fasta file for GATK
echo "Prepping reference geonome fasta file for GATK....."

# Create sequence dictionary using Picard Tools.
# the following command produces a SAM-style header file describing the contents of our fasta file.
$PICARD CreateSequenceDictionary \
REFERENCE=$REFERENCE \
OUTPUT="${REFERENCE%.fa}.dict"

echo "created sequence dictionary ${Reference%.fa}.dict for the reference genome."

echo "indexing the reference genome...."

# Create the fasta index file.
# The index file describes byte offset in the fasta file for each contig. It is a text file with one record
# per line for each of the fasta contigs. Each record is of the type -
# contig, size, location, basePerLine, bytesPerLine
$SAMTOOLS faidx $REFERENCE

echo "Reference genome is now ready for GATK."




###############################################################
## Summary Statistics

$PICARD MeanQualityByCycle \
INPUT=$BAM \
CHART_OUTPUT=mean_quality_by_cycle.pdf \
OUTPUT=read_quality_by_cycle.txt \
REFERENCE_SEQUENCE=$REFERENCE \
VALIDATION_STRINGENCY=LENIENT


$PICARD QualityScoreDistribution \
INPUT=$BAM \
CHART_OUTPUT=mean_quality_overall.pdf \
OUTPUT=read_quality_overall.txt \
REFERENCE_SEQUENCE=$REFERENCE \
VALIDATION_STRINGENCY=LENIENT

$PICARD CollectWgsMetrics \
INPUT=$BAM OUTPUT=stats_picard.txt \
REFERENCE_SEQUENCE=$REFERENCE \
MINIMUM_MAPPING_QUALITY=20 \
MINIMUM_BASE_QUALITY=20 \
VALIDATION_STRINGENCY=LENIENT


#############################################################
# Mark duplicate reads.

echo "mark the duplicates in the bam file."

$PICARD MarkDuplicates INPUT=$BAM OUTPUT="${BAM%.bam}_dups_marked.bam" \
METRICS_FILE="${BAM%.bam}_dups_metrics.txt" REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT

echo "index the dup-marked bam file,"${BAM%.bam}_dups_marked.bam" "

$SAMTOOLS index "${BAM%.bam}_dups_marked.bam"


BAM_FILE="${BAM%.bam}_dups_marked.bam"

#############################################################
## GATK Data Pre-Processing

# Step 1 - Local realignment around indels.
# Create a target list of intervals to be realigned.

echo "Creating a target list of intervals to be realigned...."

$GATK \
-T RealignerTargetCreator \
-R $REFERENCE \
-I $BAM_FILE \
-o "${BAM%.bam}_target_intervals.list"

# do the local realignment.
echo "local realignment..."

$GATK \
-T IndelRealigner \
-R $REFERENCE \
-I $BAM_FILE \
-targetIntervals "${BAM%.bam}_target_intervals.list" \
-o "${BAM%.bam}_realigned_reads.bam"

echo "indexing the realigned bam file..."

# Create a new index file.
$SAMTOOLS index "${BAM%.bam}_realigned_reads.bam"






# Step 2 - Base recalibration (fixes them so they better reflect the probability of mismatching the genome).
# Analyze patterns of covariation in the sequence.
#echo "base recalibration...."

#$GATK \
#-T BaseCalibrator \
#-R $REFERENCE \
#-I "${BAM%.bam}_realigned_reads.bam" \
#-knownSites $KNOWN \
#-o "${BAM%.bam}_recal_data.table"


# Do a second pass to analyze covariation remaining after recalibration
#echo "base recalibration - part 2..."

#$GATK \
#-T BaseCalibrator \
#-R $REFERENCE \
#-I "${BAM%.bam}_realigned_reads.bam" \
#-knownSites $KNOWN \
#-BQSR "${BAM%.bam}_recal_data.table" \
#-o "${BAM%.bam}_post_recal_data.table"


# Generate before/after plots
#$GATK \
#-T AnalyzeCovariates \
#-R $REFERENCE \
#-before "${BAM%.bam}_recal_data.table" \
#-after "${BAM%.bam}_post_recal_data.table" \
#-plot "${BAM%.bam}_recalibration_plots.pdf"


# Apply recalibration to the sequence data.

#echo "recalibrating the sequence data.."

#$GATK \
#-T PrintReads \
#-R $REFERENCE \
#-I "${BAM%.bam}_realigned_reads.bam" \
#-BQSR "${BAM%.bam}_recal_data.table" \
#-o "${BAM%.bam}_recal_reads.bam"



ln -s "${BAM%.bam}_realigned_reads.bam" "${BAM%.bam}_recal_reads.bam" # Skipping recal step as there're no known sites
ln -s "${BAM%.bam}_realigned_reads.bam.bai" "${BAM%.bam}_recal_reads.bam.bai" 

###########################################################################
# GATK Variant Calling -  HaplotypeCaller
# Set -nct, outmode, emit_thresh, call_threh,

outmode="EMIT_ALL_CONFIDENT_SITES"
emit_thresh=20	#Threshold for tagging possible variants
call_thresh=30	#Threshold for tagging _good_ variants
hetrate=0.03	#Popgen heterozygosity rate (that is, for any two random chrom in pop, what is rate of mismatch).               Human is ~0.01, so up maize to ~0.03
minBaseScore=20	#Minimum Phred base score to count a base (20 = 0.01 error, 30=0.001 error, etc)

echo "calling variants...."

$GATK \
-T HaplotypeCaller \
-R $REFERENCE \
-I "${BAM%.bam}_recal_reads.bam" \
--emitRefConfidence GVCF \
--variant_index_type LINEAR \
--variant_index_parameter 128000 \
-hets $hetrate \
-mbq $minBaseScore \
-stand_emit_conf $emit_thresh \
-stand_call_conf $call_thresh \
-out_mode $outmode \
-nct 15 \
-o "${BAM%.bam}_output.raw.snps.indels.g.vcf"
