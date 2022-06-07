###########################################################################
# Bash script main.sh (v1.0) to call other scripts
# How to run: nohup bash 'main.sh' > 'log.txt' 2>&1 &
# Last update: 2022_06_07
# Author: Michelle Almeida da Paz
###########################################################################

#! /usr/bin/env bash
# Exit the script if any statement returns a non-true return value
set -o errexit
# Use by default the error status of the last item in a pipeline
set -o pipefail
# Exit the script if use an uninitialised variable
set -o nounset

error_exit() {
	line=$1
	shift 1
	echo "ERROR: non zero return code from line: $line - $@" 1>&2
	exit 1
}

help_function() {
	echo "";
	echo "Usage bash $0";
	echo -e "\tspecies\thg38 (Homo sapiens) or mm10 (Mus musculus)";
	echo -e "\titerations\tnumber_of_iterations (e.g. 100)";
	echo -e "\talpha\tthreshold (e.g. 0.05)";
	echo -e "\tenrichment\tthreshold (e.g. 1.0)";
	echo -e "\tfilter\t0 (no - do not filter high signal regions) or 1 (yes - filter high signal regions)";
	exit 1; # Exit the script after printing help message
}

## Set path folders
# Local work directory
LOCALDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )";
WORKDIR="$LOCALDIR";
[ -d $WORKDIR  ] || error_exit "$LINENO" "Cannot find $WORKDIR!";
# Scripts directory
CODEDIR="$WORKDIR/scripts";
[ -d $CODEDIR  ] || error_exit "$LINENO" "Cannot find $CODEDIR!";
# Data directory
DATASET="$WORKDIR/bam";
[ -d $DATASET  ] || error_exit "$LINENO" "Cannot find $DATASET!";
# Output directory
OUTDIR="$WORKDIR/results";
[ -d $OUTDIR  ] || error_exit "$LINENO" "Cannot find $OUTDIR!";

## Set path files
# Path to dataset file
PATH_FILE="$WORKDIR/references/path_dataset.csv";
# Path to control/sample names file
INFO="$WORKDIR/references/control_sample.txt";
[ -f $INFO  ] || error_exit "$LINENO" "Cannot find $INFO!";
# Path to parameters file
PARAMETERS="$WORKDIR/references/parameters.txt";
[ -f $PARAMETERS  ] || error_exit "$LINENO" "Cannot find $PARAMETERS!";

## Set parameters
# Set species
SPECIES=`grep "species" $PARAMETERS | cut -f2`;
if [[ ! $SPECIES =~ ^(hg38|mm10)$ ]]; then
	help_function
	exit 1
fi
# Set iterations
int='^[0-9]+$';
ITER=`grep "iterations" $PARAMETERS | cut -f2`;
if ! [[ "$ITER" =~ $int ]]; then
	help_function
	exit 1
fi
# Set threshold alpha
float='^[0-9]+\.?[0-9]*$';
ALPHA=`grep "alpha" $PARAMETERS | cut -f2`;
if ! [[ "$ALPHA" =~ $float ]]; then
	help_function
	exit 1
fi
# Set threshold enrichment
ENRICHMENT=`grep "enrichment" $PARAMETERS | cut -f2`;
if ! [[ "$ENRICHMENT" =~ $float ]]; then
	help_function
	exit 1
fi
# Set filter (no/yes)
FILTER_HIGH_SIGNAL=`grep "filter" $PARAMETERS | cut -f2`;
if [[ ! $FILTER_HIGH_SIGNAL =~ ^(0|1)$ ]]; then
	help_function
	exit 1
fi

## Path to repeat file
REPEATS="$WORKDIR/repeats/rmsk_$SPECIES.bed";
[ -f $REPEATS  ] || error_exit "$LINENO" "Cannot find $REPEATS!";

## Print message
echo "Running..."
echo "Species: $SPECIES";
echo "Number of iterations: $ITER";
echo "Thresholds: alpha $ALPHA and log2FC $ENRICHMENT";
echo "Filter high signal regions (0-no, 1-yes): $FILTER_HIGH_SIGNAL";

## Preprocess dataset
if [ 2 -gt 1 ]; then

	if [ -f $PATH_FILE ]; then
		rm $PATH_FILE;
	fi

	for BAMFILE in $DATASET/*.bam; do 
		BAMFULLNAME="${BAMFILE##*/}";
		SAMPLENAME="${BAMFULLNAME%.*}";
		INTERM="$OUTDIR/$SAMPLENAME/"`basename $SAMPLENAME`"_interm.txt";
		PRE_BED="$OUTDIR/$SAMPLENAME/"`basename $SAMPLENAME`"_pre.bed";
		FILTERED_BED="$OUTDIR/$SAMPLENAME/"`basename $SAMPLENAME`"_filtered.bed";
		BED="$OUTDIR/$SAMPLENAME/"`basename $SAMPLENAME`".bed";
		DICT="$OUTDIR/$SAMPLENAME/"`basename $SAMPLENAME`".dict";
		NR_READS_TO_SUBSAMPLE=20000000;

		# Create an output folder for each sample
		if ! [[ -d $OUTDIR/$SAMPLENAME ]]; then
			echo "Creating output folder for $SAMPLENAME..."
			mkdir $OUTDIR/$SAMPLENAME;
		fi
		# Remove CHR_M, CHR_VAR and unmapped reads
		echo "# Filtering chrM, chr_var, and unmapped reads...";
		samtools view $BAMFILE | awk '{if($3 ~ "chr"){if($3 !~ "_"){if($3 !~ "chrM"){printf "%s\t%d\t%s\t%s\n", $3, $4-1, $1, $10;}}}}' > $INTERM;
		head -n3 $INTERM;
		# Calculate the read length average
		echo "# Calculating the read length average...";
		read_len=`(cat $INTERM | head -n10000; true) | awk 'BEGIN{total=0;count=0;}{if(length($4) > 1){total+=length($4); count++;}}END{print int(total/count);}'`;
		echo "# Read length of $SAMPLENAME: $read_len";
		# Reformat to BED file
		echo "# Reformating to BED file...";
		cat $INTERM | awk -v rl=$read_len '{printf "%s\t%d\t%d\t%s\n", $1, $2, $2+rl, $3;}' | sort -k1,1 -k2,2n > $PRE_BED;
		head -n3 $PRE_BED;
		echo "# Output: $PRE_BED";
		# Filter high signal regions
		if [ $FILTER_HIGH_SIGNAL == 1 ]; then
			WINDOWS=100;
			STEPS=50;
			PERCENTILE=0.99997;
			GENOME="$WORKDIR/references/$SPECIES.genome";
			GENOME_IN_WINDOWS="$WORKDIR/references/${SPECIES}_${WINDOWS}_${STEPS}.bed";
			SAMPLE_IN_WINDOWS="$OUTDIR/$SAMPLENAME/"`basename $SAMPLENAME`"_${WINDOWS}_${STEPS}.bed";

			# Systematic sliding window approach
			echo "Creating windows of size $WINDOWS and steps of size $STEPS...";
			bedtools makewindows -g $GENOME -w $WINDOWS -s $STEPS > $GENOME_IN_WINDOWS;
			echo "# Output: $GENOME_IN_WINDOWS";
			# Intersect and count for each window
			echo "Intersecting genome in windows file with BED file and counting the number of reads for each window";
			bedtools intersect -a $GENOME_IN_WINDOWS -b $PRE_BED -c > $SAMPLE_IN_WINDOWS;
			echo "# Output: $SAMPLE_IN_WINDOWS";
			# Identify high peaks
			echo "Identifying high peak regions...";
			if ! [[ -d $OUTDIR/$SAMPLENAME/R_output ]]; then
				echo "Creating folder for R output $SAMPLENAME..."
				mkdir $OUTDIR/$SAMPLENAME/R_output;
			fi
			PATH_IN_R="$OUTDIR/$SAMPLENAME/R_output";
			INPUT_IN_R="$SAMPLE_IN_WINDOWS";
			R --slave --no-save --no-restore --no-environ --silent --args $PERCENTILE $PATH_IN_R $INPUT_IN_R $SPECIES < $CODEDIR/filter_highcov.R;
			# Merge overlapping regions
			echo "Merging overlapping high peak regions...";
			for CHR_COUNT in $PATH_IN_R/*.bed; do
				CHR_COUNT_FULLNAME="${CHR_COUNT##*/}";
				CHR_COUNT_NAME="${CHR_COUNT_FULLNAME%.*}";
				CHR_COUNT_OUT="$PATH_IN_R/"`basename $CHR_COUNT_NAME`"_merged.bed";
				bedtools merge -i $CHR_COUNT > $CHR_COUNT_OUT;
			done
			# Create high Signal Region File
			HIGH_SIGNAL="$PATH_IN_R/${SAMPLENAME}_high_signal_region_${WINDOWS}_${STEPS}.bed";
			# Concatenate all merged windows of all chromosomes
			echo "Concatenate all merged windows of all chromosomes...";
			cat $PATH_IN_R/*_merged.bed > $PATH_IN_R/${SAMPLENAME}_high_signal_region_${WINDOWS}_${STEPS}.bed
			# Sort high signal region bed file
			echo "Sort high signal region bed file...";
			sort-bed --max-mem 16G $PATH_IN_R/${SAMPLENAME}_high_signal_region_${WINDOWS}_${STEPS}.bed > $PATH_IN_R/${SAMPLENAME}_high_signal_region_${WINDOWS}_${STEPS}_sorted.bed
			input1="$REPEATS";
			input2="${input1##*/}";
			input3="${input2%.txt}";
			ANNOTATION="$WORKDIR/repeats/"`basename $input3`"_filtered_grouped.bed";
			# Annotate the high signal regions in regards of repeat sequences
			echo "Annotate the high signal regions in regards of repeat sequences...";
			bedmap --echo --echo-map --echo-overlap-size $ANNOTATION $PATH_IN_R/${SAMPLENAME}_high_signal_region_${WINDOWS}_${STEPS}_sorted.bed > $PATH_IN_R/${SAMPLENAME}_high_signal_region_${WINDOWS}_${STEPS}_annotation.bed
			cat $PATH_IN_R/${SAMPLENAME}_high_signal_region_${WINDOWS}_${STEPS}_annotation.bed | grep "|chr" > $PATH_IN_R/${SAMPLENAME}_high_signal_region_${WINDOWS}_${STEPS}_final_annotation.bed
			# Filter out high signal regions
			echo "# Filter out high signal regions...";
			bedtools intersect -a $PRE_BED -b $HIGH_SIGNAL -v > $FILTERED_BED;
			head -n3 $FILTERED_BED;
			echo "# Output: $FILTERED_BED";
		else
			cat $PRE_BED > $FILTERED_BED;
		fi
		# Calculate the number of reads
		echo "# Calculating the number of reads for $SAMPLENAME...";
		ACTUAL_NR_READS=`(cat $FILTERED_BED | cut -f4) | awk '{c[$0]++} END {for (line in c) print c[line], line}' | wc -l`;
		echo "# Read number of $SAMPLENAME: $ACTUAL_NR_READS";
		if [ $ACTUAL_NR_READS -gt $NR_READS_TO_SUBSAMPLE ]; then
			echo "# Subsampling FILTERED_BED file...";
			perl $CODEDIR/subsample.pl $FILTERED_BED $NR_READS_TO_SUBSAMPLE > $BED;
			head -n3 $BED;
			echo "# Output: $BED";
			ACTUAL_NR_READS=$NR_READS_TO_SUBSAMPLE;
			echo "# New library size of $SAMPLENAME: $ACTUAL_NR_READS";
		else
			cat $FILTERED_BED > $BED;
			head -n3 $BED;
			echo "# Output: $BED";
		fi
		# Create DICT file
		echo "# Creating DICT file..."; 
		cat $BED | cut -f1,4 > $DICT
		head -n3 $DICT;
		echo "# Output: $DICT";
		# Create path file 
		echo "Feeding PATH_FILE file..."; 
		echo -e "$BAMFILE;$ACTUAL_NR_READS;$read_len" >> $PATH_FILE;
		# Remove temporary files
		rm $INTERM;
		rm $PRE_BED;
		rm $FILTERED_BED;
	done
fi

## Count for TE families/subfamilies
if [ 2 -gt 1 ]; then
	for BAMFILE in `cut -d ';' -f1 $PATH_FILE | grep -P -v "^#"`; do
		[ -f $BAMFILE  ] || error_exit "$LINENO" "Cannot find $BAMFILE!";
		BAMFULLNAME="${BAMFILE##*/}";
		SAMPLENAME="${BAMFULLNAME%.*}";
		BED="$OUTDIR/$SAMPLENAME/"`basename $SAMPLENAME`".bed";
		[ -f $BED  ] || error_exit "$LINENO" "Cannot find $BED!";
		INTERSECT="$OUTDIR/$SAMPLENAME/"`basename $SAMPLENAME`"_intersected.txt";
		COUNTS="$OUTDIR/$SAMPLENAME/"`basename $SAMPLENAME`"_counts.txt";
		nr_reads=`grep "$BAMFILE" $PATH_FILE | cut -d ';' -f2`;
		read_len=`grep "$BAMFILE" $PATH_FILE | cut -d ';' -f3`;
		echo "Analysis $SAMPLENAME... $nr_reads reads";
		SECONDS=0;
		# Intersect BED file and FILTERED REPEATS
		[ -f $REPEATS  ] || error_exit "$LINENO" "Cannot find $REPEATS!";
		echo "# Intersecting BED file and REPEATS...";
		bedmap --echo --echo-map --echo-overlap-size $REPEATS $BED > $INTERSECT;
		head -n3 $INTERSECT;
		# Count
		echo "# Counting for TE families/subfamilies...";
		perl $CODEDIR/count.pl $read_len $BED $INTERSECT > $COUNTS;
		head -n3 $COUNTS;

		if (( $SECONDS > 3600 )) ; then
			let "hours=SECONDS/3600"
			let "minutes=(SECONDS%3600)/60"
			let "seconds=(SECONDS%3600)%60"
			echo "Completed in $hours hour(s), $minutes minute(s) and $seconds second(s)" 
		elif (( $SECONDS > 60 )) ; then
			let "minutes=(SECONDS%3600)/60"
			let "seconds=(SECONDS%3600)%60"
			echo "Completed in $minutes minute(s) and $seconds second(s)"
		else
			echo "Calculations completed for $BAMFILE in $SECONDS seconds"
		fi
	done
fi

## Run T3E
if [ 2 -gt 1 ]; then
	for CONTROL in `cut -f1 $INFO | grep -P -v "^#"`; do
		CONTROL_BED="$OUTDIR/$CONTROL/"`basename $CONTROL`".bed";
		CONTROL_COUNTS="$OUTDIR/$CONTROL/"`basename $CONTROL`"_counts.txt";
		# Calculate input-based probability distribution
		if ! [[ -d $OUTDIR/$CONTROL/probabilities ]]; then
			READ_LEN=`grep "$CONTROL" $PATH_FILE | cut -d ';' -f3`;
			echo "Calculating input-based probabilities of... $CONTROL";
			mkdir $OUTDIR/$CONTROL/probabilities;
			FOLDER_PROBABILITY_CONTROL="$OUTDIR/$CONTROL/probabilities";
			python3 -u $CODEDIR/probabilities.py --control $CONTROL_BED --readlen $READ_LEN --species $SPECIES --outputfolder $FOLDER_PROBABILITY_CONTROL;
		fi
		FOLDER_PROBABILITY_CONTROL="$OUTDIR/$CONTROL/probabilities";
		SAMPLES=`grep "$CONTROL" $INFO | cut -f2`;
		IFS=',' read -ra SAMPLES_ARRAY <<< "$SAMPLES"
		for SAMPLE in "${SAMPLES_ARRAY[@]}"; do
			echo "Running T3E to sample $SAMPLE";
			SAMPLE_DICT="$OUTDIR/$SAMPLE/"`basename $SAMPLE`".dict";
			[ -f $REPEATS  ] || error_exit "$LINENO" "Cannot find $REPEATS!";
			[ -f $SAMPLE_DICT  ] || error_exit "$LINENO" "Cannot find $SAMPLE_DICT!";
			[ -f $CONTROL_BED  ] || error_exit "$LINENO" "Cannot find $CONTROL_BED!";
			[ -f $CONTROL_COUNTS  ] || error_exit "$LINENO" "Cannot find $CONTROL_COUNTS!";
			READ_LEN=`grep "$SAMPLE" $PATH_FILE | cut -d ';' -f3`;
			# Run Simulations T3E and Count for TE families/subfamilies
			python3 -u $CODEDIR/t3e.py --repeat $REPEATS --sample $SAMPLE_DICT --readlen $READ_LEN --control $CONTROL_BED --controlcounts $CONTROL_COUNTS --probability $FOLDER_PROBABILITY_CONTROL --iter $ITER --species $SPECIES --outputfolder $OUTDIR/$SAMPLE/ --outputprefix $SAMPLE > $WORKDIR/log_$SAMPLE.txt;
			BACKGROUND="$OUTDIR/$SAMPLE/"`basename $SAMPLE`"_background.txt";
			[ -f $BACKGROUND  ] || error_exit "$LINENO" "Cannot find $BACKGROUND!";
			SAMPLE_COUNTS="$OUTDIR/$SAMPLE/"`basename $SAMPLE`"_counts.txt";
			[ -f $SAMPLE_COUNTS  ] || error_exit "$LINENO" "Cannot find $SAMPLE_COUNTS!";
			# Calculate enrichment for TE families/subfamilies
			echo "Enrichment of sample $SAMPLE";
			python3 -u $CODEDIR/enrichment.py --background $BACKGROUND --signal $SAMPLE_COUNTS --iter $ITER --alpha $ALPHA --enrichment $ENRICHMENT --outputfolder $OUTDIR/$SAMPLE/ --outputprefix $SAMPLE;
		done
	done
	exit 0;
fi
