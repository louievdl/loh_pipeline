#!/bin/bash
# usage: bash submit_array.sh -g hg19 -r /SAN/colcc/alex_work/samples_analysis -c "" -s configuration/2021-06_icgc_cancer_vs_normal.txt

# some defaults
GENOME="hg19" # default; must be hg19 or grch38
CHR_SUFFIX="_nochr" # default; must be either "_nochr" OR "_chr"; genome files (fasta, bed) are of the form grch38_nochr*, hg19_chr*, etc
RESULTS_BASE_DIR=/SAN/colcc/alex_work/samples_analysis
SAMPLE_METADATA=configuration/2021-06_icgc_cancer_vs_normal.txt

# set variables from arguments
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" # https://stackoverflow.com/questions/59895/how-can-i-get-the-source-directory-of-a-bash-script-from-within-the-script-itsel
PPLN_BASE_DIR=$SCRIPT_DIR/..
while getopts g:r:s: flag
do
    case "${flag}" in
        g) GENOME=${OPTARG};;            # 'hg19' or 'grch38'
        c) CHR_SUFFIX=${OPTARG};;      # '_nochr' or '' if chromosomes are specified without or with 'chr' respectively; used to choose appropriately named genomic bed and fasta files
        r) RESULTS_BASE_DIR=${OPTARG};;  # base directory for analysis results
        s) SAMPLE_METADATA=${OPTARG};;   # sample metadata file
    esac
done
echo "Genome: $GENOME";
echo "chr in file name: $CHR_SUFFIX";
echo "Results location: $RESULTS_BASE_DIR";
echo "Sample metadata: $SAMPLE_METADATA";

if [ ! -d $RESULTS_BASE_DIR ]
then
    mkdir $RESULTS_BASE_DIR
fi

# check for finished jobs, with lohhla outputs
# identify pt ids (subdirs within RESULTS_BASE_DIR) with completed figures
# use found pt ids to search first col in configuration/xx.txt; create file of grep patterns; empty file ok; do not repeat analysis in case of success
SAMPLE_METADATA_TMP=$SAMPLE_METADATA".tmp"
SUCCESS_IDS_TMP=$SAMPLE_METADATA".success"
ls $RESULTS_BASE_DIR/*/lohhla/results/Figures/*.HLA.pdf | \
    perl -ne 'chomp; s/.lohhla.results.Figures..*$//; s/^.*\///; print "^$_\t\n"' > $SUCCESS_IDS_TMP
grep -Ev -f $SUCCESS_IDS_TMP $SAMPLE_METADATA > $SAMPLE_METADATA_TMP
rm $SUCCESS_IDS_TMP
# note, upon lohhla failure, Figures dirs will not exist, so no deletion necessary before rerunning job.

# calculate number of tasks required
SAMPLES_LINES=($(cat $SAMPLE_METADATA_TMP | cut -f 1)) # includes a header
NUM_LINES=${#SAMPLES_LINES[@]} # array length
TASK_COUNT=$((NUM_LINES-1))

echo "GENOME=$GENOME SAMPLE_METADATA=$SAMPLE_METADATA_TMP TASK_COUNT=$TASK_COUNT "
#SGE_TASK_ID=2 && . ./job.sh # extra . means run-in-this-env
#exit 0

# qsub the array job
# Request 1 core, 5 hour runtime (esp if 'chr' must be cleared from bam files), 20GB RAM, job array, 2 tasks run recurrently
qsub -cwd -pe smp 1 -l h_rt=5:0:0,h_vmem=20G,tmem=20G -t 1-$TASK_COUNT -tc 2 -v GENOME=$GENOME,CHR_SUFFIX=$CHR_SUFFIX,RESULTS_BASE_DIR=$RESULTS_BASE_DIR,SAMPLE_METADATA=$SAMPLE_METADATA_TMP job.sh

