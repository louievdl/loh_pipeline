#!/bin/bash
# usage: bash submit_array.sh -g hg19 -r /SAN/colcc/alex_work/samples_analysis -c nochr -s configuration/2021-06_icgc_cancer_vs_normal.txt
#
# or specify a SGE_TASK_ID, in which case run interactively rather than submitting, in which case
# usage: bash submit_array.sh -g hg19 -r /SAN/colcc/alex_work/samples_analysis -c nochr -s configuration/2021-06_icgc_cancer_vs_normal.txt -t 57 >stdout.txt 2>stderr.txt &

# some defaults
GENOME="hg19" # default; must be hg19 or grch38
CHR_SUFFIX="nochr" # default; must be either "nochr" OR "chr"; genome files (fasta, bed) are of the form grch38_nochr*, hg19_chr*, etc
RESULTS_BASE_DIR=/SAN/colcc/alex_work/samples_analysis
SAMPLE_METADATA=configuration/2021-06_icgc_cancer_vs_normal.txt
SGE_TASK_ID=-1

# set variables from arguments
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" # https://stackoverflow.com/questions/59895/how-can-i-get-the-source-directory-of-a-bash-script-from-within-the-script-itsel
PPLN_BASE_DIR=$SCRIPT_DIR/../..
while getopts g:c:r:s:t: flag
do
    case "${flag}" in
        g) GENOME=${OPTARG};;            # 'hg19' or 'grch38'
        c) CHR_SUFFIX=${OPTARG};;        # 'chr' or 'nochr' if chromosomes are specified with or without 'chr' respectively; used to choose appropriately named genomic bed and fasta files
        r) RESULTS_BASE_DIR=${OPTARG};;  # base directory for analysis results
        s) SAMPLE_METADATA=${OPTARG};;   # sample metadata file
        t) SGE_TASK_ID=${OPTARG};;       # run a specific task (line number in the file, after removal of successfully completed jobs)
    esac
done
echo "Genome: $GENOME";
echo "chr in file name: $CHR_SUFFIX";
echo "Results location: $RESULTS_BASE_DIR";
echo "Sample metadata: $SAMPLE_METADATA";
## NOTE: if underscore is placed at beginning of a text string, seems to need double quotes around it when assigning to a variable: echo works ok, but passing it to a script doesn't seem to, so avoid if possible


if [ ! -d $RESULTS_BASE_DIR ]
then
    mkdir $RESULTS_BASE_DIR
fi

# if a job is specified (only to be done if running interactively), analyse chosen sample and exit the script
if [ $SGE_TASK_ID -gt -1 ]
then
    #export GENOME=$GENOME
    #export CHR_SUFFIX=$CHR_SUFFIX
    #export RESULTS_BASE_DIR=$RESULTS_BASE_DIR
    #export SAMPLE_METADATA=$SAMPLE_METADATA
    #export SGE_TASK_ID=$SGE_TASK_ID
    #export PPLN_BASE_DIR=$PPLN_BASE_DIR
    . ./job.sh # extra . means run-in-this-env
    exit 0
fi

# otherwise, check for finished jobs, with lohhla outputs
SAMPLE_METADATA_TMP=$SAMPLE_METADATA".tmp"
SUCCESS_IDS_TMP=$SAMPLE_METADATA".success"
if [ `ls $RESULTS_BASE_DIR/*/lohhla/results/Figures/*.HLA.pdf | grep -c ^` -gt 0 ]
then
    # identify pt ids (subdirs within RESULTS_BASE_DIR) with completed figures
    # use found pt ids to search first col in configuration/xx.txt; create file of grep patterns; empty file ok; do not repeat analysis in case of success
    ls $RESULTS_BASE_DIR/*/lohhla/results/Figures/*.HLA.pdf | \
       perl -ne 'chomp; s/.lohhla.results.Figures..*$//; s/^.*\///; print "^$_\t\n"' > $SUCCESS_IDS_TMP
    grep -Ev -f $SUCCESS_IDS_TMP $SAMPLE_METADATA > $SAMPLE_METADATA_TMP
    rm $SUCCESS_IDS_TMP
else
    cp $SAMPLE_METADATA $SAMPLE_METADATA_TMP
fi

# calculate number of tasks required
SAMPLES_LINES=($(cat $SAMPLE_METADATA_TMP | cut -f 1)) # includes a header
NUM_LINES=${#SAMPLES_LINES[@]} # array length
TASK_COUNT=$((NUM_LINES-1))

echo "GENOME=$GENOME SAMPLE_METADATA=$SAMPLE_METADATA_TMP RESULTS_BASE_DIR=$RESULTS_BASE_DIR TASK_COUNT=$TASK_COUNT SGE_TASK_ID=$SGE_TASK_ID "

# qsub the array job
# Request 1 core, 8 hour runtime, 8GB RAM, job array, 10 tasks run concurrently; tell nodes conda location
qsub -cwd -pe smp 1 -l h_rt=8:0:0,h_vmem=8G,tmem=8G -t 1-$TASK_COUNT -tc 10 -v GENOME=$GENOME,CHR_SUFFIX=$CHR_SUFFIX,RESULTS_BASE_DIR=$RESULTS_BASE_DIR,SAMPLE_METADATA=$SAMPLE_METADATA_TMP,PPLN_BASE_DIR=$PPLN_BASE_DIR,CONDA_EXE=$CONDA_EXE job.sh

