#!/bin/bash
# usage: bash submit_array.sh -g hg19 -r /SAN/colcc/alex_work/samples_analysis -c nochr -s configuration/2021-06_icgc_cancer_vs_normal.txt -m 16 -h 10 -S 400
#
# or specify a SGE_TASK_ID, if running an interactive shell, in which case
# usage: bash submit_array.sh -g hg19 -r /SAN/colcc/alex_work/samples_analysis -c nochr -s configuration/2021-06_icgc_cancer_vs_normal.txt -t 57 >stdout.txt 2>stderr.txt &

# some defaults
GENOME="hg19" # default; must be hg19 or grch38
CHR_SUFFIX="nochr" # default; must be either "nochr" OR "chr"; genome files (fasta, bed) are of the form grch38_nochr*, hg19_chr*, etc
RESULTS_BASE_DIR=/SAN/colcc/alex_work/samples_analysis
SAMPLE_METADATA=configuration/2021-06_icgc_cancer_vs_normal.txt
SGE_TASK_ID=-1
REQUEST_GB_MEM=16
REQUEST_HRS_RUN=16
REQUEST_GB_SCRATCH=-1 # request explicitly with -S if scratch space required

# set variables from arguments
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" # https://stackoverflow.com/questions/59895/how-can-i-get-the-source-directory-of-a-bash-script-from-within-the-script-itsel
PPLN_BASE_DIR=$SCRIPT_DIR/../..
while getopts g:c:r:s:m:h:S:t: flag
do
    case "${flag}" in
        g) GENOME=${OPTARG};;             # 'hg19' or 'grch38'
        c) CHR_SUFFIX=${OPTARG};;         # 'chr' or 'nochr' if chromosomes are specified with or without 'chr' respectively; used to choose appropriately named genomic bed and fasta files
        r) RESULTS_BASE_DIR=${OPTARG};;   # base directory for analysis results
        s) SAMPLE_METADATA=${OPTARG};;    # sample metadata file
        m) REQUEST_GB_MEM=${OPTARG};;     # request GB RAM per core
        h) REQUEST_HRS_RUN=${OPTARG};;    # request hours of runtime per job
        S) REQUEST_GB_SCRATCH=${OPTARG};; # request GB scratch space
        t) SGE_TASK_ID=${OPTARG};;        # run a specific task (line number in the file, after removal of successfully completed jobs)
    esac
done
echo "Genome: $GENOME"
echo "chr in file name: $CHR_SUFFIX"
echo "Results location: $RESULTS_BASE_DIR"
echo "Sample metadata: $SAMPLE_METADATA"
echo "Requested GB memory: $REQUEST_GB_MEM"
echo "Requested run hours: $REQUEST_HRS_RUN"
echo "Requested GB scratch: $REQUEST_GB_SCRATCH"
## NOTE: if underscore is placed at beginning of a text string, seems to need double quotes around it when assigning to a variable: echo works ok, but passing it to a script doesn't seem to, so avoid if possible


if [ ! -d $RESULTS_BASE_DIR ]
then
    mkdir $RESULTS_BASE_DIR
fi

# if a job is specified (only if running in interactive shell), analyse chosen sample and exit the script
if [ $SGE_TASK_ID -gt -1 ]
then
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

# scratch request string
if [ $REQUEST_GB_SCRATCH -gt 0 ]
then
    SCRATCH_RQ_STR=",tscratch=$REQUEST_GB_SCRATCH"G
else
    SCRATCH_RQ_STR=""
fi

echo "GENOME=$GENOME SAMPLE_METADATA=$SAMPLE_METADATA_TMP RESULTS_BASE_DIR=$RESULTS_BASE_DIR TASK_COUNT=$TASK_COUNT SGE_TASK_ID=$SGE_TASK_ID "


# qsub the array job
# Request 1 core, 10 hour runtime, 16GB RAM, job array, 10 tasks run concurrently; 400G scratch space enough for ICGC data; tell nodes conda location
# note '-pe smp N' is for multi-threaded jobs with N>1 cores and restricts jobs to multi-core hosts - unnecessary
qsub -cwd -l "h_rt="$REQUEST_HRS_RUN":0:0,h_vmem="$REQUEST_GB_MEM"G,tmem="$REQUEST_GB_MEM"G"$SCRATCH_RQ_STR -t 1-$TASK_COUNT -tc 10 -v GENOME=$GENOME,CHR_SUFFIX=$CHR_SUFFIX,RESULTS_BASE_DIR=$RESULTS_BASE_DIR,SAMPLE_METADATA=$SAMPLE_METADATA_TMP,PPLN_BASE_DIR=$PPLN_BASE_DIR,CONDA_EXE=$CONDA_EXE,REQUEST_GB_MEM=$REQUEST_GB_MEM job.sh

