#!/bin/bash
#SAMPLE_METADATA=configuration/xxx.txt
#TASK_COUNT=$((`grep -c ^ $SAMPLE_METADATA` - 1))
#qsub -cwd -l h_rt=0:10:0,h_vmem=16G,tmem=16G -t 1-$TASK_COUNT -tc 10 -v GENOME=hg19,CHR_SUFFIX=nochr,RESULTS_BASE_DIR=/SAN/colcc/alex_work/samples_analysis,SAMPLE_METADATA=$SAMPLE_METADATA,PPLN_BASE_DIR=/SAN/colcc/alex_work/pipeline,CONDA_EXE=$CONDA_EXE redo_lohhla.sh
# get TASK_COUNT from num lines in configuration/xxx.txt

# running env supplies GENOME, CHR_SUFFIX, RESULTS_BASE_DIR, SAMPLE_METADATA
# these jobs may be given only 10 min or so?

echo "starting job $JOB_ID.$SGE_TASK_ID at "`date`

#source ~/.bashrc # add condabin to PATH
echo "path $PATH ; conda $CONDA_EXE"

# in case this is a subshell, set anaconda vars
source $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh




#################################################
#
# CHOOSE SAMPLE AND CREATE LINKS TO LARGE FILES
#
#################################################

SAMPLE_IDS=($(cat $SAMPLE_METADATA | cut -f 1))
CTL_FILES=($(cat $SAMPLE_METADATA | cut -f 2))
TRT_FILES=($(cat $SAMPLE_METADATA | cut -f 3))
SAMPLE_ID=${SAMPLE_IDS[SGE_TASK_ID]}
CTL_BAM=${CTL_FILES[SGE_TASK_ID]}
TRT_BAM=${TRT_FILES[SGE_TASK_ID]}
SAMPLE_DIR="$RESULTS_BASE_DIR/$SAMPLE_ID" # compute in individual task; should be a full path to an individual's files with sanitised names
BAM_DIR="$SAMPLE_DIR/bams"
CTL_BAM_LINK=$BAM_DIR/$(basename $CTL_BAM | perl -ne 's/\./_/g; s/_bam$/.bam/; print') # sanitise bam names and link sanitised names to orig files, either in scratch or orig location
TRT_BAM_LINK=$BAM_DIR/$(basename $TRT_BAM | perl -ne 's/\./_/g; s/_bam$/.bam/; print')
CTL_BASENAME=$(basename $CTL_BAM_LINK | perl -ne 's/\.bam$//; print') # remove .bam
TRT_BASENAME=$(basename $TRT_BAM_LINK | perl -ne 's/\.bam$//; print')

[ ! -d $SAMPLE_DIR ] && mkdir $SAMPLE_DIR
[ ! -d $BAM_DIR    ] && mkdir $BAM_DIR

# remove existing links which may be misconnected
[ -L $CTL_BAM_LINK     ] && rm -f $CTL_BAM_LINK
[ -L $CTL_BAM_LINK.bai ] && rm -f $CTL_BAM_LINK.bai
[ -L $TRT_BAM_LINK     ] && rm -f $TRT_BAM_LINK
[ -L $TRT_BAM_LINK.bai ] && rm -f $TRT_BAM_LINK.bai

# use files directly via symlink - forbidden on UCL HPC
ln -s $CTL_BAM $CTL_BAM_LINK
ln -s $TRT_BAM $TRT_BAM_LINK
ln -s $CTL_BAM.bai $CTL_BAM_LINK.bai
ln -s $TRT_BAM.bai $TRT_BAM_LINK.bai



#################################################
#
# SEQUENZA COPY NUMBER ANALYSIS
#
#################################################

### SET SEQUENZA VARS AND DIRS ###

SEQUENZA_DIR=$SAMPLE_DIR/sequenza



#################################################
#
# POLYSOLVER
#
#################################################


### SET POLYSOLVER VARS AND DIRS ###

POLYSOLVER_DIR=$SAMPLE_DIR/polysolver
POLYSOLVER_OUT_DIR=$POLYSOLVER_DIR/results

# vars assumed accessible to polysolver script
HLA_DATA_DIR=$PPLN_BASE_DIR/hla_fasta
HLA_FASTA=abc_complete.fasta

POLYSOLVER_WORK_DIR=$POLYSOLVER_DIR/work



#################################################
#
# LOHHLA
#
#################################################

echo "starting LOHHLA analysis at "`date`

### SET LOHHLA VARS AND DIRS ###

CTL_BAM_SIZE=$(stat -c%s "$CTL_BAM") # rough estimates of num reads
CTL_BAM_READ_COUNT=$(($CTL_BAM_SIZE / 80))
TRT_BAM_SIZE=$(stat -c%s "$TRT_BAM")
TRT_BAM_READ_COUNT=$(($TRT_BAM_SIZE / 80))

LOHHLA_DIR=$SAMPLE_DIR/lohhla
LOHHLA_OUT_DIR=$LOHHLA_DIR/results


### RUN LOHHLA ###

conda activate lohhla

# LOHHLA is fork of bitbucket.org/mcgranahanlab/lohhla
#echo "$SHELL"
Rscript $PPLN_BASE_DIR/LOHHLA/LOHHLAscript.R \
  --genomeAssembly=$GENOME \
  --patientId=$SAMPLE_ID \
  --outputDir=$LOHHLA_OUT_DIR \
  --BAMDir=$BAM_DIR \
  --normalBAMfile=$CTL_BAM_LINK \
  --normalAlignedReads=$CTL_BAM_READ_COUNT \
  --tumorBAMfile=$TRT_BAM_LINK \
  --tumorAlignedReads=$TRT_BAM_READ_COUNT \
  --minCoverageFilter=15 \
  --hlaPath=$POLYSOLVER_OUT_DIR/winners.hla.txt \
  --HLAfastaLoc=$PPLN_BASE_DIR/hla_fasta/abc_complete.fasta \
  --CopyNumLoc=$SEQUENZA_DIR/purity_ploidy.txt \
  --plottingStep=TRUE \
  --cleanUp=TRUE \
  --HLAexonLoc=$PPLN_BASE_DIR/hla_fasta/hla.dat
#  || touch $LOHHLA_OUT_DIR/Figures/no_result.pdf # don't do this; if LOHHLA runs, worst case it puts out an empty pdf
# hla.dat is flat file showing exons; note, hla_x_nn names in this file are in standard notation
# use a copynumloc file, not FALSE.
# by supplying estimated read counts, cut down run time substantially
#--LOHHLA_loc=$PPLN_BASE_DIR/LOHHLA \ # forking from mcg, this opt unnec
[ `ls $LOHHLA_OUT_DIR/Figures/*.pdf | grep -c ^` -eq 0 ] && touch $LOHHLA_OUT_DIR/Figures/no_lohhla_result.pdf

conda deactivate
echo "finished LOHHLA analysis at "`date`
