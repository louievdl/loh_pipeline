#!/bin/bash
# note, all job attributes specified outside this script

# submission script supplies GENOME, CHR_SUFFIX, RESULTS_BASE_DIR, SAMPLE_METADATA

#source ~/.bashrc # add condabin to PATH
echo "path $PATH ; conda $CONDA_EXE"

# in case this is a subshell, set anaconda vars
source $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh

# deletion of scratch dir at task end, regardless of condition
SCRATCH_DIR=/scratch0/lvandela/$JOB_ID.$SGE_TASK_ID
function finish {
    echo "scratch contents..."
    ls -l $SCRATCH_DIR
    echo "deleting scratch..."
    rm -rf $SCRATCH_DIR
}
# Always wipe scratch
trap finish EXIT ERR INT TERM




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

USE_SCRATCH=true
if [ $USE_SCRATCH ]
then
    # use scratch space
    echo "mkdir -p $SCRATCH_DIR"
    mkdir -p $SCRATCH_DIR
    time rsync -av $CTL_BAM $TRT_BAM $CTL_BAM.bai $TRT_BAM.bai $SCRATCH_DIR/
    ln -s $SCRATCH_DIR/$(basename $CTL_BAM) $CTL_BAM_LINK
    ln -s $SCRATCH_DIR/$(basename $TRT_BAM) $TRT_BAM_LINK
    ln -s $SCRATCH_DIR/$(basename $CTL_BAM).bai $CTL_BAM_LINK.bai
    ln -s $SCRATCH_DIR/$(basename $TRT_BAM).bai $TRT_BAM_LINK.bai
    CTL_BAM=$SCRATCH_DIR/$(basename $CTL_BAM)
    TRT_BAM=$SCRATCH_DIR/$(basename $TRT_BAM)
    echo "scratch contents"
    ls -l $SCRATCH_DIR
else
    # use files directly via symlink - forbidden on UCL HPC
    ln -s $CTL_BAM $CTL_BAM_LINK
    ln -s $TRT_BAM $TRT_BAM_LINK
    ln -s $CTL_BAM.bai $CTL_BAM_LINK.bai
    ln -s $TRT_BAM.bai $TRT_BAM_LINK.bai
fi



#################################################
#
# SEQUENZA COPY NUMBER ANALYSIS
#
#################################################

### SET SEQUENZA VARS AND DIRS ###

SEQUENZA_DIR=$SAMPLE_DIR/sequenza
[ ! -d $SEQUENZA_DIR ] && mkdir $SEQUENZA_DIR


### CREATE PILEUPS, BINNED SEQZ FILES ###

conda activate sequenza

if [ $USE_SCRATCH ]
then
    CTL_PILEUP=$SCRATCH_DIR/ctl.pileup
    TRT_PILEUP=$SCRATCH_DIR/trt.pileup
    SEQZ_FILE=$SCRATCH_DIR/binned.seqz.gz
    [ -f $SEQUENZA_DIR/ctl.pileup     ] && time rsync -av $SEQUENZA_DIR/ctl.pileup     $CTL_PILEUP
    [ -f $SEQUENZA_DIR/trt.pileup     ] && time rsync -av $SEQUENZA_DIR/trt.pileup     $TRT_PILEUP
    [ -f $SEQUENZA_DIR/binned.seqz.gz ] && time rsync -av $SEQUENZA_DIR/binned.seqz.gz $SEQZ_FILE
else
    CTL_PILEUP=$SEQUENZA_DIR/ctl.pileup
    TRT_PILEUP=$SEQUENZA_DIR/trt.pileup
    SEQZ_FILE=$SEQUENZA_DIR/binned.seqz.gz
fi

if [ ! -f $SEQUENZA_DIR/ctl.pileup ]
then

    samtools mpileup --fasta-ref $PPLN_BASE_DIR/genome/$GENOME"_"$CHR_SUFFIX".fa" \
        --positions $PPLN_BASE_DIR/genome/$GENOME"_"$CHR_SUFFIX"_exons.bed" \
        --min-MQ 20 --excl-flags 256 $CTL_BAM \
        -o $CTL_PILEUP # ~68GB bam took 1h12m2s clock time to run, producing file of size 8,818,616,226; <1GB RAM

    [ $USE_SCRATCH ] && [ -f $CTL_PILEUP ] && rsync -av $CTL_PILEUP $SEQUENZA_DIR/ctl.pileup

fi

if [ ! -f $SEQUENZA_DIR/trt.pileup ]
then

    samtools mpileup --fasta-ref $PPLN_BASE_DIR/genome/$GENOME"_"$CHR_SUFFIX".fa" \
        --positions $PPLN_BASE_DIR/genome/$GENOME"_"$CHR_SUFFIX"_exons.bed" \
        --min-MQ 20 --excl-flags 256 $TRT_BAM \
        -o $TRT_PILEUP # ~68GB bam took 1h4m18s clock time to run, producing file of size 8,342,805,140; would this be faster if connections were faster?

    [ $USE_SCRATCH ] && [ -f $TRT_PILEUP ] && rsync -av $TRT_PILEUP $SEQUENZA_DIR/trt.pileup

fi

if [ ! -f $SEQUENZA_DIR/binned.seqz.gz ]
then
    sequenza-utils bam2seqz -p -n $CTL_PILEUP -t $TRT_PILEUP \
        --fasta $PPLN_BASE_DIR/genome/$GENOME"_"$CHR_SUFFIX".fa" \
        -gc $PPLN_BASE_DIR/genome/$GENOME"_"$CHR_SUFFIX".gc50Base.wig.gz" | \
        sequenza-utils seqz_binning -s - -w 100000 -o $SEQZ_FILE # pileups from two 68GB bams took 29m9s clock time to run, producing file of size 17,025,764

    [ $USE_SCRATCH ] && [ -f $SEQZ_FILE ] && rsync -av $SEQZ_FILE $SEQUENZA_DIR/binned.seqz.gz
fi


### ANALYSE SEQZ ###

if [ ! -f $SEQUENZA_DIR/purity_ploidy.txt ]
then
    Rscript -e '
        Args <- commandArgs(T);
        data.file <- Args[1];
        sample_name <- Args[2];
        library(sequenza);
        seqz.data <- read.seqz(data.file);
        seqz.data <- seqz.data[ ! seqz.data$chromosome %in% c("M", "MT"), ]
        gc.stats <- gc.sample.stats(data.file);
        gc.normal.vect <- mean_gc(gc.stats$normal);
        gc.tumor.vect <- mean_gc(gc.stats$tumor);
        seqz.data$adjusted.ratio <- round(( seqz.data$depth.tumor /  gc.tumor.vect[as.character(seqz.data$GC.percent)]) /
                                          (seqz.data$depth.normal / gc.normal.vect[as.character(seqz.data$GC.percent)]), 3);
        seqz.hom <- seqz.data$zygosity.normal == "hom";
        seqz.het <- seqz.data[!seqz.hom, ];
        breaks <- find.breaks(seqz.het, gamma = 80, kmin = 10, baf.thres = c(0, 0.5));
        seqz.data <- seqz.data[which(seqz.data$chromosome %in% breaks$chrom),];
        seg.s1 <- segment.breaks(seqz.data, breaks = breaks);
        seg.filtered <- seg.s1[(seg.s1$end.pos - seg.s1$start.pos) > 3e6, ];
        weights.seg  <- (seg.filtered$end.pos - seg.filtered$start.pos) / 1e6;
        avg.depth.ratio <- 1;
        CP <- baf.model.fit(Bf = seg.filtered$Bf,
                            depth.ratio = seg.filtered$depth.ratio,
                            weight.ratio = weights.seg,
                            weight.Bf = weights.seg,
                            sd.ratio = seg.filtered$sd.ratio,
                            sd.Bf = seg.filtered$sd.BAF,
                            avg.depth.ratio = avg.depth.ratio,
                            cellularity = seq(0.02, 1, 0.01),
                            ploidy = seq(0.5, 3, 0.05));
        confint <- get.ci(CP);
        ploidy <- confint$max.ploidy;
        cellularity <- confint$max.cellularity;
        cat("\tPloidy\ttumorPurity\ttumorPloidy\n",sample_name,"\t2\t",cellularity,"\t",ploidy,"\n",sep="");
    ' $SEQZ_FILE $TRT_BASENAME > $SEQUENZA_DIR/purity_ploidy.txt || true
    #   Ploidy  tumorPurity     tumorPloidy     
    #example_tumor_sorted    2       0.8     1.8     # not sure what Ploidy is for, need tumor* vars for lohhla
fi

conda deactivate



#################################################
#
# POLYSOLVER
#
#################################################

### SET POLYSOLVER VARS AND DIRS ###

POLYSOLVER_DIR=$SAMPLE_DIR/polysolver
POLYSOLVER_OUT_DIR=$POLYSOLVER_DIR/results
[ ! -d $POLYSOLVER_DIR ] && mkdir $POLYSOLVER_DIR
[ ! -d $POLYSOLVER_OUT_DIR ] && mkdir $POLYSOLVER_OUT_DIR

# vars assumed accessible to polysolver script
HLA_DATA_DIR=$PPLN_BASE_DIR/hla_fasta
HLA_FASTA=abc_complete.fasta

if [ $USE_SCRATCH ]
then
    POLYSOLVER_WORK_DIR=$SCRATCH_DIR/polysolver
else
    POLYSOLVER_WORK_DIR=$POLYSOLVER_DIR/work
fi
[ ! -d $POLYSOLVER_WORK_DIR ] && mkdir $POLYSOLVER_WORK_DIR

#echo "running pipeline with genome $GENOME; files-list $SAMPLE_METADATA; SGE_TASK_ID $SGE_TASK_ID; pt-id $SAMPLE_ID; control-bam $CTL_BAM; treated-bam $TRT_BAM; bam-dir $BAM_DIR; target-bam $BAM; work-dir $WORK_DIR; out-dir $OUT_DIR"


### RUN POLYSOLVER ###

if [ ! -f $POLYSOLVER_OUT_DIR/winners.hla.txt ]
then

    conda activate polysolver
    . $PPLN_BASE_DIR/polysolver_inst/bin/shell_call_hla_type $CTL_BAM Unknown 1 $GENOME STDFQ 0 $POLYSOLVER_WORK_DIR || true # slightly modified polysolver script, permitting custom HLA allele library
    conda deactivate

    if [ -f $POLYSOLVER_WORK_DIR/winners.hla.txt ]
    then
        rsync -av $POLYSOLVER_WORK_DIR/winners.hla.txt $POLYSOLVER_OUT_DIR/winners.hla.txt
        rm -r $POLYSOLVER_WORK_DIR/*
    fi
fi



#################################################
#
# LOHHLA
#
#################################################

### SET LOHHLA VARS AND DIRS ###

CTL_BAM_SIZE=$(stat -c%s "$CTL_BAM") # rough estimates of num reads
CTL_BAM_READ_COUNT=$(($CTL_BAM_SIZE / 80))
TRT_BAM_SIZE=$(stat -c%s "$TRT_BAM")
TRT_BAM_READ_COUNT=$(($TRT_BAM_SIZE / 80))

LOHHLA_DIR=$SAMPLE_DIR/lohhla
LOHHLA_OUT_DIR=$LOHHLA_DIR/results
[ ! -d $LOHHLA_DIR ] && mkdir $LOHHLA_DIR
[ ! -d $LOHHLA_DIR/results ] && mkdir $LOHHLA_DIR/results
[ ! -d $LOHHLA_DIR/results/Figures ] && mkdir $LOHHLA_DIR/results/Figures


### RUN LOHHLA ###

# use less than total available memory for sorting
SORT_MEM="$(($REQUEST_GB_MEM - 1))G"

conda activate lohhla # required for novoalign, samtools, bedtools, jellyfish

# sort input files if nec. use less than total requested memory.

TRT_BAM_SORTED=$(samtools view -H $TRT_BAM | egrep -c ^.HD..+SO.coordinate)
if [ $TRT_BAM_SORTED -eq 0 ]; then
    TRT_SORTED_BAM=$BAM_DIR/$(basename $TRT_BAM)
    samtools sort -m $SORT_MEM -T $BAM_DIR/tmp -O BAM $TRT_BAM > $TRT_SORTED_BAM
    samtools index $TRT_SORTED_BAM
    rm $TRT_BAM_LINK $TRT_BAM_LINK.bai
    ln -s $TRT_SORTED_BAM $TRT_BAM_LINK
    ln -s $TRT_SORTED_BAM.bai $TRT_BAM_LINK.bai
fi

CTL_BAM_SORTED=$(samtools view -H $CTL_BAM | egrep -c ^.HD..+SO.coordinate)
if [ $CTL_BAM_SORTED -eq 0 ]; then
    CTL_SORTED_BAM=$BAM_DIR/$(basename $CTL_BAM)
    samtools sort -m $SORT_MEM -T $BAM_DIR/tmp -O BAM $CTL_BAM > $CTL_SORTED_BAM
    samtools index $CTL_SORTED_BAM
    rm $CTL_BAM_LINK $CTL_BAM_LINK.bai
    ln -s $CTL_SORTED_BAM $CTL_BAM_LINK
    ln -s $CTL_SORTED_BAM.bai $CTL_BAM_LINK.bai
fi

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
  --HLAexonLoc=$PPLN_BASE_DIR/hla_fasta/hla.dat \
  || touch $LOHHLA_OUT_DIR/Figures/no_result.pdf
# hla.dat is flat file showing exons; note, hla_x_nn names in this file are in standard notation
# use a copynumloc file, not FALSE.
# by supplying estimated read counts, cut down run time substantially
#--LOHHLA_loc=$PPLN_BASE_DIR/LOHHLA \ # forking from mcg, this opt unnec

conda deactivate

