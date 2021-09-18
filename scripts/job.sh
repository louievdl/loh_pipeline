#!/bin/bash
# note, all job attributes specified outside this script

# submission script supplies GENOME, CHR_SUFFIX, RESULTS_BASE_DIR, SAMPLE_METADATA

source ~/.bashrc # add condabin to PATH
echo "path $PATH"



#################################################
#
# CHOOSE SAMPLE AND CREATE LINKS TO BAMS
#
#################################################

SAMPLE_IDS=($(cat $SAMPLE_METADATA | cut -f 1))
CTL_FILES=($(cat $SAMPLE_METADATA | cut -f 2))
TRT_FILES=($(cat $SAMPLE_METADATA | cut -f 3))
SAMPLE_ID=${SAMPLE_IDS[SGE_TASK_ID]}
CTL_BAM=${CTL_FILES[SGE_TASK_ID]}
TRT_BAM=${TRT_FILES[SGE_TASK_ID]}
SAMPLE_DIR="$RESULTS_BASE_DIR/$SAMPLE_ID" # compute in individual task; should be a full path to an individual's files with sanitised names
mkdir $SAMPLE_DIR
BAM_DIR="$SAMPLE_DIR/bams"
mkdir $BAM_DIR
CTL_BAM_LINK=$BAM_DIR/$(basename $CTL_BAM | perl -ne 's/\./_/g; s/_bam$/.bam/; print') # sanitise bam names and link sanitised names to orig files
TRT_BAM_LINK=$BAM_DIR/$(basename $TRT_BAM | perl -ne 's/\./_/g; s/_bam$/.bam/; print')
ln -s $CTL_BAM $CTL_BAM_LINK
ln -s $TRT_BAM $TRT_BAM_LINK
ln -s $CTL_BAM.bai $CTL_BAM_LINK.bai
ln -s $TRT_BAM.bai $TRT_BAM_LINK.bai
CTL_BASENAME=$(basename $CTL_BAM_LINK | perl -ne 's/\.bam$//; print') # remove .bam
TRT_BASENAME=$(basename $TRT_BAM_LINK | perl -ne 's/\.bam$//; print')



#################################################
#
# SEQUENZA
#
#################################################

### SET SEQUENZA VARS AND DIRS ###

SEQUENZA_DIR=$SAMPLE_DIR/sequenza
if [ ! -d $SEQUENZA_DIR ]
then
    mkdir $SEQUENZA_DIR
fi


### CREATE PILEUPS, BINNED SEQZ FILES ###

conda activate sequenza

if [ ! -f $SEQUENZA_DIR/ctl.pileup ]
then
    samtools mpileup --fasta-ref $PPLN_BASE_DIR/genome/$GENOME"$CHR_SUFFIX".fa \
        --positions $PPLN_BASE_DIR/genome/$GENOME"$CHR_SUFFIX"_exons.bed \
        --min-MQ 20 --excl-flags 256 $CTL_BAM \
        -o $SEQUENZA_DIR/ctl.pileup # ~68GB bam took 1h12m2s clock time to run, producing file of size 8,818,616,226; <1GB RAM
fi

if [ ! -f $SEQUENZA_DIR/trt.pileup ]
then
    samtools mpileup --fasta-ref $PPLN_BASE_DIR/genome/$GENOME"$CHR_SUFFIX".fa \
        --positions $PPLN_BASE_DIR/genome/$GENOME"$CHR_SUFFIX"_exons.bed \
        --min-MQ 20 --excl-flags 256 $TRT_BAM \
        -o $SEQUENZA_DIR/trt.pileup # ~68GB bam took 1h4m18s clock time to run, producing file of size 8,342,805,140; would this be faster if connections were faster?
fi

if [ ! -f $SEQUENZA_DIR/binned.seqz.gz ]
then
    sequenza-utils bam2seqz -p -n $SEQUENZA_DIR/ctl.pileup -t $SEQUENZA_DIR/trt.pileup \
        --fasta $PPLN_BASE_DIR/genome/$GENOME"$CHR_SUFFIX".fa \
        -gc $PPLN_BASE_DIR/genome/$GENOME"$CHR_SUFFIX".gc50Base.wig.gz | \
        sequenza-utils seqz_binning -s - -w 100000 -o $SEQUENZA_DIR/binned.seqz.gz # pileups from two 68GB bams took 29m9s clock time to run, producing file of size 17,025,764
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
    ' $SEQUENZA_DIR/binned.seqz.gz $TRT_BASENAME > $SEQUENZA_DIR/purity_ploidy.txt
    #   Ploidy  tumorPurity     tumorPloidy     
    #example_tumor_sorted    2       0.8     1.8     # not sure what Ploidy is for, need tumor* vars for lohhla
fi



#################################################
#
# POLYSOLVER
#
#################################################

### SET POLYSOLVER VARS AND DIRS ###

POLYSOLVER_DIR=$SAMPLE_DIR/polysolver
if [ ! -d $POLYSOLVER_DIR ]
then
    mkdir $POLYSOLVER_DIR
fi

POLYSOLVER_WORK_DIR=$POLYSOLVER_DIR/work
if [ ! -d "$POLYSOLVER_WORK_DIR" ]
then
    mkdir $POLYSOLVER_WORK_DIR
fi

POLYSOLVER_OUT_DIR=$POLYSOLVER_DIR/results
if [ ! -d $POLYSOLVER_OUT_DIR ]
then
    mkdir $POLYSOLVER_OUT_DIR
fi

#echo "running pipeline with genome $GENOME; files-list $SAMPLE_METADATA; SGE_TASK_ID $SGE_TASK_ID; pt-id $SAMPLE_ID; control-bam $CTL_BAM; treated-bam $TRT_BAM; bam-dir $BAM_DIR; target-bam $BAM; work-dir $WORK_DIR; out-dir $OUT_DIR"


### RUN POLYSOLVER ###

if [ ! -f $POLYSOLVER_OUT_DIR/winners.hla.txt ]
then

    conda activate polysolver
    $PPLN_BASE_DIR/polysolver_inst/bin/shell_call_hla_type $CTL_BAM Unknown 1 $GENOME STDFQ 0 $POLYSOLVER_WORK_DIR # slightly modified polysolver script, permitting custom HLA allele library
    conda deactivate

    if [ -f $POLYSOLVER_WORK_DIR/winners.hla.txt ]
    then
        mv $POLYSOLVER_WORK_DIR/winners.hla.txt $POLYSOLVER_OUT_DIR
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
if [ ! -d $LOHHLA_DIR ]
then
    mkdir $LOHHLA_DIR
fi

LOHHLA_WORK_DIR=$LOHHLA_DIR/work
if [ ! -d $LOHHLA_DIR/work ]
then
    mkdir $LOHHLA_DIR/work
fi

LOHHLA_OUT_DIR=$LOHHLA_DIR/results
if [ ! -d $LOHHLA_DIR/results ]
then
    mkdir $LOHHLA_DIR/results
fi
#echo "running lohhla with genome $GENOME; pt-id $SAMPLE_ID; out-dir $OUT_DIR"


### RUN LOHHLA ###

# use less than total available memory for sorting
SORT_MEM=16G

#source ~/.bashrc # put conda into path # already done
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

# lohhla is fork of slagtermaarten/LOHHLA
#echo "$SHELL"
Rscript $PPLN_BASE_DIR/LOHHLA/LOHHLAscript.R \
  --LOHHLA_loc=$PPLN_BASE_DIR/LOHHLA \
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
  --HLAexonLoc=$PPLN_BASE_DIR/hla_fasta/hla.dat & # flat file showing exons; note, hla_x_nn names in this file are in standard notation
# use a copynumloc file, not FALSE.
# by supplying estimated read counts, cut down run time substantially

conda deactivate

