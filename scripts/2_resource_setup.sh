#!/bin/bash

#################################################
#
# run this script second to set up genomic data
# supporting sequenza-polysolver-lohhla pipeline
#
#################################################

source ~/.bashrc # add condabin to PATH
echo "path $PATH"
LOH_BASE_DIR=/SAN/colcc/alex_work # TODO MOVE THIS

# get and process genomic fasta files if not present

if [ ! -d "$LOH_BASE_DIR/genome" ]
then
    mkdir $LOH_BASE_DIR/genome
fi

if [ ! -f $LOH_BASE_DIR/genome/grch38_nochr.fa ] || [ ! -f $LOH_BASE_DIR/genome/hg19_nochr.fa ]
then

    # download fasta files
    conda activate sequenza

    wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz | gunzip -c > $LOH_BASE_DIR/genome/grch38.fa
    wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz | gunzip -c > $LOH_BASE_DIR/genome/hg19.fa
    cat $LOH_BASE_DIR/genome/grch38.fa | perl -ne 's/^>chr/>/; print' > $LOH_BASE_DIR/genome/grch38_nochr.fa
    cat $LOH_BASE_DIR/genome/hg19.fa | perl -ne 's/^>chr/>/; print' > $LOH_BASE_DIR/genome/hg19_nochr.fa

    conda deactivate
fi

if [ ! -f $LOH_BASE_DIR/genome/hg19.fa.fai ] || [ ! -f $LOH_BASE_DIR/genome/hg19_nochr.fa.fai ] || [ ! -f $LOH_BASE_DIR/genome/grch38.fa.fai ] || [ ! -f $LOH_BASE_DIR/genome/grch38_nochr.fa.fai ]
then

    # index fasta files
    conda activate sequenza

    samtools faidx $LOH_BASE_DIR/genome/hg19.fa --fai-idx $LOH_BASE_DIR/genome/hg19.fa.fai
    samtools faidx $LOH_BASE_DIR/genome/hg19_nochr.fa --fai-idx $LOH_BASE_DIR/genome/hg19_nochr.fa.fai
    samtools faidx $LOH_BASE_DIR/genome/grch38.fa --fai-idx $LOH_BASE_DIR/genome/grch38.fa.fai
    samtools faidx $LOH_BASE_DIR/genome/grch38_nochr.fa --fai-idx $LOH_BASE_DIR/genome/grch38_nochr.fa.fai

    conda deactivate
fi

if [ ! -f $LOH_BASE_DIR/genome/grch38.gc50Base.wig.gz ] || [ ! -f $LOH_BASE_DIR/genome/hg19.gc50Base.wig.gz ] || [ ! -f $LOH_BASE_DIR/genome/grch38_nochr.gc50Base.wig.gz ] || [ ! -f $LOH_BASE_DIR/genome/hg19_nochr.gc50Base.wig.gz ]
then

    # Process FASTA files to produce GC Wiggle track files
    conda activate sequenza

    sequenza-utils gc_wiggle -w 50 --fasta $LOH_BASE_DIR/genome/grch38.fa -o $LOH_BASE_DIR/genome/grch38.gc50Base.wig.gz
    sequenza-utils gc_wiggle -w 50 --fasta $LOH_BASE_DIR/genome/hg19.fa -o $LOH_BASE_DIR/genome/hg19.gc50Base.wig.gz
    sequenza-utils gc_wiggle -w 50 --fasta $LOH_BASE_DIR/genome/grch38_nochr.fa -o $LOH_BASE_DIR/genome/grch38_nochr.gc50Base.wig.gz
    sequenza-utils gc_wiggle -w 50 --fasta $LOH_BASE_DIR/genome/hg19_nochr.fa -o $LOH_BASE_DIR/genome/hg19_nochr.gc50Base.wig.gz

    conda deactivate
fi

# get protein-coding exon annotation if not present; needed for sequenza copy number analysis, to be done using exons of protein coding genes

if [ ! -d "$LOH_BASE_DIR/copy_number_analysis" ]
then
    mkdir $LOH_BASE_DIR/copy_number_analysis
fi

if [ ! -d "$LOH_BASE_DIR/copy_number_analysis/bedfiles" ]
then
    mkdir $LOH_BASE_DIR/copy_number_analysis/bedfiles
fi

if [ ! -f $LOH_BASE_DIR/copy_number_analysis/bedfiles/grch38_exons_hlaregions.bed ] || [ ! -f $LOH_BASE_DIR/copy_number_analysis/bedfiles/hg19_exons_hlaregions.bed ]
then

    conda activate sequenza

    echo "6	29909037	29913661" > $LOH_BASE_DIR/copy_number_analysis/bedfiles/hg19_nochr_hlaregions.bed
    echo "6	31321649	31324964" >> $LOH_BASE_DIR/copy_number_analysis/bedfiles/hg19_nochr_hlaregions.bed
    echo "6	31236526	31239869" >> $LOH_BASE_DIR/copy_number_analysis/bedfiles/hg19_nochr_hlaregions.bed

    echo "6	29941260	29945884" > $LOH_BASE_DIR/copy_number_analysis/bedfiles/grch38_nochr_hlaregions.bed
    echo "6	31353872	31357187" >> $LOH_BASE_DIR/copy_number_analysis/bedfiles/grch38_nochr_hlaregions.bed
    echo "6	31268749	31272092" >> $LOH_BASE_DIR/copy_number_analysis/bedfiles/grch38_nochr_hlaregions.bed

    wget -O - "http://ftp.ensembl.org/pub/grch37/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz" |
    gunzip -c | egrep '\sexon\s' | grep gene_biotype..protein_coding | 
    perl -F"\t" -ane '/gene_name."([\w\.\_]+)"/; $gene = $1; print "$F[0]\t$F[3]\t$F[4]\t$gene\t1\t$F[6]\n"' | 
    sort -k 1,1 -k2,2n | bedtools merge > $LOH_BASE_DIR/copy_number_analysis/bedfiles/hg19_nochr_exons.bed

    wget -O - "http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz" |
    gunzip -c | egrep '\sexon\s' | grep gene_biotype..protein_coding | 
    perl -F"\t" -ane '/gene_name."([\w\.\_]+)"/; $gene = $1; print "$F[0]\t$F[3]\t$F[4]\t$gene\t1\t$F[6]\n"' | 
    sort -k 1,1 -k2,2n | bedtools merge > $LOH_BASE_DIR/copy_number_analysis/bedfiles/grch38_nochr_exons.bed

    cat $LOH_BASE_DIR/copy_number_analysis/bedfiles/hg19_nochr_exons.bed $LOH_BASE_DIR/copy_number_analysis/bedfiles/hg19_nochr_hlaregions.bed | 
    sort -k 1,1 -k2,2n | bedtools merge > $LOH_BASE_DIR/copy_number_analysis/bedfiles/hg19_nochr_exons_hlaregions.bed

    cat $LOH_BASE_DIR/copy_number_analysis/bedfiles/grch38_nochr_exons.bed $LOH_BASE_DIR/copy_number_analysis/bedfiles/grch38_nochr_hlaregions.bed | 
    sort -k 1,1 -k2,2n | bedtools merge > $LOH_BASE_DIR/copy_number_analysis/bedfiles/grch38_nochr_exons_hlaregions.bed

    cat $LOH_BASE_DIR/copy_number_analysis/bedfiles/hg19_nochr_exons_hlaregions.bed | perl -ne 'print "chr$_"' > $LOH_BASE_DIR/copy_number_analysis/bedfiles/hg19_exons_hlaregions.bed
    cat $LOH_BASE_DIR/copy_number_analysis/bedfiles/grch38_nochr_exons_hlaregions.bed | perl -ne 'print "chr$_"' > $LOH_BASE_DIR/copy_number_analysis/bedfiles/grch38_exons_hlaregions.bed

    conda deactivate
fi


