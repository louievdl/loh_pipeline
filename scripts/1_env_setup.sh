
#################################################
#
# purpose: build environment(s) for running
#    sequenza-polysolver-LOHHLA pipeline
# run this script first, in interactive session
#
#################################################

# typical interactive login to CoLCC
qlogin -l h_rt=3:0:0,h_vmem=8G,tmem=8G

module load python # v3.8.5 default; only necessary during initial installation

# define root dir for this analysis
LOH_BASE_DIR="/SAN/colcc/alex_work"

# get pipeline setup scripts and other files
cd $LOH_BASE_DIR
git clone https://github.com/louievdl/hla_loh_pipeline.git

# setup anaconda
# note: conda 4.8.2 (vs a newer version) was required to achieve a working set of polysolver supporting packages. Newer versions ended in error
wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
bash Anaconda3-2020.02-Linux-x86_64.sh # installed in $LOH_BASE_DIR/anaconda3
conda config --set auto_activate_base false # don't activate conda base env on startup
conda config --add channels defaults # Warning: 'defaults' already in 'channels' list, moving to the top
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels compbiocore
conda config --add channels vacation
conda config --add channels r
exit # restart shell for changes to take effect

# env for polysolver
conda create -n polysolver --file $LOH_BASE_DIR/hla_loh_pipeline/polysolver-package-list.txt # success reinstalling packages from export file, see above

# env for lohhla
conda create -n lohhla python=3 # skip conda update
conda install -n lohhla lohhla=20171108 # success; installed here as a way of obtaining supporting packages; we use an up to date version in this environment
conda install -n lohhla jellyfish=2.2.6 samtools=1.11 r-tidyverse=1.3.1 gcc=5.4.0 r-naturalsort=0.1.3 libgomp=11.1.0
conda install -n lohhla -c bioconda --force-reinstall bioconductor-genomeinfodbdata=1.2.6 # updated version contains required data; results in large exchange of packages.
# in my hands, this call results in a shipload of package changes; I checked that this reinstall worked:
#conda activate lohhla
#R --vanilla
#> library(GenomeInfoDb) # unsuccessful, R complained "no package called ‘GenomeInfoDbData’"
# after rerunning the reinstall, package loaded ok

# after lohhla, packages were installed in this order...
#conda install -n lohhla jellyfish=2.2.6 # done 20210616
#conda install -n lohhla samtools=1.11  # done 20210616, to avoid error, see https://www.biostars.org/p/455593/, overwrite samtools 1.3, which came with lohhla
#conda install -n lohhla -c r r-tidyverse=1.3.1
#conda install -n lohhla gcc=5.4.0
#conda install -n lohhla -c conda-forge r-naturalsort=0.1.3
#conda install -n lohhla -c bioconda --force-reinstall bioconductor-genomeinfodbdata=1.2.6 # GenomeInfoDb and GenomeInfoDbData are listed as installed, but GenomeInfoDbData was not present
#conda update -n lohhla -c conda-forge libgomp=11.1.0

# env for sequenza
conda create -n sequenza
conda install -n sequenza sequenza-utils=3.0.0 r-sequenza=3.0.0 samtools=1.13 bedtools=2.30.0 bcftools=1.13 wget=1.20.1
# alternate piecemeal way to install these
#conda install -n sequenza -c bioconda sequenza-utils=3.0.0
#conda install -n sequenza -c bioconda r-sequenza=3.0.0 # installs old version of samtools requiring libcrypto.so.1.0.0, from openssl
##conda update -n sequenza samtools # downgrades samtools from 1.7.1 to 1.3.1 from compbiocore channel
#conda install -n sequenza samtools=1.13
#conda install -n sequenza bedtools=2.30.0
#conda install -n sequenza bcftools=1.13
#conda install -n sequenza -c anaconda wget=1.20.1

# no script updates necessary for sequenza - use as is

# get updated polysolver
cd $LOH_BASE_DIR
git clone https://github.com/louievdl/hla-polysolver.git # 779 MB, including dbs
cd $LOH_BASE_DIR/hla-polysolver/
SRC_DIR=$LOH_BASE_DIR/hla-polysolver && PREFIX=$LOH_BASE_DIR/polysolver_inst && mkdir $PREFIX
source build.sh # put tools in right places and let polysolver know where it is; finishes in $LOH_BASE_DIR/hla-polysolver/include/strelka-upstream-v1.0.11
cd $LOH_BASE_DIR
mkdir $LOH_BASE_DIR/hla_fasta
cp -a $LOH_BASE_DIR/hla-polysolver/data/abc_complete.fasta $LOH_BASE_DIR/hla_fasta
# novoalign_header.sam is generated from abc_complete.fasta


# get updated lohhla
cd $LOH_BASE_DIR
git clone https://github.com/louievdl/LOHHLA.git
cp -a $LOH_BASE_DIR/LOHHLA/data/hla.dat $LOH_BASE_DIR/hla_fasta

