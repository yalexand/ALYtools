#usage: qsub test.sh

#!/bin/bash
#PBS -lselect=1:ncpus=1:mem=2gb
#PBS -lwalltime=0:10:0

#PBS -J 1-2

module load java
module load julia

JULIA_DIR=/work/yalexand/OPT_JULIA/JULIA

SRC_DIR=/work/yalexand/TestData
DST_DIR=/work/yalexand

fname=fluor.OME.tiff

fname=$(head -n $PBS_ARRAY_INDEX $SRC_DIR/filelist.txt | tail -1 )
 
cd $JULIA_DIR

src=$SRC_DIR/$fname
dst=$DST_DIR/reconstr_$fname
julia OPT_reconstruction_test_nix.jl $src $dst


