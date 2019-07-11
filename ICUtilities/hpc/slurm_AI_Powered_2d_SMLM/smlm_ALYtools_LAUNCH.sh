#!/bin/bash
#
#  Created by Yuriy Alexandrov
#

SRCDIR=$1
DSTDIR=$2
settings_file=$3
src_filename_must_have_token=$4
file_extension_template=$5

OMEtiffcount=0
for f in $SRCDIR/*; do
	case "$f" in
	*$file_extension_template*)
	    #echo $(basename "$f")
	    OMEtiffcount=$((OMEtiffcount + 1))
	;;
	esac
done
echo $OMEtiffcount

if (( (($OMEtiffcount == 0)) )); then
      echo "nothing to do, exiting"
      return
fi

mrt_camp=matlab-runtime

ARGS="$SRCDIR:$DSTDIR:$OMEtiffcount:$mrt_camp:$settings_file:$src_filename_must_have_token:$file_extension_template"
echo $ARGS

one=$(sbatch $HOME/smlm_ALYtools_camp_software/smlm_ALYtools_Dispatch.slurm $ARGS)

echo "launching processing job"
echo $one
 
