#!/bin/bash
#
#  Created by Yuriy Alexandrov on 26/02/2018.
#

USAGE="Usage: LAUNCH_OPT_Reconstruction src_path dst_path swap_flag"

echo "fogim queue"
QUEUE="pqfogim"

SRCDIR=$1
DSTDIR=$2
SWAP_FLAG=$3

dircount=0
OMEtiffcount=0
dircount=0
for f in $SRCDIR/*; do
	case "$f" in
	*.OME.tiff*)
	    #echo $(basename "$f")
	    OMEtiffcount=$((OMEtiffcount + 1))
	;;
	esac
if [ -d  $f ]; then
  echo "Directory -> $f"
  dircount=$((dircount + 1))    
fi
done
echo $OMEtiffcount
echo $dircount


if (( (($dircount == 0)) && (($OMEtiffcount == 0)) )); then
      echo "nothing to do, exiting"
      return
fi

if (( (($dircount != 0)) && (($OMEtiffcount != 0)) )); then
      echo "confusing source specification, exiting"
      return
fi

ARGS="$SRCDIR":"$DSTDIR":"$SWAP_FLAG":"$OMEtiffcount":"$dircount"
echo $ARGS

# WITHOUT GPU
one=$(qsub -q $QUEUE -v SETUP_ARGS=$ARGS $HOME/OPT_Reconstruction/OPT_ARRScript.pbs)

# TRY IT WITH GPU
# one=$(qsub -q $QUEUE -v SETUP_ARGS=$ARGS $HOME/OPT_Reconstruction/OPT_ARRScript_GPU.pbs)

echo "launching processing job"
echo $one






