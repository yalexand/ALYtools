#!/bin/sh

#  NodeScript.sh
#  
#
#  Created by Yuriy Alexandrov 26/02/2018.
#  The script that runs on each node

echo "Start time $(date)"

ARGS="$1"

ARRARGS=(${ARGS//:/ })
SRCDIR=${ARRARGS[0]}
DSTDIR=${ARRARGS[1]}
SWAPFLAG=${ARRARGS[2]}
OMEtiffcount=${ARRARGS[3]}
dircount=${ARRARGS[4]}

echo $SRCDIR
echo $DSTDIR
echo $SWAPFLAG
echo $OMEtiffcount
echo $dircount
echo $PBS_ARRAY_INDEX

fname="blank"

#
# CASE WHEN INPUTS ARE XY STACKS
#
if (($dircount != 0)); then
cnt=0
for f in $SRCDIR/*; do
	if [ -d  $f ]; then
	  #echo "Directory -> $f"
	  cnt=$((cnt + 1))    
		if (($PBS_ARRAY_INDEX==$cnt)); then
			fname=$f
		fi
	fi	
done
fi

#
# CASE WHEN INPUTS ARE OME.tiff files
#
if (($OMEtiffcount != 0)); then
cnt=0
for f in $SRCDIR/*; do
	case "$f" in
	*.OME.tiff*)
	    #echo $(basename "$f")
	    cnt=$((cnt + 1))
		if (($PBS_ARRAY_INDEX==$cnt)); then
			fname=$f
		fi
	;;
	esac
done
fi

echo $fname


if [ "$fname" != "blank" ]; then

  space=" " 
  ARGS_FULL="/apps/MATLAB/R2017b"$space"$fname"$space"$DSTDIR"$space"$SWAPFLAG" 
  # ARGS_FULL="/apps/MATLAB/R2017b"$space"$fname"$space"$DSTDIR"$space"$SWAPFLAG" 
  echo $ARGS_FULL

  # module load sysconfcpus/0.5
  # sysconfcpus -n 24 $HOME/OPT_Reconstruction/run_ic_OPTtools_reconstruct.sh $ARGS_FULL
  
  # THIS ALSO WORKS
  $HOME/OPT_Reconstruction/run_ic_OPTtools_reconstruct.sh $ARGS_FULL
 	

fi
