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
OMEtiffcount=${ARRARGS[2]}
mrt_camp=${ARRARGS[3]}
settings_file=${ARRARGS[4]}
src_filename_must_have_token=${ARRARGS[5]}
file_extension_template=${ARRARGS[6]}
jno=${ARRARGS[7]}
mrt=${ARRARGS[8]}

echo $SRCDIR
echo $DSTDIR
echo $OMEtiffcount
echo $jno
echo $mrt

fname="blank"

# CASE WHEN INPUTS ARE OME.tiff files
#
if (($OMEtiffcount != 0)); then
cnt=0
for f in $SRCDIR/*; do
	case "$f" in
	*$file_extension_template*)
	    #echo $(basename "$f")
	    cnt=$(($cnt + 1))
		if [ "x$jno" == "x$cnt" ]; then
			fname=$f
		fi
	;;
	esac
done
fi

echo $fname

if [ "$fname" != "blank" ]; then

  # this is the call of matlab-compiled sh file
  ARGS_FULL="$mrt $fname $DSTDIR $settings_file $src_filename_must_have_token" 
  
  echo $ARGS_FULL
  
  $HOME/smlm_ALYtools_camp_software/run_AI_Powered_2D_SMLM_reconstruct_CentOS7.sh $ARGS_FULL
 	
fi
