#!/bin/sh
#  smlm_ALYtools_Dispatch.slurm
#  
#  Created by Yuriy Alexandrov
#

### not sure how good is that - worked once for set of 3 small images
#SBATCH --nodes=1
#SBATCH --mem=128000
#SBATCH --time=2:00:00
#SBATCH --array=1-3

## Tubulin article
## #SBATCH --nodes=1
## #SBATCH --mem=0
## ##SBATCH --mem-per-cpu=200000
## #SBATCH --partition=hmem
## #SBATCH --time=2:00:00
## #SBATCH --array=1

ARGS="$1"
ARRARGS=(${ARGS//:/ })
mrt_camp=${ARRARGS[3]}

module load $mrt_camp

# add environment variables to args list
ARGS="$1:$SLURM_ARRAY_TASK_ID:$MRT"

echo $ARGS

sbatch $HOME/smlm_ALYtools_camp_software/smlm_ALYtools_NodeScript.sh $ARGS

exit
 
