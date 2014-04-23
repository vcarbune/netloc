#!/bin/bash

# Syntax: ./compare_"".sh $MOTHERDIR $SUBDIR1 $SUBDIR2

# Uses this: ./plot_combined.sh VAR FILE DIR1 VAL1 DIR2 VAL2 DIR3 VAL3 YLABEL COL TITLE SCRIPT

OUTDIR=sampling-$2-$1
mkdir combined/$OUTDIR

# Avg
for obj in ec2 gbs voi; do

titleobj=$(echo $obj | tr '[:lower:]' '[:upper:]')

./plot_combined.sh "" ${obj} \
                   $1/$2 "sampling" $1/$3 "fixed" $1/$4 "0.10" \
                   "Mass (%)" 2 \
                   "Average Mass (${titleobj})" script_multi.gnu \
                   "mass" $OUTDIR

./plot_combined.sh "" ${obj} \
                   $1/$2 "sampling" $1/$3 "fixed" $1/$4 "0.10" \
                   "Probability Distance (%)" 4 \
                   "Average Probability Distance (${titleobj})" script_multi.gnu \
                   "distance" $OUTDIR

./plot_combined.sh "" ${obj} \
                   $1/$2 "sampling" $1/$3 "fixed" $1/$4 "0.10" \
                   "Rank" 6 \
                   "Average Rank of Source Node (${titleobj})" script_multi.gnu \
                   "rank" $OUTDIR

./plot_combined.sh "" ${obj} \
                   $1/$2 "sampling" $1/$3 "fixed" $1/$4 "0.10" \
                   "nDCG" 8 \
                   "Average nDCG \@ Top 5\% (${titleobj})" script_multi.gnu \
                   "nDCG" $OUTDIR

./plot_combined.sh "" ${obj} \
                   $1/$2 "sampling" $1/$3 "fixed" $1/$4 "0.10" \
                   "Hops" 10 \
                   "Average Hops to Real Source (${titleobj})" script_multi.gnu \
                   "Hops" $OUTDIR
done
