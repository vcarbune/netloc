#!/bin/bash

# Syntax: ./compare_"".sh $MOTHERDIR $SUBDIR1 $SUBDIR2 $SUBDIR3

# Uses this: ./plot_combined_nodes.sh VAR FILE DIR1 VAL1 DIR2 VAL2 DIR3 VAL3 YLABEL COL TITLE SCRIPT

VAL1="100"
VAL2="240"
VAL3="500"

OUTDIR=size-$1+$2+$3+$4
mkdir combined/$OUTDIR

# Avg
for obj in ec2 gbs voi ec2_high rand; do

titleobj=$(echo $obj | tr '[:lower:]' '[:upper:]')

./plot_combined_nodes.sh "" ${obj} \
                   $1/$2 ${VAL1} $1/$3 ${VAL2} $1/$4 ${VAL3} \
                   "Mass (%)" 2 \
                   "Average Mass (${titleobj})" script_multi.gnu \
                   "mass" $OUTDIR

./plot_combined_nodes.sh "" ${obj} \
                   $1/$2 ${VAL1} $1/$3 ${VAL2} $1/$4 ${VAL3} \
                   "Probability Distance (%)" 4 \
                   "Average Probability Distance (${titleobj})" script_multi.gnu \
                   "distance" $OUTDIR

./plot_combined_nodes.sh "" ${obj} \
                   $1/$2 ${VAL1} $1/$3 ${VAL2} $1/$4 ${VAL3} \
                   "Rank" 6 \
                   "Average Rank of Source Node (${titleobj})" script_multi.gnu \
                   "rank" $OUTDIR

./plot_combined_nodes.sh "" ${obj} \
                   $1/$2 ${VAL1} $1/$3 ${VAL2} $1/$4 ${VAL3} \
                   "nDCG" 8 \
                   "Average nDCG \@ Top 5\% (${titleobj})" script_multi.gnu \
                   "nDCG" $OUTDIR

./plot_combined_nodes.sh "" ${obj} \
                   $1/$2 ${VAL1} $1/$3 ${VAL2} $1/$4 ${VAL3} \
                   "Hops" 10 \
                   "Average Hops to Real Source (${titleobj})" script_multi.gnu \
                   "Hops" $OUTDIR
done
