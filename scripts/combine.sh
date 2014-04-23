#!/bin/bash

# Syntax: ./plot_combined.sh $MOTHERDIR $SUBDIR1 $SUBDIR2 $SUBDIR3

# Uses this: ./plot_combined.sh VAR FILE DIR1 VAL1 DIR2 VAL2 DIR3 VAL3 YLABEL COL TITLE SCRIPT

OUTDIR=eps-$2-$1

mkdir combined/$OUTDIR

# Avg
for obj in ec2 gbs voi ec2_high rand; do

./plot_combined.sh eps ${obj} \
                   $1/$2 "0.01" $1/$3 "0.05" $1/$4 "0.10" \
                   "Average Mass (%)" 2 \
                   "Noise Impact on Average Mass - ${obj}" script_multi.gnu \
                   "mass" $OUTDIR

./plot_combined.sh eps ${obj} \
                   $1/$2 "0.01" $1/$3 "0.05" $1/$4 "0.10" \
                   "Average Distance in Probability (%)" 4 \
                   "Noise Impact on Average Pb Distance - ${obj}" script_multi.gnu \
                   "distance" $OUTDIR

./plot_combined.sh eps ${obj} \
                   $1/$2 "0.01" $1/$3 "0.05" $1/$4 "0.10" \
                   "Average Rank of Source Node" 6 \
                   "Noise Impact on Average Rank - ${obj}" script_multi.gnu \
                   "rank" $OUTDIR

./plot_combined.sh eps ${obj} \
                   $1/$2 "0.01" $1/$3 "0.05" $1/$4 "0.10" \
                   "Average nDCG @ Top 5%" 8 \
                   "Noise Impact on nDCG - ${obj}" script_multi.gnu \
                   "nDCG" $OUTDIR

./plot_combined.sh eps ${obj} \
                   $1/$2 "0.01" $1/$3 "0.05" $1/$4 "0.10" \
                   "Average Hops to Source Node" 10 \
                   "Noise Impact on Hops - ${obj}" script_multi.gnu \
                   "Hops" $OUTDIR
done
