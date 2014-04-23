#!/bin/bash

# Usage:
# ./plot_combined.sh VAR FILE DIR1 VAL1 DIR2 VAL2 DIR3 VAL3 YLABEL COL TITLE SCRIPT TYPE

N1=100
N2=239
N3=500

for type in barabasi forest erdos; do
	echo $3/${type}_$2
	echo $5/${type}_$2
	echo $7/${type}_$2
	NNAME=$(echo $type | sed 's/^\s*./\U&\E/g')
	gnuplot -e "file1='$3/${type}${N1}_$2'; val1='$4'; \
		    file2='$5/${type}${N2}_$2'; val2='$6'; \
		    file3='$7/${type}${N3}_$2'; val3='$8'; \
		    outfile='combined/${14}/${type}_${1}_${2}_${13}'; \
		    yl='${9}'; var='$1'; col=${10}; \
		    plottitle='${11} (${NNAME})'" \
	${12} # SCRIPT
done

