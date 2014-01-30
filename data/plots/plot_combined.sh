#!/bin/bash

# Usage:
# ./plot_combined.sh VAR FILE DIR1 VAL1 DIR2 VAL2 DIR3 VAL3 YLABEL TITLE SCRIPT

for type in barabasi forest erdos; do
	gnuplot -e "file1='$3/${type}_$2'; ; val1='$4'; \
		    file2='$5/${type}_$2'; val2='$6'; \
		    file3='$7/${type}_$2'; val3='$8'; \
		    outfile='combined/${type}_${1}_${2}'; \
		    yl='${9}'; var='$1'; \
		    plottitle='${10} (${type}, N=100)'" \
	${11} # SCRIPT
done

