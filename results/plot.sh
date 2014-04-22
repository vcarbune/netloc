#!/bin/bash

echo "Plotting Logs in $1"

NODES=$2

for networkType in barabasi erdos forest; do
        DEBUG=""
	NNAME=$(echo $networkType | sed 's/^\s*./\U&\E/g')

	# These only include plots with source estimation only.
	gnuplot -e "var='Mass'; col=2; dir='$1'; network='${networkType}${NODES}'; plottitle='Average Mass (${NNAME}$DEBUG)';" objective_plot.gnu
	gnuplot -e "var='Difference'; col=4; dir='$1'; network='${networkType}${NODES}'; plottitle='Distance in Probability to Source (${NNAME}$DEBUG)';" objective_plot.gnu

	# The ones below include EPFL data too.
	gnuplot -e "var='Rank'; col=6; dir='$1'; network='${networkType}${NODES}'; plottitle='Average Rank of Real Source (${NNAME}$DEBUG)';" default_yerror_plot.gnu
	gnuplot -e "var='Hops'; col=10; dir='$1'; network='${networkType}${NODES}'; plottitle='Hops to Real Source (${NNAME}$DEBUG)';" default_yerror_plot.gnu
	gnuplot -e "var='nDCG'; col=8; dir='$1'; network='${networkType}${NODES}'; plottitle='nDCG \@ Top 5\% (${NNAME}$DEBUG)';" default_yerror_plot.gnu
	gnuplot -e "var='Localization'; col=12; dir='$1'; network='${networkType}${NODES}'; plottitle='Localization Probability (${NNAME}$DEBUG)';" default_plot.gnu
done
