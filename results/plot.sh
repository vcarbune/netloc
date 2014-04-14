#!/bin/bash

echo "Plotting Logs in $1"

NODES=$2

for networkType in barabasi erdos forest; do
        DEBUG="$3"
	gnuplot -e "var='Mass'; col=2; dir='$1'; network='${networkType}${NODES}'; plottitle='Average Mass (${networkType} $DEBUG)';" objective_plot.gnu
	gnuplot -e "var='Difference'; col=4; dir='$1'; network='${networkType}${NODES}'; plottitle='Distance in Probability to Source (${networkType} $DEBUG)';" objective_plot.gnu

	# The ones below include EPFL data too.
	gnuplot -e "var='Rank'; col=6; dir='$1'; network='${networkType}${NODES}'; plottitle='Average Rank (${networkType} $DEBUG)';" default_yerror_plot.gnu
	gnuplot -e "var='Hops'; col=10; dir='$1'; network='${networkType}${NODES}'; plottitle='Distance in Hops to Source (${networkType} $DEBUG)';" default_yerror_plot.gnu
	gnuplot -e "var='nDCG'; col=8; dir='$1'; network='${networkType}${NODES}'; plottitle='nDCG @ 5% (${networkType} $DEBUG)';" default_yerror_plot.gnu
	gnuplot -e "var='Localization'; col=12; dir='$1'; network='${networkType}${NODES}'; plottitle='Localization Probability (${networkType} $DEBUG)';" default_plot.gnu
done
