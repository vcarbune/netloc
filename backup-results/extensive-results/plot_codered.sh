#!/bin/bash

echo "Plotting Logs in $1"

for networkType in p11558; do
        DEBUG=""
	gnuplot -e "var='Mass'; col=2; dir='$1'; network='${networkType}'; plottitle='Average Mass (${networkType} $DEBUG)';" objective_plot.gnu
	gnuplot -e "var='Difference'; col=4; dir='$1'; network='${networkType}'; plottitle='Distance in Probability to Source';" objective_plot.gnu

	# The ones below include EPFL data too.
	gnuplot -e "var='Rank'; col=6; dir='$1'; network='${networkType}'; plottitle='Average Rank';" default_yerror_plot.gnu
	gnuplot -e "var='Hops'; col=10; dir='$1'; network='${networkType}'; plottitle='Distance in Hops to Source';" default_yerror_plot.gnu
	gnuplot -e "var='nDCG'; col=8; dir='$1'; network='${networkType}'; plottitle='nDCG \@ 5%';" default_yerror_plot.gnu
	gnuplot -e "var='Localization'; col=12; dir='$1'; network='${networkType}'; plottitle='Localization Probability';" default_plot.gnu
done
