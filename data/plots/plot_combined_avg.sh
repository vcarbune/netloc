#!/bin/bash

for type in barabasi forest erdos; do
	gnuplot -e "file1='set4/${type}_avg_test'; ; val1='0.05'; \
		    file2='set6_eps_0.025/${type}_avg_test'; val2='0.025'; \
		    file3='set7_eps_0.100/${type}_avg_test'; val3='0.100'; \
		    outfile='${type}_eps_avgcombined_test'; \
		    yl='avg delta of observers (%)'; var='Eps'; \
		    plottitle='Observers Needed (${type}, N=100)'" \
	script_avgdiff_multi.gnu
done

for type in barabasi forest erdos; do
	gnuplot -e "file1='set4/${type}_avg_pb'; ; val1='0.05'; \
		    file2='set6_eps_0.025/${type}_avg_pb'; val2='0.025'; \
		    file3='set7_eps_0.100//${type}_avg_pb'; val3='0.100'; \
		    outfile='${type}_eps_avgcombined_pb'; \
		    yl='avg delta of mass (%)'; var='Eps'; \
		    plottitle='Avg Mass Concentration (Solution Cluster) (${type}, N=100)'" \
	script_avgdiff_multi.gnu
done

for type in barabasi forest erdos; do
	gnuplot -e "file1='set4/${type}_avg_test'; ; val1='4'; \
		    file2='set8_gensteps_2/${type}_avg_test'; val2='2'; \
		    file3='set9_gensteps_8/${type}_avg_test'; val3='8'; \
		    outfile='${type}_steps_avgcombined_test'; \
		    yl='avg delta of observers (%)'; var='Steps'; \
		    plottitle='Observers Needed (${type}, N=100)'" \
	script_avgdiff_multi.gnu
done

for type in barabasi forest erdos; do
	gnuplot -e "file1='set4/${type}_avg_pb'; ; val1='4'; \
		    file2='set8_gensteps_2/${type}_avg_pb'; val2='2'; \
		    file3='set9_gensteps_8/${type}_avg_pb'; val3='8'; \
		    outfile='${type}_steps_avgcombined_pb'; \
		    yl='avg delta of mass (%)'; var='Steps'; \
		    plottitle='Avg Mass Concentration (Solution Cluster) (${type}, N=100)'" \
	script_avgdiff_multi.gnu
done

