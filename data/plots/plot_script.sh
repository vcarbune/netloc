#!/bin/bash

rm -f *.png

# Mass concentration in solution cluster
gnuplot -e "filename='forest_pb'; plottitle='Forest Fire (MaxObservers=18%, N=100, Sources=20)'" script_pb.gnu
gnuplot -e "filename='barabasi_pb'; plottitle='Barabasi-Albert (MaxObservers=18%, N=100, Sources=20)'" script_pb.gnu
gnuplot -e "filename='erdos_pb'; plottitle='Erdos-Renyi (MaxObservers=32%, N=100, Sources=20)'" script_pb.gnu

# Observers needed to identify correct solution
gnuplot -e "filename='forest_test'; plottitle='Forest Fire (MinMass=25%, N=100)'" script_tst.gnu
gnuplot -e "filename='barabasi_test'; plottitle='Barabasi Albert (MinMass=25%, N=100, Sources=20)'" script_tst.gnu
gnuplot -e "filename='erdos_test'; plottitle='Erdos-Reniy (MinMass=25%, N=100, Sources=20)'" script_tst.gnu

# Average plots with yerrorbars
gnuplot -e "filename='forest_avg_pb'; plottitle='Forest Fire (MaxObservers=18%, N=100, Sources=20)'" script_avg_pb.gnu
gnuplot -e "filename='forest_avg_test'; plottitle='Forest Fire (MinMass=25%, N=100, Sources=20)'" script_avg_tst.gnu
gnuplot -e "filename='barabasi_avg_pb'; plottitle='Barabasi Albert (MaxObservers=18%, N=100, Sources=20)'" script_avg_pb.gnu
gnuplot -e "filename='barabasi_avg_test'; plottitle='Barabasi Albert (MinMass=25%, N=100, Sources=20)'" script_avg_tst.gnu
gnuplot -e "filename='erdos_avg_pb'; plottitle='Erdos-Reniy (MaxObservers=32%, N=100, Sources=20)'" script_avg_pb.gnu
gnuplot -e "filename='erdos_avg_test'; plottitle='Erdos-Reniy (MinMass=25%, N=100, Sources=20)'" script_avg_tst.gnu

# Only yerrorbars to identify decrease in dimensionality
gnuplot -e "filename='forest_avg_pb'; plottitle='Standard Error of Mass Concentration Plot (Forest Fire)'" script_yerror.gnu
gnuplot -e "filename='forest_avg_test'; plottitle='Standard Error of Min Observers Plot (Forest Fire)'" script_yerror.gnu
gnuplot -e "filename='barabasi_avg_pb'; plottitle='Standard Error of Mass Concentration Plot (Barabasi Albert)'" script_yerror.gnu
gnuplot -e "filename='barabasi_avg_test'; plottitle='Standard Error of Min Observers Plot (Barabasi Albert)'" script_yerror.gnu
gnuplot -e "filename='erdos_avg_pb'; plottitle='Standard Error of Mass Concentration Plot (Erdos-Reniy)'" script_yerror.gnu
gnuplot -e "filename='erdos_avg_test'; plottitle='Standard Error of Min Observers Plot (Erdos-Reniy)'" script_yerror.gnu
