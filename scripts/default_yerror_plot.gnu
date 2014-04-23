set title sprintf("%s", plottitle) font "arial,25"
set terminal svg enhanced size 1024,768
set output sprintf("%s/%s_%s.svg", dir, network, var)

set xlabel "Observers (%)" font "arial,25"
set ylabel sprintf("%s", var) font "arial,25"

set key right center font "arial, 20"
set xrange [-0.05:1.05]

set xtics font "arial, 25"
set ytics font "arial, 25"

plot \
	sprintf("%s/%s_gbs.log", dir, network) using 1:col title "GBS" with lines lc 1, \
	sprintf("%s/%s_gbs.log", dir, network) using 1:col:col+1 notitle with yerrorbars lc 1, \
	sprintf("%s/%s_ec2.log", dir, network) using 1:col title "EC2" with lines lc 3 lt 1 lw 3, \
	sprintf("%s/%s_ec2.log", dir, network) using 1:col:col+1 notitle with yerrorbars lc 3, \
	sprintf("%s/%s_ec2_high.log", dir, network) using 1:col title "HD" with lines lc 9, \
	sprintf("%s/%s_ec2_high.log", dir, network) using 1:col:col+1 notitle with yerrorbars lc 9, \
	sprintf("%s/%s_voi.log", dir, network) using 1:col title "VOI" with lines lc 4, \
	sprintf("%s/%s_voi.log", dir, network) using 1:col:col+1 notitle with yerrorbars lc 4, \
	sprintf("%s/%s_rand.log", dir, network) using 1:col title "RAND" with lines lc 2, \
	sprintf("%s/%s_rand.log", dir, network) using 1:col:col+1 notitle with yerrorbars lc 2, \
	sprintf("%s/%s_epfl_high.log", dir, network) using 1:col title "MLE + HD" with lines lc 6 lt 0 lw 2, \
	sprintf("%s/%s_epfl_high.log", dir, network) using 1:col:col+1 notitle with yerrorbars lc 6 lt 0 lw 2, \
	sprintf("%s/%s_epfl_ec2.log", dir, network) using 1:col title "MLE + EC2" with lines lc 7 lt 0 lw 2, \
	sprintf("%s/%s_epfl_ec2.log", dir, network) using 1:col:col+1 notitle with yerrorbars lc 7 lt 0 lw 2
