set title sprintf("%s", plottitle)
set terminal png enhanced size 1024,768
set output sprintf("%s_%s.png", network, var)

set xlabel "Observers (%)"
set ylabel sprintf("%s", var)

set key left top
set xrange [-0.05:1.05]

plot \
	sprintf("%s_gbs.log", network) using 1:col title "GBS" with lines lc 1, \
	sprintf("%s_ec2.log", network) using 1:col title "EC2" with lines lc 3 lt 1 lw 3, \
	sprintf("%s_ec2_high.log", network) using 1:col title "HD" with lines lc 9, \
	sprintf("%s_voi.log", network) using 1:col title "VOI" with lines lc 4, \
	sprintf("%s_rand.log", network) using 1:col title "RAND" with lines lc 2, \
	sprintf("%s_epfl_high.log", network) using 1:col title "MLE + HD" with lines lc 6 lt 0 lw 2, \
	sprintf("%s_epfl_ec2.log", network) using 1:col title "MLE + EC2" with lines lc 7 lt 0 lw 2
