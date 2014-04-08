set title sprintf("%s", plottitle)
set terminal png enhanced size 1024,768
set output sprintf("%s_%s.png", network, var)

set xlabel "Observers (%)"
set ylabel sprintf("%s", var)

set key left top
set xrange [-0.05:1.05]

plot \
	sprintf("%s_ec2.log", network) using 1:col title "EC2" with lines, \
	sprintf("%s_ec2_high.log", network) using 1:col title "HD" with lines, \
	sprintf("%s_gbs.log", network) using 1:col title "GBS" with lines, \
	sprintf("%s_voi.log", network) using 1:col title "VOI" with lines, \
	sprintf("%s_rand.log", network) using 1:col title "RAND" with lines, \
	sprintf("%s_epfl_high.log", network) using 1:col title "MLE + HD" with lines, \
	sprintf("%s_epfl_ec2.log", network) using 1:col title "MLE + EC2" with lines
