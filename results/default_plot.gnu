set title sprintf("%s", plottitle)
set terminal png enhanced size 1024,768
set output sprintf("%s/%s_%s.png", dir, network, var)

set xlabel "Observers (%)"
set ylabel sprintf("%s", var)

set key left top
set xrange [-0.05:1.05]

plot \
	sprintf("%s/%s_ec2.log", dir, network) using 1:col title "EC2" with lines, \
	sprintf("%s/%s_gbs.log", dir, network) using 1:col title "GBS" with lines, \
	sprintf("%s/%s_voi.log", dir, network) using 1:col title "VOI" with lines, \
	sprintf("%s/%s_rand.log", dir, network) using 1:col title "RAND" with lines, \
	sprintf("%s/%s_epfl_high.log", dir, network) using 1:col title "EPFL " with lines, \
	sprintf("%s/%s_epfl_ec2.log", dir, network) using 1:col title "EPFL + EC2" with lines
	
