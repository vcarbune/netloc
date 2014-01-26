set title sprintf("%s", plottitle)
set terminal png
set output sprintf("%s.png", filename)
set xlabel "ln(cluster size)"
set ylabel "mass in solution cluster (%)"
set key left top
set xrange [2.8:9.5]
plot for [i=2:21] \
	sprintf("%s.dat", filename) using (log($1)):i notitle with linespoints

set output sprintf("%s_ip.png", filename)
set xlabel "ln(cluster size)"
set ylabel "identification probability (%)"
plot sprintf("%s.dat", filename) using (log($1)):22 notitle with linespoints
