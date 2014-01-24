set title sprintf("%s", plottitle)
set terminal png
set output sprintf("%s.png", filename)
set xlabel "Log2(Cluster Size)"
set ylabel "Observers Needed(%)"
set key left top
set xrange [2.8:9.5]
plot for [i=2:21] \
	sprintf("%s.dat", filename) using (log($1)):i notitle with linespoints
