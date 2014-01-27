set title sprintf("%s", plottitle)
set terminal png
set output sprintf("%s.png", filename)
set xlabel "ln(cluster size)"
set ylabel "avg delta of observers (%)"
set key left top
set xrange [2.8:9.5]
plot \
	sprintf("%s.dat", filename) using (log($1)):2 notitle with lines, \
	sprintf("%s.dat", filename) using (log($1)):2:3 notitle with yerrorbars
