set title "ForestFire Network (N=500)"
set terminal png
set output sprintf("%s.png", filename)
set xlabel "Cluster Size"
set ylabel "Percentage(%)"
set key left top
plot sprintf("%s.dat", filename) using 1:2 title "Tests(%)" with lines linetype 1, \
	sprintf("%s.dat", filename) using 1:2:3 notitle with yerrorbars linetype 1, \
	sprintf("%s.dat", filename) using 1:4 title "Identification(%)" with lines
