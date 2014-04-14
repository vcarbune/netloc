set title sprintf("%s", plottitle)
set terminal png
set output sprintf("%s.png", outfile)
set xlabel "Observers (%)"
set ylabel sprintf("%s", yl)
set key left top
set xrange [2.8:9.5]
plot \
	sprintf("%s.dat", file1) using (log($1)):2 title var."=".val1 with lines linecolor rgb "red", \
	sprintf("%s.dat", file1) using (log($1)):2:3 notitle with yerrorbars linecolor rgb "red", \
	sprintf("%s.dat", file2) using (log($1)):2 title var."=".val2 with lines linecolor rgb "blue", \
	sprintf("%s.dat", file2) using (log($1)):2:3 notitle with yerrorbars linecolor rgb "blue", \
	sprintf("%s.dat", file3) using (log($1)):2 title var."=".val3 with lines linecolor rgb "green", \
	sprintf("%s.dat", file3) using (log($1)):2:3 notitle with yerrorbars linecolor rgb "green"
