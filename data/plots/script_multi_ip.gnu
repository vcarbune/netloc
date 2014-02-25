set terminal png
set output sprintf("%s.png", outfile)

set title sprintf("%s", plottitle)
set xlabel "ln(cluster size)"
set ylabel sprintf("%s", yl)
set key left top
set xrange [2.8:9.5]
plot \
	sprintf("%s.dat", file1) using (log($1)):22 title var."=".val1 with lines linecolor rgb "red", \
	sprintf("%s.dat", file2) using (log($1)):22 title var."=".val2 with lines linecolor rgb "blue", \
	sprintf("%s.dat", file3) using (log($1)):22 title var."=".val3 with lines linecolor rgb "green"
