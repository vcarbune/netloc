set title sprintf("%s", plottitle) font "arial,25"
set terminal svg enhanced size 1024,768
set output sprintf("%s.svg", outfile)

set xrange [-0.05:1.05]
set xlabel "Observers (%)" font "arial,25"

set ylabel sprintf("%s", yl) font "arial,25"
set key right top font "arial, 20"

set xtics font "arial, 25"
set ytics font "arial, 25"

plot \
	sprintf("%s.log", file1) using 1:col title var."".val1 with lines linecolor rgb "red", \
	sprintf("%s.log", file1) using 1:col:col+1 notitle with yerrorbars linecolor rgb "red", \
	sprintf("%s.log", file2) using 1:col title var."".val2 with lines linecolor rgb "blue", \
	sprintf("%s.log", file2) using 1:col:col+1 notitle with yerrorbars linecolor rgb "blue", \
	sprintf("%s.log", file3) using 1:col title var."".val3 with lines linecolor rgb "green", \
	sprintf("%s.log", file3) using 1:col:col+1 notitle with yerrorbars linecolor rgb "green"
