set title sprintf("%s", plottitle)
set terminal png
set output sprintf("%s_yerror.png", filename)
set xlabel "ln(cluster size)"
set ylabel "standard error (%)"
set key left top
set xrange [2.8:9.5]
plot sprintf("%s.dat", filename) using (log($1)):3 notitle with lines
