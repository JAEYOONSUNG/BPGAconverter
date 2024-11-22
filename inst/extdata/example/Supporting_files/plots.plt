set style data dots
set style function lines
set xzeroaxis lt 6 lc rgb "grey" lw 3
set terminal unknown
set xrange [0.95:]
f(x) = a*x**b
b = 1
a = 10000
fit f(x) 'pan_box.txt' using 6:7 via a, b
plot 'pan_box.txt' using 6:7 with dots lc rgb"black" title '', \
f(x) lc rgb "white"  lw 5 title ''


f1(x) = c*exp(d*x)
d = -1
c = 10000
fit f1(x) 'core_box.txt' using 6:7 via c, d
replot 'core_box.txt' using 6:7 with dots lc rgb"black" title '', \
f1(x) lc rgb"white"  lw 5 title ''

set terminal pdf enhanced size 8,6
set offset graph 0.0, graph 0.1, graph 0.1, graph 0.1
set output 'Core_Pan_Plot.pdf'
set title "Core and Pan Genome Plot" font "arial bold, 20" textcolor rgb"blue"
set xlabel "Number of Genomes\n\n" font "arial bold, 18" textcolor rgb"black"
set ylabel "\nNumber of Gene Families\n" font "arial bold, 18" textcolor rgb"black"
set key outside horizontal bottom title "" font "arial bold, 14" textcolor rgb"black"

set bars 3.0
set boxwidth  -2
set style fill solid border -1
set xtics 1 font "arial bold, 16" textcolor rgb"black"
set ytics 1000 font "arial bold, 16" textcolor rgb"black"


replot 'pan_box.txt' using 6:2:1:5:4 with candlesticks lc rgb"cyan" title 'Pan Genome' whiskerbars, \
''         using 6:3:3:3:3 with candlesticks lt 5 title ''
 
unset output

set output 'Core_Pan_Plot.pdf'

replot 'core_box.txt' using 6:2:1:5:4 with candlesticks lc rgb"pink"  title 'Core Genome' whiskerbars, \
''         using 6:3:3:3:3 with candlesticks lt 5 lc rgb"black" title '                                               Median Values'


unset output

exit gnuplot;

