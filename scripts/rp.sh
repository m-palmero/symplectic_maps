set terminal pngcairo size 2000,2000 font 'Times,30'
set size square
unset colorbox
set autoscale xfix
set autoscale yfix
set palette maxcolors 2
set palette defined (0 'white', 1 'black')
set output 'rp.png'
plot "rp.dat" w p pt 7 ps 1.0 lc 'black' notitle