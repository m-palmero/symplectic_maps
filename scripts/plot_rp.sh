set terminal pngcairo size 2500,2500 font 'Times,30'
set size square
unset colorbox
#unset xtics
#unset ytic
set ytics 0,2500,10000
set xtics 0,2500,10000
set autoscale xfix
set autoscale yfix
set palette maxcolors 2
set palette defined (0 'white', 1 'black')
set output 'RP_high_rr.png'
plot 'coord_rp_high_rr.dat' w p pt 7 ps 0.25 lc 'blue' noti
