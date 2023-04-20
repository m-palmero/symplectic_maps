reset
set terminal pngcairo size 1500,1500 font 'Times,15'
set colors classic
set size square
unset border
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
unset key
unset ytics
unset xtics

set autoscale yfix 
set autoscale xfix

outfile = sprintf("ps.png")
set output outfile
#plot 'phase_space.dat' w p pt 7 ps 0.3 lc 'grey10' noti, 'initial_conditions.dat' w p pt 2 ps 4 lw 3 lc 'blue' noti
plot 'phase_space.dat' w p pt 7 ps 0.3 lc 'grey10' noti