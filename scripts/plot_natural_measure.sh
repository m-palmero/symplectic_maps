
set terminal png enhanced size 2500,2500 font "../misc/cmr10.ttf" 50  
set colors classic  
set size square  
set log cb  
set format cb "10^{%L}"  
unset xtics  
unset ytics 
unset cbtics 
unset colorbox
unset border
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
set palette defined (0 'white', 1 'medium-blue', 2 'medium-blue', 3 'web-green', 4'gold', 5 'orange-red', 6 'red', 7 'dark-red')
#set border 15 back lw 3  	 
#set xrange [0:6.28]  
set autoscale yfix
set autoscale xfix
#set yrange[0.0:0.4]  
set cbrange [1e-9:1e-4] 
#set cbrange [8e-9:1e-4]

# #set output "escape_measure_0.0133400.png"
# #plot 'escape_measure_0.0133400.dat' u 1:2:3 w image noti, 'initial_conditions.dat' w p pt 7 lc 'black' noti 

# set output "escape_measure_0.0183600.png"
# plot 'escape_measure_0.0183600.dat' u 1:2:3 w image noti, 'initial_conditions.dat' w p pt 7 lc 'black' noti   

# set output "escape_measure_0.5930000_.png"
# plot '../data/escape_measure_0.5930000.dat' u 1:2:3 w image noti

set output "escape_measure_0.6056000_.png"
plot '../data/escape_measure_0.6056000.dat' u 1:2:3 w image noti

#################### histogram #######################

# reset
# set terminal pngcairo size 800,800 font 'Helvetica,15'

# width = 1000.0;
# hist(x,width)=width*floor(x/width)+width/2.0
# set boxwidth width*0.9
# set style fill solid 0.5
# set tics out nomirror
# set xlabel "Escape time"
# set ylabel "Number of orbits"
# set format x "%1.0e"

# outfile = sprintf("figures/histogram.png")
# set output outfile
# plot 'results/histogram.dat' u (hist($1,width)):(1.0) smooth freq w boxes lc rgb "green" notitle
