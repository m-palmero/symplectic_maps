# reset
# set terminal pngcairo size 2000,2000 font 'Times,15'
# set colors classic
# set size square

# unset key
# unset ytics
# unset xtics

# set autoscale yfix 
# set autoscale xfix

# outfile = sprintf("phase_space_orbits.png")
# set output outfile
#plot 'results/standard_map/phase_space_full.dat' w p pt 7 ps 0.5 lc 'grey-10' noti, 'orbit_center_island.dat' w p pt 7 ps 1.5 lc 'blue' noti, 'orbit_p6_island.dat' w p pt 7 ps 1.5 lc 'magenta' noti, 'orbit_local_chaos.dat' w p pt 7 ps 1.5 lc 'forest-green' noti, 'orbit_chaos.dat' w p pt 7 ps 1.75 lc 'red' noti
#plot 'results/standard_map/phase_space_full.dat' w p pt 7 ps 0.5 lc 'grey-10' noti, 'orbit_intermit.dat' w p pt 7 ps 2 lc 'blue' noti, 'orbit_normal_chaos.dat' w p pt 7 ps 2 lc 'red' noti, 'initial_condition_intermit.dat' w p pt 2 ps 3.0 lw 2 lc 'blue' noti, 'initial_condition_normal_chaos.dat' w p pt 2 ps 4.0 lw 3 lc 'red' noti

## Zoom
reset
set terminal pngcairo size 2000,2000 font 'Times,15'
set colors classic
set size square

unset key
unset ytics
unset xtics

set xrange[3.8:4.8]
set yrange[5.2:5.7]

outfile = sprintf("phase_space_zoom_orbits.png")
set output outfile

plot 'phase_space.dat' w p pt 7 ps 0.5 lc 'grey-10' noti, 'orbit_intermit.dat' w p pt 7 ps 4 lc 'blue' noti, 'orbit_normal_chaos.dat' w p pt 7 ps 4 lc 'red' noti