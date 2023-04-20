reset
set terminal pngcairo size 1000,1000 font 'Helvetica,15'

n=250 #number of intervals
min=0. #min value
max=6.28 #max value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set boxwidth width*0.9
set style fill solid 0.5 # fill style

unset key
#unset xtics

#set xrange [0:6.28]

outfile = sprintf("test.png")
set output outfile
plot "results/escape_test.dat" u (hist($2,width)):(1.0) smooth freq w boxes lc rgb "black" notitle
