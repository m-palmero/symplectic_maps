#define GNUPLOT "gnuplot --persist"

const char *data_path = "../data/";

void plot_gnuplot_general(double parameter, int iter, char file_name[])
{	
	FILE *gp;
	gp = popen(GNUPLOT, "w");
	// Setting terminal
	fprintf(gp, "set terminal png enhanced size 2500,2500 font \"../misc/cmr10.ttf\" 50 \n");
	fprintf(gp, "set colors classic \n");
	fprintf(gp, "set size square \n");
	// Setting title
	fprintf(gp, "set title \"Phase Space Test\\n{/*0.7 Control Parameter=%.4f; Max Iterations =%d;}\" \n", parameter, iter); 
	// Setting key (legend)
	//fprintf(gp, "set key font \"../misc/cmr10.ttf,35\" \n");
	//fprintf(gp, "set key outside \n");
	//fprintf(gp, "set key horiz \n");
	//fprintf(gp, "set key at screen 0.83,0.89 \n");
	// Setting frame
	// fprintf(gp, "unset xtics \n");
	// fprintf(gp, "unset ytics \n");
	// fprintf(gp, "unset cbtics \n");
	// fprintf(gp, "unset colorbox \n");
	// fprintf(gp, "unset border \n");
	// fprintf(gp, "set lmargin 0 \n");
	// fprintf(gp, "set rmargin 0 \n");
	// fprintf(gp, "set tmargin 0 \n");
	// fprintf(gp, "set bmargin 0 \n");	
	fprintf(gp, "set border 15 back lw 3 \n");
	fprintf(gp, "set autoscale yfix \n");
	fprintf(gp, "set autoscale xfix \n");
	//fprintf(gp, "deltax = 0.25 \n");
	//fprintf(gp, "deltay = 0.2 \n");
	//fprintf(gp, "set xrange[0:0+deltax] \n");
	//fprintf(gp, "set yrange[3.1415926-deltay:3.1415926+deltay] \n");;
	// Setting file's name
	fprintf(gp, "set output \'%s\' \n", file_name);
	// Plotting
	fprintf(gp, "plot '%sphase_space_%.4f.dat' w p pt 7 ps 0.5 lc 'grey10' noti,  '%sinitial_conditions_%.4f.dat' w p pt 7 ps 2 lw 3 lc 'forest-green' noti \n", data_path, parameter, data_path, parameter);
	// Plotting with iteration heat map
	//fprintf(gp, "plot 'phase_space.dat' u 1:2:3 w p pt 7 ps 0.5 palette noti \n");
	//fprintf(gp, "plot 'phase_space.dat' u 1:2:3 w p pt 7 ps 0.5 palette noti \n");
	fclose(gp);
}

void plot_gnuplot_boozer(double parameter, int iter, char filename[])
{	
	printf("Plotting...\n");
	FILE *gp;
	gp = popen(GNUPLOT, "w");
	// Setting terminal
	fprintf(gp, "set terminal png enhanced size 2500,2500 font \"../misc/cmr10.ttf\" 50 \n");
	fprintf(gp, "set colors classic \n");
	fprintf(gp, "set size square \n");
	// fprintf(gp, "unset xtics \n");
	// fprintf(gp, "unset ytics \n");
	// fprintf(gp, "set lmargin 0 \n");
	// fprintf(gp, "set rmargin 0 \n");
	// fprintf(gp, "set tmargin 0 \n");
	// fprintf(gp, "set bmargin 0 \n");	
	// Setting title
	fprintf(gp, "set title \"Single-null divertor tokamak map\\n{/*0.8 Control Parameter=%.7f}\" \n", parameter); 
	// Setting key (legend)
	// fprintf(gp, "set key font \"../misc/cmr10.ttf,35\" \n");
	// fprintf(gp, "set key outside \n");
	// fprintf(gp, "set key horiz \n");
	// fprintf(gp, "set key at screen 0.83,0.89 \n");
	// Setting frame
	fprintf(gp, "set border 15 back lw 3 \n");
	//fprintf(gp, "deltax = 0.0006 \n");
	//fprintf(gp, "deltay = 0.0002 \n");
	//fprintf(gp, "set xrange[0-deltax:0+deltax] \n");
	//fprintf(gp, "set yrange[0.997129-deltay:0.997129+deltay] \n");
	//fprintf(gp, "set autoscale yfix \n");
	//fprintf(gp, "set autoscale xfix \n");
	fprintf(gp, "set xrange [-0.002:0.002] \n");
	fprintf(gp, "set yrange[0.996:1.0001] \n");
	//fprintf(gp, "set xrange [-0.8:0.8] \n");
	//fprintf(gp, "set yrange[-0.6:1.4] \n");
	// Setting file's name
	fprintf(gp, "set output \'%s\' \n", filename);
	// Plotting
	//fprintf(gp, "plot 'phase_space_bw_%.3f.dat' w p pt 7 ps 0.7 lc 'grey10' noti, 'phase_space_fw_%.3f.dat' w p pt 7 ps 0.7 lc 'grey10' noti, 'upo_test.dat' w p pt 2 ps 5 lw 3 lc 'red' noti, 'unstable_manifold_%.3f.dat' w l lw 3 lc 'orange-red' noti, 'stable_manifold_%.3f.dat' w l lw 3 lc 'web-blue' noti, 'initial_conditions.dat' w p pt 7 ps 2 lc 'forest-green' noti\n", parameter, parameter, parameter,parameter);
	//fprintf(gp, "plot 'phase_space_fw_%.3f.dat' w p pt 7 ps 0.7 lc 'grey10' noti, 'upo_test.dat' w p pt 2 ps 5 lw 3 lc 'red' noti, 'stable_manifold_%.3f.dat' w l lw 2 lc 'orange-red' noti \n", parameter, parameter);
	//fprintf(gp, "plot '%sphase_space_inner_%.4f.dat' w p pt 7 ps 0.7 lc 'blue' noti, '%sphase_space_outer_%.4f.dat' w p pt 7 ps 0.7 lc 'red' noti \n",data_path, parameter, data_path, parameter);
	//fprintf(gp, "plot '%sphase_space_invariant_%.4f.dat' w p pt 7 ps 0.5 lc 'red' noti, '%sphase_space_%.4f.dat' w p pt 7 ps 1 lc 'grey10' noti, '%s/boozer/manifold_mag_saddle.dat' w l lw 2 lc 'red' noti \n",data_path, parameter, data_path, parameter, data_path);
	fprintf(gp, "plot '%sphase_space_%.7f.dat' w p pt 7 ps 0.5 lc 'grey10' noti, '%sinitial_conditions_%.7f.dat' w p pt 7 ps 1.5 lc 'web-green' noti \n",data_path, parameter, data_path, parameter);
	//fprintf(gp, "plot '%sphase_space_inner_%.4f.dat' w p pt 7 ps 0.7 lc 'blue' noti, \n",data_path, parameter);
	//fprintf(gp, "plot '%ssurvived_%.6f.dat' w p pt 7 ps 0.5 noti \n",data_path, parameter);
	//fprintf(gp, "plot '%sescape_%.4f.dat' w p pt 7 ps 0.5 lc 'grey10' noti \n",data_path, parameter);
	fclose(gp);
	printf("Plot completed!\n");
}

void plot_gnuplot_magnetic_footprint(double parameter, int iter, char filename[])
{	
	printf("Plotting...\n");
	FILE *gp;
	gp = popen(GNUPLOT, "w");
	// Setting terminal
	fprintf(gp, "set terminal png enhanced size 2500,2500 font \"../misc/cmr10.ttf\" 50 \n");
	fprintf(gp, "set colors classic \n");
	fprintf(gp, "set size square \n");
	//fprintf(gp, "unset xtics \n");
	//fprintf(gp, "unset ytics \n");
	// fprintf(gp, "set lmargin 0 \n");
	// fprintf(gp, "set rmargin 0 \n");
	// fprintf(gp, "set tmargin 0 \n");
	// fprintf(gp, "set bmargin 0 \n");	
	// Setting title
	fprintf(gp, "set title \"Single-null divertor tokamak map\\n{/*0.8 Control Parameter=%.3f}\" \n", parameter); 
	fprintf(gp, "set autoscale yfix \n");
	fprintf(gp, "set autoscale xfix \n");
	fprintf(gp, "set xrange [0:1] \n");
	fprintf(gp, "set xtics 0,0.5,1 \n");
	fprintf(gp, "set yrange[0:0.0014] \n");
	//fprintf(gp, "set log cb \n");
	fprintf(gp, "set cbrange[1:1000] \n");
	//fprintf(gp, "set palette cubehelix start 3 cycles 0.5 saturation 3 \n");
	fprintf(gp, "set palette defined (0 '#352a87',1 '#0363e1',2 '#1485d4',3 '#06a7c6',4 '#38b99e',5 '#92bf73',6 '#d9ba56',7 '#fcce2e',8 '#f9fb0e') \n");
	// Setting file's name
	fprintf(gp, "set output \'%s\' \n", filename);
	// Plotting
	fprintf(gp, "plot '%smag_footprint_%.4f.dat' u 1:2:3 w p pt 7 ps 1 palette noti \n",data_path, parameter);
	fclose(gp);
	printf("Plot completed!\n");
}

void plot_gnuplot_ullmann(double parameter, int iter, char filename[])
{	
	printf("Plotting...\n");
	FILE *gp;
	gp = popen(GNUPLOT, "w");
	// Setting terminal
	fprintf(gp, "set terminal png enhanced size 2500,2500 font \"../misc/cmr10.ttf\" 50 \n");
	fprintf(gp, "set colors classic \n");
	fprintf(gp, "set size square \n");
	// Setting title
	fprintf(gp, "set title \"Ergodic magnetic limiter map\\n{/*0.8 Control Parameter=%.7f}\" \n", parameter); 
	// Setting key (legend)
	// fprintf(gp, "set key font \"../misc/cmr10.ttf,35\" \n");
	// fprintf(gp, "set key outside \n");
	// fprintf(gp, "set key horiz \n");
	// fprintf(gp, "set key at screen 0.83,0.89 \n");
	// Setting frame
	fprintf(gp, "set border 15 back lw 3 \n");
	//fprintf(gp, "deltax = 0.0006 \n");
	//fprintf(gp, "deltay = 0.0002 \n");
	//fprintf(gp, "set xrange[0-deltax:0+deltax] \n");
	//fprintf(gp, "set yrange[0.997129-deltay:0.997129+deltay] \n");
	//fprintf(gp, "set xrange [0:6.28] \n");
	//fprintf(gp, "set yrange[0:0.4] \n");
	// Setting file's name
	fprintf(gp, "set output \'%s\' \n", filename);
	// Plotting
	fprintf(gp, "plot '%sphase_space_%.6f.dat' w p pt 7 ps 0.7 lc 'grey10' noti \n",data_path, parameter);
	fclose(gp);
	printf("Plot completed!\n");
}

void plot_gnuplot_fum(double parameter, int iter, char file_name[])
{	
	FILE *gp;
	gp = popen(GNUPLOT, "w");
	// Setting terminal
	fprintf(gp, "set terminal png enhanced crop size 2500,2500 font \"../misc/cmr10.ttf\" 50 \n");
	fprintf(gp, "set colors classic \n");
	fprintf(gp, "set size square \n");
	// Setting title
	fprintf(gp, "set title \"Simplified Fermi-Ulam model\\n{/*0.7 Control Parameter=%.3f; Max Iterations =%d;}\" \n", parameter, iter); 
	// Setting key (legend)
	//fprintf(gp, "set key font \"../misc/cmr10.ttf,35\" \n");
	//fprintf(gp, "set key outside \n");
	//fprintf(gp, "set key horiz \n");
	//fprintf(gp, "set key at screen 0.83,0.89 \n");
	// Setting frame
	// fprintf(gp, "unset xtics \n");
	// fprintf(gp, "unset ytics \n");
	// fprintf(gp, "unset cbtics \n");
	// fprintf(gp, "unset colorbox \n");
	// fprintf(gp, "unset border \n");
	// fprintf(gp, "set lmargin 0 \n");
	// fprintf(gp, "set rmargin 0 \n");
	// fprintf(gp, "set tmargin 0 \n");
	// fprintf(gp, "set bmargin 0 \n");	
	fprintf(gp, "set border 15 back lw 3 \n");
	
	fprintf(gp, "set xrange[0:6.28] \n");
	fprintf(gp, "set xtics (\"0\" 0, \"{/Symbol p}\" 3.14, \"2{/Symbol p}\" 6.28) \n");
	fprintf(gp, "set xlabel \"{/Symbol f}\" \n");

	fprintf(gp, "set yrange[0:0.07] \n");
	fprintf(gp, "set ytics (\"0\" 0, \"0.035\" 0.035, \"0.07\" 0.07) \n");
	fprintf(gp, "set ylabel \"V\" rotate by 0\n");

	// Setting file's name
	fprintf(gp, "set output \'%s\' \n", file_name);
	// Plotting
	fprintf(gp, "plot '%sphase_space_%.3f.dat' w p pt 7 ps 0.5 lc 'grey-10' noti,  '%sinitial_conditions_%.3f.dat' w p pt 7 ps 3 lw 4 lc 'forest-green' noti \n", data_path, parameter, data_path, parameter);
	// Plotting with iteration heat map
	//fprintf(gp, "plot 'phase_space.dat' u 1:2:3 w p pt 7 ps 0.5 palette noti \n");
	//fprintf(gp, "plot 'phase_space.dat' u 1:2:3 w p pt 7 ps 0.5 palette noti \n");
	fclose(gp);
}

void plot_gnuplot_histogram(char file_name[], double parameter)
{
	FILE *gp;
	gp = popen(GNUPLOT, "w");
	// Setting terminal
	fprintf(gp, "set terminal png enhanced size 2500,2500 font \"../misc/cmr10.ttf\" 50 \n");
	fprintf(gp, "set colors classic \n");
	fprintf(gp, "set size square \n");
	fprintf(gp, "n = 250 \n"); //number of intervals
	fprintf(gp, "min = 0.0 \n");
	fprintf(gp, "max = 1.0 \n");
	fprintf(gp, "width = abs((max - min) / n)\n");
	fprintf(gp, "hist(x,width)=width*floor(x/width)+width/2.0\n");
	//fprintf(gp, "set style data histogram \n");
	//fprintf(gp, "set style histogram clustered gap 0.0\n");
	//fprintf(gp, "set xrange[0:1] \n");
	//fprintf(gp, "set yrange[0:1000] \n");;
	fprintf(gp, "set boxwidth width*0.9 \n");
	fprintf(gp, "set style fill solid 0.5\n");
	fprintf(gp, "set output \'%s\' \n", file_name);
	//fprintf(gp, "set autoscale yfix \n");
	//fprintf(gp, "set autoscale xfix \n");
	fprintf(gp, "plot '%sescape_%.3f.dat' u (hist($1,width)) smooth freq w boxes lc 'black' noti \n", data_path, parameter);
	fclose(gp);
} 

void plot_gnuplot_time_series(char title[])
{	
	FILE *gp;
	gp = popen(GNUPLOT, "w");
	fprintf(gp, "set terminal png enhanced size 2500,2500 font \"../misc/cmr10.ttf\" 50 \n");
	fprintf(gp, "set colors classic \n");
	fprintf(gp, "set size square \n");
	//fprintf(gp, "unset key \n");
	//fprintf(gp, "unset xtics \n");
	//fprintf(gp, "unset ytics \n");
	fprintf(gp, "set autoscale yfix \n");
	fprintf(gp, "set autoscale xfix \n");
	//fprintf(gp, "set format y '' \n");
	//fprintf(gp, "set cbtics 0,500,1000 \n");
	//fprintf(gp, "set palette maxcolors 4\n");
	//fprintf(gp, "set xrange[0:6.28] \n");
	//fprintf(gp, "set yrange[0:1.0] \n");
	fprintf(gp, "set output '%s.png'\n",title);
	fprintf(gp, "plot 'orbit.dat' u 3 w l lw 1.5 lc 'black' noti\n");
	//fprintf(gp, "pause -1");
	fclose(gp);
}

void plot_gnuplot_RP_matrix(char title[])
{	
	FILE *gp;
	gp = popen(GNUPLOT, "w");
	
	fprintf(gp, "set terminal png enhanced size 2500,2500 font \"../misc/cmr10.ttf\" 50 \n");
	fprintf(gp, "set size square \n");
	fprintf(gp, "unset key \n");
	fprintf(gp, "unset xtics \n");
	fprintf(gp, "unset ytics \n");
	fprintf(gp, "unset colorbox \n");
	fprintf(gp, "set autoscale yfix \n");
	fprintf(gp, "set autoscale xfix \n");
	fprintf(gp, "set palette maxcolors 2\n");
	fprintf(gp, "set palette defined (0 'white', 1 'black')\n");
	fprintf(gp, "set output '%s.png'\n",title);
	fprintf(gp, "plot 'recurrence_plot.dat' matrix w image\n");
	//fprintf(gp, "pause -1");
	fclose(gp);
}

void plot_gnuplot_RP_coordinates(char filename[], double parameter)
{	
	FILE *gp;
	gp = popen(GNUPLOT, "w");
	fprintf(gp, "set terminal png enhanced size 2500,2500 font \"../misc/cmr10.ttf\" 50 \n");
	fprintf(gp, "set size square \n");
	//fprintf(gp, "unset key \n");
	//fprintf(gp, "unset xtics \n");
	//fprintf(gp, "unset ytics \n");
	fprintf(gp, "unset colorbox \n");
	fprintf(gp, "set autoscale yfix \n");
	fprintf(gp, "set autoscale xfix \n");
	//fprintf(gp, "set format x '' \n");
	//fprintf(gp, "set format y '' \n");
	//fprintf(gp, "set grid lw 4 \n");
	//fprintf(gp, "set xtics 0,500,1000 \n");
	//fprintf(gp, "set ytics 0,500,1000 \n");
	//fprintf(gp, "set log x\n");
	//fprintf(gp, "set log y\n");
	fprintf(gp, "set palette maxcolors 2\n");
	fprintf(gp, "set palette defined (0 'white', 1 'black')\n");
	fprintf(gp, "set output '%s'\n",filename);
	fprintf(gp, "plot '%srp.dat' w p pt 7 ps 1 lc 'black' noti\n", data_path);
	//fprintf(gp, "plot 'coord_rp.dat' w p pt 7 ps 1.0 lc 'black' notitle\n");
	//fprintf(gp, "plot 'coord_rp.dat' w p pt 7 ps 0.8 lc 'black' notitle\n");
	//fprintf(gp, "pause -1");
	fclose(gp);
}

//void plot_gnuplot_RPs_title(char filename[], int n, double eps, int noise_level, double noise_intensity)
void plot_gnuplot_RPs_title(char filename[], int n, double eps)
{	
	FILE *gp;
	gp = popen(GNUPLOT, "w");
	fprintf(gp, "set terminal png enhanced size 2500,2500 font \"../misc/cmr10.ttf\" 50 \n");
	fprintf(gp, "set size square \n");
	//fprintf(gp, "unset key \n");
	//fprintf(gp, "unset xtics \n");
	//fprintf(gp, "unset ytics \n");
	fprintf(gp, "unset colorbox \n");
	//fprintf(gp, "set autoscale yfix \n");
	//fprintf(gp, "set autoscale xfix \n");
	//fprintf(gp, "set xrange [1:%d] \n", n);
	//fprintf(gp, "set xrange [1:%d] \n", n);
	fprintf(gp, "set xrange [1:2000] \n");
	fprintf(gp, "set yrange [1:2000] \n");
	//fprintf(gp, "set format x '' \n");
	//fprintf(gp, "set format y '' \n");
	//fprintf(gp, "set grid lw 4 \n");
	fprintf(gp, "set border 15 back lw 3 \n");
	//fprintf(gp, "set xtics 0,500,1000 \n");
	//fprintf(gp, "set ytics 0,500,1000 \n");
	//fprintf(gp, "set xtics 0,1000,5000 \n");
	//fprintf(gp, "set ytics 0,1000,5000 \n");
	fprintf(gp, "set palette maxcolors 2\n");
	fprintf(gp, "set palette defined (0 'white', 1 'black')\n");
	// Setting title
	//fprintf(gp, "set title \"Test Recurrence Plot - Standard Map + Noise\\n{/*0.7 Distance Threshold = %.2f; Noise level = %d; Noise intensity = %.3f}\" \n", eps, noise_level, noise_intensity);
	fprintf(gp, "set title \"Test RPs normalized\"\n");
	fprintf(gp, "set output '%s'\n",filename);
	fprintf(gp, "plot 'rp.dat' w p pt 7 ps 0.25 lc 'black' notitle\n");
	//fprintf(gp, "plot 'coord_rp.dat' w p pt 7 ps 0.8 lc 'black' notitle\n");
	//fprintf(gp, "pause -1");
	fclose(gp);
}

void plot_gnuplot_line_y_parameters(char filename[], double paramater)
{
	printf("Plotting...\n");
	FILE *gp;
	gp = popen(GNUPLOT, "w");
	// Setting terminal
	fprintf(gp, "set terminal png enhanced background rgb 'black' size 2500,2500 font \"../misc/cmr10.ttf\" 50 \n");
	fprintf(gp, "set colors classic \n");
	fprintf(gp, "set size square \n");
	fprintf(gp, "set log cb \n");
	fprintf(gp, "set format cb '10^{%%L}' \n");
	fprintf(gp, "unset xtics \n");
	fprintf(gp, "unset ytics \n");
	fprintf(gp, "set border 15 back lw 3 \n");
	fprintf(gp, "set autoscale yfix \n");
	fprintf(gp, "set autoscale xfix \n");
	fprintf(gp, "set xrange [0.53:0.55] \n");
	fprintf(gp, "set yrange[0.999:1.0] \n");

	
	//fprintf(gp, "set xrange [-0.8:0.8] \n");
	//fprintf(gp, "set yrange[-0.6:1.4] \n");
	// Setting file's name
	fprintf(gp, "set output \'%s\' \n", filename);
	// Plotting
	fprintf(gp, "plot '%sline_x.dat' w p pt 7 ps 0.5 lc 'white' noti\n",data_path);
	fclose(gp);
	printf("Plot completed!\n");
}

void plot_gnuplot_transient_escape_measure(char filename[], double parameter)
{	
	printf("Plotting...\n");
	FILE *gp;
	gp = popen(GNUPLOT, "w");
	// Setting terminal
	fprintf(gp, "set terminal png enhanced size 2500,2500 font \"../misc/cmr10.ttf\" 50 \n");
	fprintf(gp, "set colors classic \n");
	fprintf(gp, "set size square \n");
	fprintf(gp, "set log cb \n");
	fprintf(gp, "set format cb '10^{%%L}' \n");
	fprintf(gp, "unset xtics \n");
	fprintf(gp, "unset ytics \n");
	//fprintf(gp, "set palette cubehelix start 3 cycles 0.5 saturation 3 \n");
	//fprintf(gp, "set palette cubehelix start 3 cycles 0.5 saturation 1 negative\n");
	fprintf(gp, "set palette defined (0 'white', 1 'medium-blue', 2 'medium-blue', 3 'web-green', 4'gold', 5 'orange-red', 6 'red', 7 'dark-red')\n");
	// Setting title
	//fprintf(gp, "set title \"Single-null divertor tokamak map\\n{/*0.8 Control Parameter=%.3f}\" \n", parameter); 
	// Setting key (legend)
	// fprintf(gp, "set key font \"../misc/cmr10.ttf,35\" \n");
	// fprintf(gp, "set key outside \n");
	// fprintf(gp, "set key horiz \n");
	// fprintf(gp, "set key at screen 0.83,0.89 \n");
	// Setting frame
	fprintf(gp, "set border 15 back lw 3 \n");
	//fprintf(gp, "deltax = 0.0006 \n");
	//fprintf(gp, "deltay = 0.0002 \n");
	//fprintf(gp, "set xrange[0-deltax:0+deltax] \n");
	//fprintf(gp, "set yrange[0.997129-deltay:0.997129+deltay] \n");
	fprintf(gp, "set autoscale yfix \n");
	fprintf(gp, "set autoscale xfix \n");
	// fprintf(gp, "set xrange [0:6.28] \n");
	// fprintf(gp, "set yrange[0.0:0.4] \n");
	// fprintf(gp, "set cbrange [7e-8:3e-5]\n");
	fprintf(gp, "set cbrange [1e-9:1e-4]\n");
	
	//fprintf(gp, "set xrange [-0.8:0.8] \n");
	//fprintf(gp, "set yrange[-0.6:1.4] \n");
	// Setting file's name
	fprintf(gp, "set output \'%s\' \n", filename);
	// Plotting
	fprintf(gp, "plot '%sescape_measure_%.7f.dat' u 1:2:3 w image noti, '%sinitial_conditions_%.7f.dat' lc rgb 'white' w d noti \n",data_path, parameter, data_path, parameter);
	fclose(gp);
	printf("Plot completed!\n");
}