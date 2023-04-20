# include "dynamics.h"
// ------------------- Main Function ---------------------
# include "phase_space.c"
# include "escape.c"
# include "time_series.c"
# include "recurrence.c" 
# include "net_transitivity.c" 
# include "natural_measure.c" 
# include "escape_measure.c"
# include "winding_number.c"
# include "ensemble_statistics.c"

// ------------------ Global Parameters ------------------
char title[T], filename1[T], filename2[T], filename3[T], filename4[T], filename5[T], filename6[T];
int iter = 1e5; // Number of iterations
int nx = 1e2; // Number of initial conditions on x
int ny = 1e2; // Number of initial conditions on y

bool standard_map = false;
bool fermi_ulam_map = false;
bool boozer_map = false; 
bool ullmann_map = true;

MAP std_map;
MAP booz_map;
MAP ullm_map;
MAP sfu_map;
// -------------------------------------------------------

int main ()
{
    bool phase_space_analysis = false; 
    bool escape_analysis = false;
    bool time_series_analysis = false;
    bool recurrence_analysis = false;
    bool network_transitivity_analysis = false;
    bool natural_measure_analysis = false;
    bool escape_measure_analysis = true;
    bool winding_number = false;
    bool ensemble_statiscs = false;

    srand(time(NULL));

    if (standard_map)
    {
        std_map.parameter = 2.4;
        //double initial_conditions[] = {1.69236864663, 3.141592653589786};
        //double initial_conditions[] = {0.512345678987654, M_PI};
        // Center island:
        //double initial_conditions[] = {4.02196, 3.42837};
        // Period 6 island:
        //double initial_conditions[] = {4.76081, 3.70200};
        // Chaotic local region
        //double initial_conditions[] = {4.60996, 3.11897};
        // Chaotic global region
        //double initial_conditions[] = {0.5123, 3.1415};
        // Intermitency
        //double initial_conditions[] = {0.51234567899871, M_PI};
        //double delta = 1e-14;
        //double initial_conditions[] = {0.51234567899870, M_PI};
        double initial_conditions[] = {3.09680,  4.20571};
        std_map.x0_center = initial_conditions[1,0]; 
        std_map.y0_center = initial_conditions[0,1];
        std_map.y0_min = 0.0;
        std_map.y0_max = 2 * M_PI;
        std_map.x0_min = 0.0;
        std_map.x0_max = 2 * M_PI;
        std_map.x_space_min = 0.0;
        std_map.x_space_max = 2 * M_PI;
        std_map.y_space_min = 0.0;
        std_map.y_space_max = 2 * M_PI;
        //Zoom "intermitency" region
        //std_map.y0_min = 3.14;
        //std_map.y0_max = 4.5;
        //std_map.x0_min = 3.8;
        //std_map.x0_max = 4.8;
        //std_map.y0_min = 5.2;
        //std_map.y0_max = 5.7;
        //std_map.x_space_min = 3.8;
        //std_map.x_space_max = 4.8;
        //std_map.y_space_min = 5.2;
        //std_map.y_space_max = 5.7;
    }

    if (boozer_map)
    {
        //booz_map.parameter = 0.59300;
        booz_map.parameter = 0.6056;
        //Period 1 UPO (Magnetic saddle)
        //double initial_conditions[] = {3.8760479886142222e-16, 9.9999999999999967e-01};
        //Period 28 UPO (Ensemble 1 on Transient Analysis)
        //double initial_conditions[] = {-2.630556363014680e-14, 9.971293753820584e-01};
        //Period 57 UPO (Ensemble 2 on Transient Analysis)
        //double initial_conditions[] = {1.390380865995411e-14, 9.974541094372927e-01};
        //Period 29 UPO (Ensemble 3 on Transient Analysis)
        //double initial_conditions[] = {6.123834771568420e-04, 9.979545378272051e-01};
        //Period 30 UPO (Ensemble 4 on Transient Analysis)
        //double initial_conditions[] = {-9.515852739853789e-04, 9.978405500505392e-01};
        //double initial_conditions[] = {0.0,  0.996148};
        double initial_conditions[] = {0.0,  0.999315};
        booz_map.x0_center = initial_conditions[1,0]; 
        booz_map.y0_center = initial_conditions[0,1];
        booz_map.x0_min = -0.001;
        booz_map.x0_max = 0.001;
        booz_map.y0_min = 0.995;
        booz_map.y0_max = 0.999;
        booz_map.x_space_min = -0.002;
        booz_map.x_space_max = 0.002;
        booz_map.y_space_min = 0.995;
        booz_map.y_space_max = 1.0;
    }
    
    if (ullmann_map)
    {
        //ullm_map.parameter = 0.02006;
        ullm_map.parameter = 0.01836;
        //double initial_conditions[] = {1.724770290893, 0.3080149418021545};
        // Invariant spaning curve:
        //double initial_conditions[] = {1.99963, 0.363568};
        // Invariant island:
        //double initial_conditions[] = {4.24029, 0.337457};
        // Region around the cantori: 
        //double initial_conditions[] = {1.65195, 0.236157};
        // Around the period 7 UPO:
        double initial_conditions[] = {0.0, 0.30};
        ullm_map.x0_center = initial_conditions[1,0]; 
        ullm_map.y0_center = initial_conditions[0,1];
        ullm_map.x0_min = 0.0;
        ullm_map.x0_max = 2 * M_PI;
        ullm_map.y0_min = 0.0;
        ullm_map.y0_max = 0.4;
        ullm_map.x_space_min = 0.0;
        ullm_map.x_space_max = 2 * M_PI;
        ullm_map.y_space_min = 0.0;
        ullm_map.y_space_max = 0.4;
    }

    if (fermi_ulam_map)
    {
        sfu_map.parameter = 0.003;
        double initial_conditions[] = {5.92931,  0.0580766};
        //Low RR
        //double initial_conditions[] = {3.2202110883642980e+00, 1e-2};
        // Chaos
        //double initial_conditions[] = {M_PI, 1e-4};
        // Stickiness
        //double initial_conditions[] = {3.88005, 0.0351867};
        // Streched island 
        //double initial_conditions[] = {2.28586, 0.0406070};
        // Near FISC
        //double initial_conditions[] = {0.0, 0.065};
        sfu_map.x0_center = initial_conditions[1,0]; 
        sfu_map.y0_center = initial_conditions[0,1];
        sfu_map.x0_min = 0.0;
        sfu_map.x0_max = 2 * M_PI;
        sfu_map.y0_min = 0.0;
        sfu_map.y0_max = 0.07;
        sfu_map.x_space_min = 0.0;
        sfu_map.x_space_max = 2 * M_PI;
        sfu_map.y_space_min = 0.0;
        sfu_map.y_space_max = 0.1;
    }
    
    // ----------------- Phase Space analysis -------------------
    if (phase_space_analysis)
    {
        sprintf(filename1, "initial_conditions.dat");
        sprintf(filename2, "phase_space.dat");

        int lapis = 10; // Auxiliary for write less points
        double ic_box_size = 2e-12; // size of box of initial conitions
        bool localized_ensemble = false; //True if you want (nx*ny) i.c. in a localized spot

        if(standard_map) phase_space(&std_map, localized_ensemble, ic_box_size, lapis);
        if(boozer_map) phase_space(&booz_map, localized_ensemble, ic_box_size, lapis);
        if(ullmann_map) phase_space(&ullm_map, localized_ensemble, ic_box_size, lapis);
        if(fermi_ulam_map) phase_space(&sfu_map, localized_ensemble, ic_box_size, lapis);

        // // Ploting Phase Space:
        //sprintf(title, "phase_space_test");
        //plot_gnuplot_phase_space(title);
    }
    // ----------------------------------------------------------

    // ------------------- Escape analysis ----------------------
    if (escape_analysis)
    {
        sprintf(filename1, "escape.dat");
        
        //Localized ensemble for escape statistics:
        double ic_box_size = 1e-5;
        //Escape condition
        double escape_x = 0.0;
        double escape_y = 1.0;
        bool condition(double *v){return (*v > escape_y);};
        
        if(standard_map) escape(&std_map, condition, ic_box_size);
        if(boozer_map) escape(&booz_map, condition, ic_box_size);
        if(ullmann_map) escape(&ullm_map, condition, ic_box_size);
    }
    // ----------------------------------------------------------

    // ------------------ Time Series analysis ------------------
    if (time_series_analysis)
    {
        sprintf(filename1, "orbit.dat");
        sprintf(filename2, "initial_condition.dat");

        bool windowed = false;

        int ti = 1e0; // Auxiliary "time-variable" for writing at each "ti"
        int n_min, n_max;
        
        if(windowed)
        {
            n_min = 450; // Number of iteration for a windowed analysis (in a range n_min to n_max)
            n_max = 750; // Number of iteration for a windowed analysis (in a range n_min to n_max)
            if(standard_map) time_series(&std_map, ti, n_min, n_max);
            if(boozer_map) time_series(&booz_map, ti, n_min, n_max);
            if(ullmann_map) time_series(&ullm_map, ti, n_min, n_max);
            if(fermi_ulam_map) time_series(&sfu_map, ti, n_min, n_max);
        }
        else 
        {
            n_min = 0;
            n_max = 0;
            if(standard_map) time_series(&std_map, ti, n_min, n_max);
            if(boozer_map) time_series(&booz_map, ti, n_min, n_max);
            if(ullmann_map) time_series(&ullm_map, ti, n_min, n_max);
            if(fermi_ulam_map) time_series(&sfu_map, ti, n_min, n_max);
        }
        
        // // Ploting Time Series:
        //sprintf(title, "ts");
        //plot_gnuplot_time_series(title);

    }
    // ----------------------------------------------------------

    // ------------------ Recurrence analysis -------------------
    if (recurrence_analysis)
    {

        sprintf(filename1, "recurrence_plot.dat");
        sprintf(filename2, "eps_dependencies.dat");
        sprintf(filename3, "coord_rp.dat");

        double eps = 0.05;
        double eps_min, eps_max, multi;

        bool varying_eps = false;
        
        if (varying_eps)
        {
            eps_min = 1e-7;
            eps_max = 1e-1;
            multi = 1.8;//Multiplier between [1,2] to have "good" #pts/decade in log
            
            if(standard_map) recurrence_plot(&std_map, eps_min, eps_max, multi, eps);
            if(boozer_map) recurrence_plot(&booz_map, eps_min, eps_max, multi, eps);
            if(ullmann_map) recurrence_plot(&ullm_map, eps_min, eps_max, multi, eps);
            if(fermi_ulam_map) recurrence_plot(&sfu_map, eps_min, eps_max, multi, eps);
        }
        else
        {
            eps_min = 0;
            eps_max = 0;
            
            if(standard_map) recurrence_plot(&std_map, eps_min, eps_max, multi, eps);
            if(boozer_map) recurrence_plot(&booz_map, eps_min, eps_max, multi, eps);
            if(ullmann_map) recurrence_plot(&ullm_map, eps_min, eps_max, multi, eps);
            if(fermi_ulam_map) recurrence_plot(&sfu_map, eps_min, eps_max, multi, eps);
        }
        
    
        // // Ploting Recurrence Plot:
        //sprintf(title, "RP_matrix_test");
        //plot_gnuplot_RP_matrix(title);
        //sprintf(title, "RP");
        //plot_gnuplot_RP_coordinates(title, parameter);
    }
    // ----------------------------------------------------------

    // -------------- Net. Transitivity analysis ----------------
    if (network_transitivity_analysis)
    {

        sprintf(filename1, "adjacency_matrix.dat");
        sprintf(filename2, "eps_dependencies.dat");
        sprintf(filename3, "transitivity.dat");

        double eps = 0.01;
        int window = 50;

        if(standard_map) network_transitivity(&std_map, eps, window);
        if(boozer_map) network_transitivity(&booz_map, eps, window);
        if(ullmann_map) network_transitivity(&ullm_map, eps, window);
        if(fermi_ulam_map) network_transitivity(&sfu_map, eps, window);

    }
    // ----------------------------------------------------------

    // ----------------- Natural Measure analysis ---------------
    if (natural_measure_analysis)
    {
        sprintf(filename1, "natural_measure.dat");
        
        int size = 512;
        double ic_box_size = 1e-5; // size of box of initial conitions
        
        if(standard_map) natural_measure(&std_map, size, ic_box_size);
        if(boozer_map) natural_measure(&booz_map, size, ic_box_size);
        if(ullmann_map) natural_measure(&ullm_map, size, ic_box_size);
    }
    // ----------------------------------------------------------

    // ---------------- Escape Measure analysis -----------------
    if (escape_measure_analysis)
    {
        sprintf(filename1, "escape_measure.dat");	
        sprintf(filename2, "initial_conditions.dat");
        sprintf(filename3, "histogram.dat");
        sprintf(filename4, "histogram_measure.dat");
        sprintf(filename5, "trapped_orbits.dat");
        sprintf(filename6, "complexity.dat");
        
        int size = 1024;
        double ic_box_size = 1e-5; // size of box of initial conitions
        
        // Escape condition defined at line 115 of escape_measure.c
        if(standard_map) escape_measure(&std_map, size, ic_box_size);
        if(boozer_map) escape_measure(&booz_map, size, ic_box_size);
        if(ullmann_map) escape_measure(&ullm_map, size, ic_box_size);

        // // Ploting Transient escape measure:
        char plot_name[T];
	    sprintf(plot_name, "escape_measure_%.4f.png", ullm_map.parameter);
        plot_gnuplot_transient_escape_measure(plot_name);
    }
    // ----------------------------------------------------------

    // ---------------- Winding Number analysis -----------------
    if (winding_number)
    {
        sprintf(filename1, "winding_number_test.dat");
        sprintf(filename2, "reference.dat");
        if(standard_map) compute_winding_number(&std_map);
        if(boozer_map) compute_winding_number(&booz_map);
        if(ullmann_map) compute_winding_number(&ullm_map);
    }
    // ----------------------------------------------------------

    // ---------------- Ensemble statistics analysis ------------
    if (ensemble_statiscs)
    {
        sprintf(filename1, "ensemble_statistics_test.dat");
        sprintf(filename2, "test.dat");
        if(standard_map) evolve_ensemble_statistics(&std_map);
        if(boozer_map) evolve_ensemble_statistics(&booz_map);
        if(ullmann_map) evolve_ensemble_statistics(&ullm_map);
        if(fermi_ulam_map) evolve_ensemble_statistics(&sfu_map);
    }
    // ----------------------------------------------------------


    return 0;   
}

// int main()
// {
//     char plot_name[T];
//     sprintf(plot_name, "escape_measure_%.4f.png", ullm_map.parameter);
//     plot_gnuplot_transient_escape_measure(plot_name);

//     return 0; 
// }