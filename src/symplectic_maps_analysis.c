# include "symplectic_maps_analysis.h"
// ------------------- Main Function ---------------------
# include "time_series.c"
# include "phase_space.c"
# include "escape_measure.c"
# include "recurrence.c"
// ------------------ Global Parameters ------------------
char data_file1[T], data_file2[T], data_file3[T], data_file4[T], data_file5[T], data_file6[T];
char plot_file[T];
char *ic_label;
int iter = 2e4; // Number of iterations
int nx = 1e2; // Number of initial conditions on x
int ny = 1; // Number of initial conditions on y

bool standard_map = false;
bool boozer_map = false; 
bool ullmann_map = true;
bool fermi_ulam_map = false;

bool phase_space_analysis = false; 
bool escape_measure_analysis = false;
bool recurrence_analysis = true;

MAP std_map;
MAP booz_map;
MAP ullm_map;
MAP sfu_map;
// -------------------------------------------------------

int maps_details()
{
    if (standard_map)
    {
        std_map.parameter = 2.4;
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
    }

    if (boozer_map)
    {
        double initial_conditions[] = {0.0, 0.0};
        
        int ic_set = 7;

        if(phase_space_analysis)
        {
            switch (ic_set)
            {
                case 1:
                    booz_map.parameter = 0.59300; //Local low escape rate
                    change(initial_conditions, 0.0006685, 0.99774007);
                    asprintf(&ic_label, "ensemble");
                    break;
                case 2:
                    booz_map.parameter = 0.6056; //Local high escape rate
                    change(initial_conditions, 0.0000001, 0.996982879);
                    asprintf(&ic_label, "ensemble");
                    break;
            }
        }

        if(recurrence_analysis)
        {
            switch (ic_set)
            {
                case 1:
                    change(initial_conditions, 0.0,  0.997169);
                    asprintf(&ic_label, "1st_ic");
                    break;
                case 2:
                    change(initial_conditions,-0.000287546,  0.997260);
                    asprintf(&ic_label, "2nd_ic");
                    break;
                case 3:
                    change(initial_conditions, 9.22709e-05,  0.997832);
                    asprintf(&ic_label, "3rd_ic");
                    break;
                case 4:
                    change(initial_conditions, -0.000426523,  0.998183);
                    asprintf(&ic_label, "4th_ic");
                    break;
                case 5:
                    booz_map.parameter = 0.59300; //Local low escape rate
                    change(initial_conditions, 0.0006685, 0.99774007);
                    asprintf(&ic_label, "ensemble");
                    break;
                case 6:
                    booz_map.parameter = 0.6056; //Local high escape rate
                    change(initial_conditions, 0.0000001, 0.997340171);
                    asprintf(&ic_label, "ensemble");
                    break;
                case 7: 
                    booz_map.parameter = 0.5930;
                    change(initial_conditions, 0.0006685, (0.99774007 - 2.1763257862517094e-11));
                    asprintf(&ic_label, "special");
                    break;
            }
        
        }
        
        booz_map.x0_center = initial_conditions[1,0]; 
        booz_map.y0_center = initial_conditions[0,1];
        
        booz_map.x0_min = -0.002;
        booz_map.x0_max = 0.002;
        booz_map.y0_min = 0.9965;
        booz_map.y0_max = 0.998;

        booz_map.x_space_min = -0.0025;
        booz_map.x_space_max = 0.0025;
        booz_map.y_space_min = 0.995;
        booz_map.y_space_max = 1.0;
    }

    if (ullmann_map)
    {
        double initial_conditions[] = {0.0, 0.0};
        
        int ic_set = 5;
        
        if(recurrence_analysis)
        {
            switch (ic_set)
            {
                case 1:
                    ullm_map.parameter = 0.015; //Standard
                    change(initial_conditions, 1.99963, 0.363568);
                    asprintf(&ic_label, "1st_ic");
                    break;
                case 2:
                    ullm_map.parameter = 0.015; //Standard
                    change(initial_conditions,  3.15057,  0.321440);
                    asprintf(&ic_label, "2nd_ic");
                    break;
                case 3:
                    ullm_map.parameter = 0.015; //Standard
                    change(initial_conditions, 3.58934,  0.332776);
                    asprintf(&ic_label, "3rd_ic");
                    break;
                case 4:
                    ullm_map.parameter = 0.015; //Standard  
                    change(initial_conditions, 1.65195, 0.236157);
                    asprintf(&ic_label, "4th_ic");
                    break;
                case 5:
                    ullm_map.parameter = 0.0133400; //Low escape rate
                    change(initial_conditions, 0.0145190817, 0.3101191903);  
                    asprintf(&ic_label, "ensemble");
                    break;
                case 6:
                    ullm_map.parameter = 0.0183600; //High escape rate
                    change(initial_conditions, 0.0117797591, 0.3111198483);
                    asprintf(&ic_label, "ensemble");
                    break;
            }
        }
        ullm_map.x0_center = initial_conditions[1,0]; 
        ullm_map.y0_center = initial_conditions[0,1];
        ullm_map.x0_min = 0.0;
        ullm_map.x0_max = 2 * M_PI;
        ullm_map.y0_min = 0.2;
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

}

int calculate_orbit_from_ic()
{

    if (boozer_map)
    {   
        sprintf(data_file1, "%sorbit_%.7f_%s.dat", data_path, booz_map.parameter, ic_label);
        time_series(&booz_map);
    }

    if (ullmann_map)
    {   
        sprintf(data_file1, "%sorbit_%.7f_%s.dat", data_path, ullm_map.parameter, ic_label);
        time_series(&ullm_map);
    }
}

int calculate_phase_spaces()
{
    if (phase_space_analysis)
    {
        int lapis = 1; // Auxiliary for write less points
        double ic_box_size = 1e-10; // size of box of initial conitions
        bool localized_ensemble = true; //True if you want (nx*ny) i.c. in a localized spot
        bool normalized_phase_space = false; //True if you want (nx*ny) i.c. in a localized spo

        if(standard_map) 
        {
            sprintf(data_file1, "%sinitial_conditions_%.4f.dat", data_path, std_map.parameter);
            sprintf(data_file2, "%sphase_space_%.4f.dat", data_path, std_map.parameter);
            phase_space(&std_map, localized_ensemble, normalized_phase_space, ic_box_size, lapis);
        }

        if(boozer_map) 
        {
            sprintf(data_file1, "%sinitial_conditions_%.7f.dat", data_path, booz_map.parameter);
            sprintf(data_file2, "%sphase_space_%.7f.dat", data_path, booz_map.parameter);
            phase_space(&booz_map, localized_ensemble, normalized_phase_space, ic_box_size, lapis);
        }
        
        if(ullmann_map) 
        {
            sprintf(data_file1, "%sinitial_conditions_%.7f.dat", data_path, ullm_map.parameter);
            sprintf(data_file2, "%sphase_space_%.7f.dat", data_path, ullm_map.parameter);    
            phase_space(&ullm_map, localized_ensemble, normalized_phase_space, ic_box_size, lapis);
        }

        if(fermi_ulam_map) 
        {
            sprintf(data_file1, "%sinitial_conditions_%.4f.dat", data_path, sfu_map.parameter);
            sprintf(data_file2, "%sphase_space_%.4f.dat", data_path, sfu_map.parameter);
            phase_space(&sfu_map, localized_ensemble, normalized_phase_space, ic_box_size, lapis);
        }
        
    }
}

int plot_phase_spaces()
{
    if(phase_space_analysis)
    {
        if (boozer_map)
        {
            sprintf(plot_file, "phase_space_%.7f.png", booz_map.parameter);
            plot_gnuplot_boozer(booz_map.parameter, iter, plot_file);
        }

        if(ullmann_map) 
        {
	        sprintf(plot_file, "phase_space_%.7f.png", ullm_map.parameter);
            plot_gnuplot_ullmann(ullm_map.parameter, iter, plot_file);
        }

    }
}

int calculate_transient_escape_measure()
{
    if (escape_measure_analysis)
    {

        int size = 1024;
        double ic_box_size = 1e-5; // size of box of initial conitions
        
        if(standard_map) 
        {
            sprintf(data_file1, "%sescape_measure_%.4f.dat", data_path, std_map.parameter);
            sprintf(data_file2, "%sinitial_conditions_%.4f.dat", data_path, std_map.parameter);
            sprintf(data_file3, "%shistogram_%.4f.dat", data_path, std_map.parameter);
            sprintf(data_file4, "%shistogram_measure_%.4f.dat", data_path, std_map.parameter);
            bool condition(double *v){return (*v > 1.0);};
            escape_measure(&std_map, size, ic_box_size, condition);
        }

        if(boozer_map) 
        {
            sprintf(data_file1, "%sescape_measure_%.7f_.dat", data_path, booz_map.parameter);
            sprintf(data_file2, "%sinitial_conditions_%.7f_.dat", data_path, booz_map.parameter);
            sprintf(data_file3, "%shistogram_%.7f_.dat", data_path, booz_map.parameter);
            sprintf(data_file4, "%shistogram_measure_%.7f_.dat", data_path, booz_map.parameter);
            bool condition(double *v){return (*v > 1.0);};
            escape_measure(&booz_map, size, ic_box_size, condition);
        }

        if(ullmann_map) 
        {
            sprintf(data_file1, "%sescape_measure_%.7f.dat", data_path, ullm_map.parameter);
            sprintf(data_file2, "%sinitial_conditions_%.7f.dat", data_path, ullm_map.parameter);
            sprintf(data_file3, "%shistogram_%.7f.dat", data_path, ullm_map.parameter);
            sprintf(data_file4, "%shistogram_measure_%.7f.dat", data_path, ullm_map.parameter);
            bool condition(double *v){return (*v < 0.0);};
            escape_measure(&ullm_map, size, ic_box_size, condition);
        }
        
        if(fermi_ulam_map)
        {
            sprintf(data_file1, "%sescape_measure_%.4f.dat", data_path, sfu_map.parameter);
            sprintf(data_file2, "%sinitial_conditions_%.4f.dat", data_path, sfu_map.parameter);
            sprintf(data_file3, "%shistogram_%.4f.dat", data_path, sfu_map.parameter);
            sprintf(data_file4, "%shistogram_measure_%.4f.dat", data_path, sfu_map.parameter);
            bool condition(double *v){return (*v > 1.0);};
            escape_measure(&sfu_map, size, ic_box_size, condition);
        }

    }
}

int plot_escape_measure()
{
    if(escape_measure_analysis)
    {
        if (boozer_map)
        {
            sprintf(plot_file, "escape_measure_%.7f.png", booz_map.parameter);
            plot_gnuplot_transient_escape_measure(plot_file, booz_map.parameter);
        }

        if(ullmann_map) 
        {
	        sprintf(plot_file, "escape_measure_%.7f.png", ullm_map.parameter);
            plot_gnuplot_transient_escape_measure(plot_file, ullm_map.parameter);
        }

    }
}

int calculate_recurrence_plots()
{
    if (recurrence_analysis)
    {
        //double eps = 1e-2;
        double eps = 0.01;
        int iter_effect_size = iter;

        if(boozer_map) 
        {
            sprintf(data_file1, "%srec_rate_%.7f_%s.dat", data_path, booz_map.parameter, ic_label);
            sprintf(data_file2, "%srec_plot_%.7f_%s.dat", data_path, booz_map.parameter, ic_label);
            recurrence_plot(&booz_map, iter_effect_size, eps);
        }

        if(ullmann_map) 
        {
            sprintf(data_file1, "%srec_rate_%.7f_%s.dat", data_path, ullm_map.parameter, ic_label);
            sprintf(data_file2, "%srec_plot_%.7f_%s.dat", data_path, ullm_map.parameter, ic_label);
            recurrence_plot(&ullm_map, iter_effect_size, eps);
        }
    
    }
}

int plot_recurrence_plots()
{
    if(recurrence_analysis)
    {
        if (boozer_map)
        {
            sprintf(plot_file, "rp_%.7f_%s.png", booz_map.parameter,ic_label);
            plot_gnuplot_RP_coordinates(plot_file, booz_map.parameter);
        }

        if(ullmann_map) 
        {
	        sprintf(plot_file, "rp_%.7f_%s.png", ullm_map.parameter, ic_label);
            plot_gnuplot_RP_coordinates(plot_file, ullm_map.parameter);
        }

    }
}

int calculate_recurrence_rate_ensemble()
{
    bool ensemble_in_x, ensemble_in_y;
    if (nx > 1) ensemble_in_x = true;
    else ensemble_in_x = false;
    if (ny > 1) ensemble_in_y = true;
    else ensemble_in_y = false;

    if (recurrence_analysis)
    {
        
        double ensemble_size = 1e-10;
        double multiplier_standard_deviation = 2.5;

        if(boozer_map)
        {
            double eps = 0.01;
            int iter_effect_size = 1e4;

            if(ensemble_in_x)
            {
                sprintf(data_file1, "%srr_%.7f_%s_x.dat", data_path, booz_map.parameter, ic_label);
                sprintf(data_file2, "%sstat_info_%.7f_%s_x.dat", data_path, booz_map.parameter, ic_label);
                sprintf(data_file3, "%saverage_rr_%.7f_%s_x.dat", data_path, booz_map.parameter, ic_label);
                sprintf(data_file4, "%slower_limit_%.7f_%s_x.dat", data_path, booz_map.parameter, ic_label);
                sprintf(data_file5, "%supper_limit_%.7f_%s_x.dat", data_path, booz_map.parameter, ic_label);
                recurrence_rate_ensemble(&booz_map, ensemble_size, iter_effect_size, multiplier_standard_deviation, eps);
            }

            if(ensemble_in_y)
            {
                sprintf(data_file1, "%srr_%.7f_%s_y.dat", data_path, booz_map.parameter, ic_label);
                sprintf(data_file2, "%sstat_info_%.7f_%s_y.dat", data_path, booz_map.parameter, ic_label);
                sprintf(data_file3, "%saverage_rr_%.7f_%s_y.dat", data_path, booz_map.parameter, ic_label);
                sprintf(data_file4, "%slower_limit_%.7f_%s_y.dat", data_path, booz_map.parameter, ic_label);
                sprintf(data_file5, "%supper_limit_%.7f_%s_y.dat", data_path, booz_map.parameter, ic_label);
                recurrence_rate_ensemble(&booz_map, ensemble_size, iter_effect_size, multiplier_standard_deviation, eps);
            }

        }

        if(ullmann_map)
        {
            double eps = 0.05;
            int iter_effect_size = iter;

            if(ensemble_in_x)
            {
                sprintf(data_file1, "%srr_%.7f_%s_x.dat", data_path, ullm_map.parameter, ic_label);
                sprintf(data_file2, "%sstat_info_%.7f_%s_x.dat", data_path, ullm_map.parameter, ic_label);
                sprintf(data_file3, "%saverage_rr_%.7f_%s_x.dat", data_path, ullm_map.parameter, ic_label);
                sprintf(data_file4, "%slower_limit_%.7f_%s_x.dat", data_path, ullm_map.parameter, ic_label);
                sprintf(data_file5, "%supper_limit_%.7f_%s_x.dat", data_path, ullm_map.parameter, ic_label);
                recurrence_rate_ensemble(&ullm_map, ensemble_size, iter_effect_size, multiplier_standard_deviation, eps);
            }

            if(ensemble_in_y)
            {
                sprintf(data_file1, "%srr_%.7f_%s_y.dat", data_path, ullm_map.parameter, ic_label);
                sprintf(data_file2, "%sstat_info_%.7f_%s_y.dat", data_path, ullm_map.parameter, ic_label);
                sprintf(data_file3, "%saverage_rr_%.7f_%s_y.dat", data_path, ullm_map.parameter, ic_label);
                sprintf(data_file4, "%slower_limit_%.7f_%s_y.dat", data_path, ullm_map.parameter, ic_label);
                sprintf(data_file5, "%supper_limit_%.7f_%s_y.dat", data_path, ullm_map.parameter, ic_label);
                recurrence_rate_ensemble(&ullm_map, ensemble_size, iter_effect_size, multiplier_standard_deviation, eps);
            }
        }
    }

}

int calculate_sticky_orbits()
{
    bool ensemble_in_x, ensemble_in_y;
    if (nx > 1) ensemble_in_x = true;
    else ensemble_in_x = false;
    if (ny > 1) ensemble_in_y = true;
    else ensemble_in_y = false;

    if (recurrence_analysis)
    {
        double eps = 0.01;
        double ensemble_size = 1e-10;
        double multiplier_standard_deviation = 3.0;
        int iter_effect_size = 1e4;

        if(boozer_map)
        {
            if(ensemble_in_x)
            {
                sprintf(data_file1, "%sstiky_orbits_%.7f_%s_x.dat", data_path, booz_map.parameter, ic_label);
                find_sticky_orbits(&booz_map, ensemble_size, iter_effect_size, multiplier_standard_deviation, eps);
            }

            if(ensemble_in_y)
            {
                sprintf(data_file1, "%sstiky_orbits_%.7f_%s_y.dat", data_path, booz_map.parameter, ic_label);
                find_sticky_orbits(&booz_map, ensemble_size, iter_effect_size, multiplier_standard_deviation, eps);
            }

        }
    }
}

int main ()
{
    double startTime, endTime, timeElapsed;
    srand(time(NULL));
    startTime = (float)clock() / CLOCKS_PER_SEC;   

    maps_details();    
    
    //calculate_phase_spaces();
    // // plot_phase_spaces();
    
    //calculate_orbit_from_ic();

    // //calculate_transient_escape_measure();
    // //plot_escape_measure();

    //calculate_recurrence_plots();
    //plot_recurrence_plots();

    calculate_recurrence_rate_ensemble();
    //calculate_sticky_orbits();
    
    endTime = (float)clock() / CLOCKS_PER_SEC;
    timeElapsed = endTime - startTime;
    printf("Time Elapsed  = %1.2f min\n", timeElapsed / 60);

    return 0;   
}
