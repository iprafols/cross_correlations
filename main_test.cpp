/**
 main_correlation.cpp
 Purpose: Compute the cross correlation of the Lyman-alpha forest and quasars. Future versions should compute the cross-correlation of two species in general
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 06/09/2015
 */

// libraries used
#include <iostream>
#include <time.h>
#include <vector>
////////

// classes used
#include "input.h"
#include "lya_auto_interpolation_map.h"
#include "plots_object.h"
////////

// functions used
#include "defines.h"
////////

#include "typedefs.h"


int main(int argc, char *argv[]){
    /**
     EXPLANATION:
     Test the lya-autocorrelation
     
     INPUTS:
     input_file[optional] - a file containing the input settings
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Input
     LyaAutoInterpolationMap
     
     FUNCITONS USED:
     */
    
    // load time control variables
    time_t start_time,end_time;
    time(&start_time);
    
    // load global variables and plot object
    std::cout << "Initializing variables" << std::endl;
    string input_filename = "";
    if (argc > 1){
        input_filename += argv[1];
    }
    Input input(input_filename);
    size_t flag_verbose_main = input.flag_verbose_main();
    const PlotsObject kPlots(input);
    
    // compute the 1D lyman-alpha auto-correlation
    std::vector<LyaAutoInterpolationMap> lya_auto;
    for (size_t i = 0; i <= input.pixels_separation(); i++){
        lya_auto.push_back(LyaAutoInterpolationMap(input, i));
    }
    
    // plot the lyman-alpha auto-correlation
    kPlots.PlotLyaAuto(lya_auto, true);
    
    
    return 0;
}