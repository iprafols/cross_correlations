/**
 lya_auto_interpolation_table.cpp
 Purpose: This files contains the body for the functions defined in lya_auto_interpolation_table.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 09/30/2014
 */

#include "lya_auto_interpolation_map.h"

LyaAutoInterpolationMap::LyaAutoInterpolationMap(const Input& input){
    /**
     EXPLANATION:
     Cosntructs a LyaAutoInterpolationMap instance and initializes its variables
     
     INPUTS:
     input - object of type Input
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Input
     InterpolationMap
     
     FUNCITONS USED:
     NONE
     */
    std::string line;
    int bin;
    double pi, value;
    
    // auxiliar variables
    int num_bins_lya_auto_correlation = input.num_bins_lya_auto_correlation();
    double step_lya_auto_correlation = input.step_lya_auto_correlation();
    
    std::ifstream file (input.lya_auto_correlation().c_str());
    if (file.is_open()){
        while (getline(file,line)){
            std::vector<std::string> cols = Split(line," ");
            bin = atoi(cols[0].c_str());
            
            // checking if the bin has r_perp = 0
            if (bin % num_bins_lya_auto_correlation == 0){
                pi = double(bin/num_bins_lya_auto_correlation)*step_lya_auto_correlation;
                value = atof(cols[1].c_str());
                
                // add to interpolation map
                interpolation_map_[pi] = value;
            }
        }
    }
    else{
        std::cout << "Error: in LyaAutoInterpolationMap::LyaAutoInterpolationMap : Could not read file" << std::endl;
        std::cout << input.lya_auto_correlation() << std::endl;
    }
}