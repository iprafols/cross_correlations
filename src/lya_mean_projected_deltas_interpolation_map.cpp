/**
 lya_auto_interpolation_table.cpp
 Purpose: This files contains the body for the functions defined in lya_auto_interpolation_table.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 09/30/2014
 */

#include "lya_mean_projected_deltas_interpolation_map.h"

LyaMeanProjectedDeltasInterpolationMap::LyaMeanProjectedDeltasInterpolationMap(const Input& input){
    /**
     EXPLANATION:
     Cosntructs a LyaMeanProjectedDeltasInterpolationMap instance and initializes its variables
     
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
    bool header = true;
    
    std::ifstream file (input.lya_projection_correction().c_str());
    if (file.is_open()){
        while (getline(file,line)){
            std::vector<std::string> cols = Split(line," ");
            
            if (header){
                header = false;
                continue;
            }
            interpolation_map_[atof(cols[0].c_str())] = atof(cols[1].c_str());
            
        }
        
    }
    else{
        std::cout << "Error: in LyaMeanProjectedDeltasInterpolationMap::LyaMeanProjectedDeltasInterpolationMap : Could not read file" << std::endl;
        std::cout << input.lya_projection_correction() << std::endl;
    }
}