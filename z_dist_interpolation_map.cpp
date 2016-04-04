/**
 z_dist_interpolation_table.cpp
 Purpose: This files contains the body for the functions defined in z_dist_interpolation_table.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 09/30/2014
 */

#include "z_dist_interpolation_map.h"

ZDistInterpolationMap::ZDistInterpolationMap(const Input& input){
    /**
     EXPLANATION:
     Cosntructs a ZDistInterpolationMap instance and initializes its variables
     
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
    double interpolation_step,z_max_interpolation,z_min_interpolation;
    
    z_max_interpolation = input.z_max_interpolation();
    z_min_interpolation = input.z_min_interpolation();
    interpolation_step = (z_max_interpolation-z_min_interpolation)/double(input.num_points_interpolation());
    
    // setting initial and auxiliar variables
    double z = 0,dist = 0; // setting initial values for the integral
    double aux = input.c()/100.0; // =c/H0*h (auxiliar variable to speed up the computation)
    double wm = input.wm(); // Omega_matter (auxiliar variable to speed up the computation)
    double wv = 1.0-wm; // Omega_vacuum (auxiliar variable to speed up the computation)
    double z_plus1 = 1.0+interpolation_step*0.5; // =1+z at mid_interval (auxiliar variable to speed up the computation)
    
    // integrating from z = 0 to z_max_interpolation; midpoint rule
    while (z <= z_max_interpolation) {
        // integration step
        z += interpolation_step;
        z_plus1 += interpolation_step;
        dist += aux/sqrt(wv+wm*z_plus1*z_plus1*z_plus1)*interpolation_step; // in Mpc/h
        
        // if redshift is between z_min_interp and z_max_interp, saves values in the map
        if (z >= z_min_interpolation) { 
            interpolation_map_[z] = dist;
        }
    }
}