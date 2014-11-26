/**
 interpolation_table.cpp
 Purpose: This files contains the body for the functions defined in interpolation_table.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 09/30/2014
 */

#include "interpolation_map.h"

InterpolationMap::InterpolationMap(const Input& input){
    /**
     EXPLANATION:
     Cosntructs a InterpolationMap instance and initializes its variables
     
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
    double z_plus1 = 1.0-interpolation_step*0.5; // =1+z at mid_interval (auxiliar variable to speed up the computation)
    
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

double InterpolationMap::interpolation_map(std::map<double,double>::iterator it) const{
    /**
     EXPLANATION:
     Access function for interpolation_map_
     
     INPUTS:
     it - pointer to the element to be recovered
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     InterpolationMap
     
     FUNCITONS USED:
     NONE
     */
    if (it == interpolation_map_.end()){
        return _BAD_DATA_;
    }
    else{
        return (*it).second;
    }
}

double InterpolationMap::LinearInterpolation(const double& z) const{
    /**
     EXPLANATION:
     Compute the distances corresponding to the given redshift by using linear interpolation
     
     INPUTS:
     z - a double with the redshift value
     
     OUTPUTS:
     dist - the associated distance
     
     CLASSES USED:
     InterpolationMap
     
     FUNCITONS USED:
     NONE
     */
    
    std::map<double,double>::const_iterator it,it2;

    it = interpolation_map_.lower_bound(z);
    // no interpolation needed if distance is already computed for the given redshift
    if ((*it).first == z){
        return (*it).second;
    } else{
        it2 = interpolation_map_.lower_bound(z);
        it2 --;
        return ((*it).second-(*it2).second)/((*it).first-(*it2).first)*(z-(*it).first)+(*it).second;
    }
}