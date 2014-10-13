/**
 lya_pixel.cpp
 Purpose: This files contains the body for the functions defined in lya_pixel.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 09/18/2014
 */

#include "lya_pixel.h"

LyaPixel::LyaPixel(const double& loglam, const double& lya_wl, const double& forest, const double& weight){
    /**
     EXPLANATION:
     Cosntructs a LyaPixel instance
     
     INPUTS:
     loglam - a double with the logarithm of the wavelength value (in Angstroms)
     lya_wl - a double with the wavelength of the lyman-alpha line (in Angstroms)
     forest - a double with the normalized flux in the ly-alpha forest
     weight - a double with the weight
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaPixel
     
     FUNCITONS USED:
     NONE
     */
    
    forest_ = forest;
    weight_ = weight;
    z_ = pow(10, loglam) / lya_wl - 1.0;

    
}

void LyaPixel::SetDistance(const InterpolationMap& redshift_distance_map){
    /**
     EXPLANATION:
     Sets the distance to object
     
     INPUTS:
     redshif_distance_map - a InterpolationMap instance with the redshift-distance relation
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaPixel
     
     FUNCITONS USED:
     NONE
     */
    
    dist_ = redshift_distance_map.LinearInterpolation(z_);

}