/**
 lya_pixel.cpp
 Purpose: This files contains the body for the functions defined in lya_pixel.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 09/18/2014
 */

#include "lya_pixel.h"

LyaPixel::LyaPixel(double bad_data){
    /**
     EXPLANATION:
     Cosntructs a bad_data LyaPixel instance
     
     INPUTS:
     bad_data - a double valued _BAD_DATA_
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaPixel
     
     FUNCITONS USED:
     NONE
     */
    if (bad_data != _BAD_DATA_){
        std::cout << "Error while initializing a LyaPixel 'bad data' instance" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    
    delta_ = _BAD_DATA_;
    weight_ = _BAD_DATA_;
    z_ = _BAD_DATA_;
    
}

LyaPixel::LyaPixel(const double& loglam, const double& lya_wl, const double& delta, const double& weight, const bool loglambda){
    /**
     EXPLANATION:
     Cosntructs a LyaPixel instance
     
     INPUTS:
     loglam - a double with the logarithm of the wavelength value (in Angstroms)
     lya_wl - a double with the wavelength of the lyman-alpha line (in Angstroms)
     delta - a double with the ly-alpha delta field
     weight - a double with the weight
     loglambda - a boolean specifying if loglam is given as the logarithm of the of the wavelength (true) or the wavelength itself (false)
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaPixel
     
     FUNCITONS USED:
     NONE
     */
    
    delta_ = delta;
    weight_ = weight;
    loglam_ = loglam
    if (loglambda){
        z_ = pow(10, loglam) / lya_wl - 1.0;
    }
    else{
        z_ = loglam / lya_wl - 1.0;
    }
    
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