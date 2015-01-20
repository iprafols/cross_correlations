/**
 pair.cpp
 Purpose: This files contains the body for the functions defined in pair.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 01/20/2015
 */

#include "pair.h"

Pair::Pair(double bad_data){
    /**
     EXPLANATION:
     Cosntructs a bad_data Pair instance
     
     INPUTS:
     bad_data - a double valued _BAD_DATA_
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Pair
     
     FUNCITONS USED:
     NONE
     */
    if (bad_data != _BAD_DATA_){
        std::cout << "Error while initializing a Pair 'bad data' instance" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    
    spectrum_ra_ = _BAD_DATA_;
    spectrum_dec_ = _BAD_DATA_;
    pixel_number_ = _BAD_DATA_INT_;
    pixel_dist_ = _BAD_DATA_;
    pixel_weight_ = _BAD_DATA_;
}

Pair::Pair(const double& spectrum_ra, const double& spectrum_dec, const int& pixel_number, const double& pixel_dist, const double& pixel_weight){
    /**
     EXPLANATION:
     Cosntructs a Pair instance
     
     INPUTS:
     spectrum_ra - spectrum's right ascension (in radians)
     spectrum_dec - spectrum's declination (in radians)
     pixel_number - pixel number
     pixel_dist - distance to pixel
     pixel_weight - pixel's weight
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Pair
     
     FUNCITONS USED:
     NONE
     */
    spectrum_ra_ = spectrum_ra;
    spectrum_dec_ = spectrum_dec;
    pixel_number_ = pixel_number;
    pixel_dist_ = pixel_dist;
    pixel_weight_ = pixel_weight;
    
}
