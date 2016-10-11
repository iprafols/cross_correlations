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

    obj_plate_ = _BAD_DATA_INT_;
    obj_num_ = _BAD_DATA_INT_;
    spec_plate_ = _BAD_DATA_INT_;
    spec_fiber_ = _BAD_DATA_INT_;
    spec_mjd_ = _BAD_DATA_INT_;
    pixel_number_ = _BAD_DATA_INT_;
    pixel_delta_ = _BAD_DATA_;
    pixel_weight_ = _BAD_DATA_;
    pixel_z_ = _BAD_DATA_;
}

Pair::Pair(const int& obj_plate, const int& obj_num, const int& spec_plate, const int& spec_fiber, const int& spec_mjd, const int& pixel_number, const double& pixel_delta, const double& pixel_z, const double& pixel_weight){
    /**
     EXPLANATION:
     Cosntructs a Pair instance
     
     INPUTS:
     obj_plate - plate the object is found in
     obj_num - number of object in that plate list
     spec_plate - spectrum's plate
     spec_fiber - spectrum's fiber
     spec_mjd - spectrum's MJD
     pixel_number - pixel number
     pixel_delta - pixel's delta field
     pixel_z - pixel's redshift
     pixel_weight - pixel's weight
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Pair
     
     FUNCITONS USED:
     NONE
     */
    
    obj_plate_ = obj_plate;
    obj_num_ = obj_num;
    spec_plate_ = spec_plate;
    spec_fiber_ = spec_fiber;
    spec_mjd_ = spec_mjd;
    pixel_number_ = pixel_number;
    pixel_delta_ = pixel_delta;
    pixel_weight_ = pixel_weight;
    pixel_z_ = pixel_z;
    
}
