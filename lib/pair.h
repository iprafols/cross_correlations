/**
 pair.h
 Purpose: This file defines the class Pair. This class contains the properties of a object-spectrum pair
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 06/17/2014
 
 */

#ifndef _Pair_h
#define _Pair_h

// libraries needed
#include <cstdlib>
#include <iostream>
////////

// classes needed
#include "sphere_point.h"
////////

// functions needed
////////

#include "defines.h"

class Pair{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs an empty object
    Pair(){};
    
    // constructs a "bad data" Pair
    Pair(double bad_data);
    
    // constructs object and initializes its variables
    Pair(const int& obj_plate, const int& object_num, const int& spec_plate, const int& spec_fiber, const int& spec_mjd, const int& pixel_number, const double& pixel_delta, const double& pixel_z, const double& pixel_weight);
    
    // -------------------------------------------------------------
    // access methods
    
    // object plate
    int obj_plate() const {return obj_plate_;}
    
    // object num in that plate's list
    int obj_num() const {return obj_num_;}
    
    // spectrum plate
    int spec_plate() const {return spec_plate_;}
    
    // spectrum fiber
    int spec_fiber() const {return spec_fiber_;}
    
    // spectrum mjd
    int spec_mjd() const {return spec_mjd_;}
    
    // pixel delta field
    double pixel_delta() const {return pixel_delta_;}
    
    // pixel number
    int pixel_number() const {return pixel_number_;}
    
    // pixel weight
    double pixel_weight() const {return pixel_weight_;}
    
    // pixel redshift
    double pixel_z() const {return pixel_z_;}
    
    // -------------------------------------------------------------
    // other methods
    
    
    
    
    
private:
    
    // object plate
    int obj_plate_;
    
    // object num in that plate's list
    int obj_num_;
    
    // spectrum plate
    int spec_plate_;
    
    // spectrum fiber
    int spec_fiber_;
    
    // spectrum mjd
    int spec_mjd_;
    
    // pixel delta field
    double pixel_delta_;
    
    // pixel number
    int pixel_number_;
    
    // pixel weight
    double pixel_weight_;
    
    // pixel redshift
    double pixel_z_;

    
};


#endif
