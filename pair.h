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
    Pair(const double& spectrum_ra, const double& spectrum_dec, const int& pixel_number, const double& pixel_dist, const double& pixel_weight);
    
    // -------------------------------------------------------------
    // access methods
    
    // distance to pixel
    double pixel_dist() const {return pixel_dist_;}
    
    // pixel number
    int pixel_number() const {return pixel_number_;}
    
    // pixel weight
    double pixel_weight() const {return pixel_weight_;}
    
    // spectrum angle
    SpherePoint spectrum_angle() const {return spectrum_angle_;}

    // -------------------------------------------------------------
    // other methods
    
    
    
    
    
private:
    // distance to pixel
    double pixel_dist_;
    
    // pixel number
    int pixel_number_;
    
    // pixel weight
    double pixel_weight_;
    
    // spectrum angle
    SpherePoint spectrum_angle_;
    

    
};


#endif
