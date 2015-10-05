/**
 lya_pixel.h
 Purpose: This file defines the class LyaPixel. This class contains the variables necessary to store a lyman-alpha pixel. 
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 09/18/2014
 
 */

#ifndef _LyaPixel_h
#define _LyaPixel_h



// libraries needed
#include <cmath>
////////

// classes needed
#include "interpolation_map.h"
////////

// functions needed
////////

#include "defines.h"

class LyaPixel{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs a "bad data" LyaPixel
    LyaPixel(double bad_data);
    
    // constructs object and initializes its variables
    LyaPixel(const double& loglam, const double& lya_wl, const double& forest, const double& weight, const bool loglambda = true);
    
    // -------------------------------------------------------------
    // access methods
    
    // access method for dist_
    double dist() const {return dist_;}
    
    // access method for forest_
    double forest() const {return forest_;}
    
    // access method for weight_
    double weight() const {return weight_;}
    
    // access method for z_
    double z() const {return z_;}
    
    // -------------------------------------------------------------
    // other methods
    
    // set the distance to object
    void SetDistance(const InterpolationMap& redshift_distance_map);
    
    
    
private:
    // distance to pixel
    double dist_;
    
    // normalized flux in the ly-alpha forest
    double forest_;
    
    // weight
    double weight_;
    
    // redshift
    double z_;
    
};    


#endif
