/**
 astro_object.h
 Purpose: This file defines the class AstroObject. This class contains the basic properties of any astronomical object and also some useful computation related with this properties
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 06/17/2014
 
 */

#ifndef _AstroObject_h
#define _AstroObject_h

// libraries needed
////////

// classes needed
#include "interpolation_map.h"
#include "sphere_point.h"
////////

// functions needed
////////


class AstroObject{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs an empty object
    AstroObject(){};
    
    // constructs object and initializes its variables
    AstroObject(double& ra, double& dec, const int& plate, const int& fiber, const int& mjd, const double& z, const bool radians = true);

    // -------------------------------------------------------------
    // access methods
    
    // access functions for angle_
    SpherePoint angle() const {return angle_;}
    
    // access function for dist_
    double dist() const {return dist_;}
    
    // access function for fiber_
    int fiber() const {return fiber_;}
    
    // access function for mjd_
    int mjd() const {return mjd_;}
    
    // access function for plate_
    int plate() const {return plate_;}
        
    // access functions for z_
    double z() const {return z_;}
    
    // -------------------------------------------------------------
    // other methods
    
    // set the distance to object
    virtual void SetDistance(const InterpolationMap& redshift_distance_map);

    
    
    
    
protected:
    // angular position 
    SpherePoint angle_;
    
    // distance to object;
    double dist_;
    
    // fiber of the object
    int fiber_;
    
    // MJD of the object
    int mjd_;
    
    // plate of the object
    int plate_;
    
    // redshift of the object
    double z_;

    
};    


#endif
