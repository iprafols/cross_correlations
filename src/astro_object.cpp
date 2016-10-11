/**
 astro_object.cpp
 Purpose: This files contains the body for the functions defined in astro_object.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 06/17/2014
 */

#include "astro_object.h"

AstroObject::AstroObject(double bad_data){
    /**
     EXPLANATION:
     Cosntructs a bad_data AstroObject instance
     
     INPUTS:
     bad_data - a double valued _BAD_DATA_
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     SpherePoint
     
     FUNCITONS USED:
     NONE
     */
    if (bad_data != _BAD_DATA_){
        std::cout << "Error while initializing a AstroObject 'bad data' instance" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    z_ = _BAD_DATA_;
    dist_ = _BAD_DATA_;
}

AstroObject::AstroObject(double& ra, double& dec, const int& plate, const int& fiber, const int& mjd, const double& z, const bool radians){
    /**
     EXPLANATION:
     Cosntructs a AstroObject instance
     
     INPUTS:
     ra - astronomical object's right ascension (in radians)
     dec - astronomical object's declination (in radians)
     plate - astronomical object's plate
     fiber - astronomical object's fiber
     mjd - astronomical object's Modified Jullian Date
     z - astronomical object's redshift
     radians - a boolean that specifies if the angles are given in radias - defaul = True
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     SpherePoint
     
     FUNCITONS USED:
     NONE
     */
    
    if (not radians){
        ra *= acos(-1)/180.0;
        dec *= acos(-1)/180.0;
    }
    SpherePoint angle(ra, dec);
    angle_ = angle;
    plate_ = plate;
    fiber_ = fiber;
    mjd_ = mjd;
    z_ = z;
    dist_ = _BAD_DATA_;
}


void AstroObject::SetDistance(const InterpolationMap& redshift_distance_map){
    /**
     EXPLANATION:
     Sets the distance to object
     
     INPUTS:
     redshif_distance_map - a InterpolationMap instance with the redshift-distance relation
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     
     FUNCITONS USED:
     NONE
     */
    
    dist_ = redshift_distance_map.LinearInterpolation(z_);
    
}
