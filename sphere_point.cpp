/**
 sphere_point.cpp
 Purpose: This files contains the body for the functions defined in sphere_point.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 09/17/2014
 */

#include "sphere_point.h"

SpherePoint::SpherePoint(const double& ra, const double& dec){
    /**
     EXPLANATION:
     Cosntructs a SpherePoint instance
     
     INPUTS:
     ra - astronomical object's right ascension (in radians)
     dec - astronomical object's declination (in radians)
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     SpherePoint
     
     FUNCITONS USED:
     NONE
     */
    
    ra_ = ra;
    dec_ = dec;
    sin_dec_ = sin(dec_);
    cos_dec_ = cos(dec_);
}

double SpherePoint::CosAngularDistance(const SpherePoint& angle){
    /**
     EXPLANATION:
     Computes the angluar distance to antoher SpherePoint instance
     
     INPUTS:
     angle - SpherePoint instance to compute the angular distance with
     
     OUTPUTS:
     d - cosine of the angular distance
     
     CLASSES USED:
     SpherePoint
     
     FUNCITONS USED:
     NONE
     */
    
    return sin_dec_*angle.sin_dec()+cos_dec_*angle.cos_dec()*cos(ra_-angle.ra());
}

SpherePoint SpherePoint::operator+ (const SpherePoint& angle){
    /**
     EXPLANATION:
     Overloads the + operator. Returns a cSpherePoint instance with the ra_ and dec_ values incremented by the ra and dec values found in angle
     
     INPUTS:
     angle - a SpherePoint instance whose angles will be added
     
     OUTPUTS:
     temp - a SpherePoint instance with ra_ and dec_ divided by n
     
     CLASSES USED:
     SpherePoint
     
     FUNCITONS USED:
     NONE
     */
    
    double ra = ra_ + angle.ra();
    double dec = dec_ + angle.dec();
    
    SpherePoint temp(ra,dec);
    
    return temp;
}

void SpherePoint::operator+= (const SpherePoint& angle){
    /**
     EXPLANATION:
     Overloads the += operator. Updates the ra_ and dec_ values by addign those in angle. Recomputes sin_dec_ and cos_dec_
     
     INPUTS:
     angle - a cSpherePoint instance whose angles will be added
     
     OUTPUTS:
     temp - a cSpherePoint instance with ra_ and dec_ divided by n
     
     CLASSES USED:
     SpherePoint
     
     FUNCITONS USED:
     NONE
     */
    ra_ += angle.ra();
    dec_ += angle.dec();
    sin_dec_ = sin(dec_);
    cos_dec_ = cos(dec_);
    
}

SpherePoint SpherePoint::operator/ (const double& div){
    /**
     EXPLANATION:
     Overloads the / operator. Returns a cSpherePoint instance with the ra_ and dec_ values divided by div
     
     INPUTS:
     div - double number by which to divide
     
     OUTPUTS:
     temp - a SpherePoint instance with ra_ and dec_ divided by n
     
     CLASSES USED:
     SpherePoint
     
     FUNCITONS USED:
     NONE
     */
    
    double ra = ra_ / div;
    double dec = dec_ / div;

    SpherePoint temp(ra, dec);
    
    return temp;
}
void SpherePoint::operator/= (const double& div){
    /**
     EXPLANATION:
     Overloads the /= operator. Updates the ra_ and dec_ values by dividing them by div. Recomputes sin_dec_ and cos_dec_
     
     INPUTS:
     div - double number by which to divide
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     SpherePoint
     
     FUNCITONS USED:
     NONE
     */
    
    ra_ /= div;
    dec_ /= div;
    sin_dec_ = sin(dec_);
    cos_dec_ = cos(dec_);    
}

std::ostream& operator<< (std::ostream& out, const SpherePoint& angle){
    /**
     EXPLANATION:
     Overloads the << operator. Return "ra_value dec_value"
     
     INPUTS:
     out - an ostream 
     angle - a SpherePoint instance to insert into out
     
     OUTPUTS:
     out - a modified version of out with angle inserted into it
     
     CLASSES USED:
     SpherePoint
     
     FUNCITONS USED:
     NONE
     */
    out << angle.ra() << " " << angle.dec();
    
    return out;
}
