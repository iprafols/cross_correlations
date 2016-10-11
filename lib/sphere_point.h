/**
 sphere_point.h
 Purpose: This file defines the class SpherePoint. This class represents a position in the sky and knows how to compute angular distances and shuch
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 09/17/2014
 
 */

#ifndef _SpherePoint_h
#define _SpherePoint_h

// libraries needed
#include <cmath>
#include <fstream>
#include <iostream>
////////

// classes needed
////////

// functions needed
////////

class SpherePoint{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs empty object
    SpherePoint(){}
    
    // constructs object and initializes its variables
    SpherePoint(const double& ra, const double& dec);
    
    // -------------------------------------------------------------
    // access methods
    
    // access function for cos_dec_
    double cos_dec() const {return cos_dec_;}
    
    // access function for dec_
    double dec() const {return dec_;}
    
    // access function for ra_
    double ra() const {return ra_;}

    // access function for cos_dec_
    double sin_dec() const {return sin_dec_;}

    
    
    // -------------------------------------------------------------
    // other methods
    
    // compute the angular distance with another cSpherePoint instance
    double CosAngularDistance(const SpherePoint& angle);
    double AngularDistance(const SpherePoint& angle) {return acos(CosAngularDistance(angle));}
    
    
    
    // -------------------------------------------------------------
    // operator overload
    
    // adition
    SpherePoint operator+ (const SpherePoint& angle);
    void operator+= (const SpherePoint& angle);
    
    // division
    SpherePoint operator/ (const double& div);
    void operator/= (const double& div);

    
    
    
private:
    
    // cosine of the declination of the object
    double cos_dec_;
    
    // declination of the object
    double dec_;

    // right ascention of the object
    double ra_; 
    
    // sinus of the declination of the object
    double sin_dec_;



};

// stream extraction
std::ostream& operator<< (std::ostream& out, const SpherePoint& angle);




#endif
