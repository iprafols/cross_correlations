/**
 plate.h
 Purpose: This file defines the class Plate. This class contains number of the plate and its mean values of right ascension and declination
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 06/26/2014
 
 */
#ifndef _Plate_h
#define _Plate_h

// libraries needed
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>

#include <CCfits>
////////

// classes needed
#include "sphere_point.h"
////////

// functions needed
////////

#include "defines.h"

class Plate{

public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs empty objects
    Plate(){};
    
    // constructs object and initializes its variables assigning zero values
    Plate(const int& plate_number);

    // constructs object and initializes its variables
    Plate(const std::string& filename);
    
    // -------------------------------------------------------------
    // access methods
    
    // access functions for angle_
    SpherePoint angle() const {return angle_;}
        
    // access function for number_of_objects_
    double number_of_objects() const {return number_of_objects_;}
    
    // access function for plate_num_
    int plate_number() const {return plate_number_;}
        
    // -------------------------------------------------------------
    // set methods
        
    // set function for plate_number
    void set_plate_number(int number) {plate_number_ = number;}
    
    // -------------------------------------------------------------
    // other methods
    
    // checks whethter or not the given plate is a neighbour plate
    bool IsNeighbour(const Plate& plate,const double& neighbours_max_distance);
    
    // normalizes the values of right ascension and declination and sets the number of objects to 0
    void Normalize();
    
    // updates the values of right ascension and declination
    void UpdateRADECValues(const Plate& plate);
    void UpdateRADECValues(const SpherePoint& angle);

    
    
    
    
protected:
    // number of the plate
    int plate_number_;
    
    
    
private:
    // angular position 
    SpherePoint angle_;
    
    // number of objects averaged (_NORM_ when already normalized)
    double number_of_objects_;
    
    
};    



#endif
