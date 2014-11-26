/**
 interpolation_table.h
 Purpose: This file defines the class InterpolationTable. This class contains the variables necessary to store an interpolation table. 
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 09/30/2014
 
 */

#ifndef _InterpolationMap_h
#define _InterpolationMap_h

// libraries needed
#include <map>
////////

// classes needed
#include "input.h"
////////

// functions needed
////////

#include "defines.h"

class InterpolationMap{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs an empty object
    InterpolationMap(){};
    
    // constructs object and initializes its variables
    InterpolationMap(const Input& input);
    
    // -------------------------------------------------------------
    // access methods
    
    // access functions for interpolation_map_
    std::map<double,double> interpolation_map() const {return interpolation_map_;}
    double interpolation_map(double first) const {return (*interpolation_map_.find(first)).second;}
    double interpolation_map(std::map<double,double>::iterator it) const;
    
    // -------------------------------------------------------------
    // other methods
    
    // Compute the distances corresponding to the given redshift by using linear interpolation
    double LinearInterpolation(const double& z) const;
    
private:
    
    std::map<double,double> interpolation_map_;
    
};


#endif
