/**
 z_dist_interpolation_table.h
 Purpose: This file defines the class ZDistInterpolationTable. This class contains the variables necessary to interpolate distances from a redshift values.
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 09/30/2014
 
 */

#ifndef _ZDistInterpolationMap_h
#define _ZDistInterpolationMap_h

// libraries needed
#include <map>
////////

// classes needed
#include "input.h"
#include "interpolation_map.h"
////////

// functions needed
////////

#include "defines.h"

class ZDistInterpolationMap: public InterpolationMap{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs object and initializes its variables
    ZDistInterpolationMap(const Input& input);
    
private:
    
    
};


#endif
