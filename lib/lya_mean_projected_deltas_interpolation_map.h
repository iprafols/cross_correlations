/**
 lya_auto_interpolation_table.h
 Purpose: This file defines the class LyaMeanProjectedDeltasInterpolationTable. This class contains the variables necessary to interpolate the average projected delta as a function of redshift.
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 01/21/2015
 
 */

#ifndef _LyaMeanProjectedDeltasInterpolationMap_h
#define _LyaMeanProjectedDeltasInterpolationMap_h

// libraries needed
#include <iostream>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
////////

// classes needed
#include "input.h"
#include "interpolation_map.h"
////////

// functions needed
#include "function_split.hpp"
////////

#include "defines.h"

class LyaMeanProjectedDeltasInterpolationMap: public InterpolationMap{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs object and initializes its variables
    LyaMeanProjectedDeltasInterpolationMap(const Input& input);
    
private:

};


#endif
