/**
 lya_auto_interpolation_table.h
 Purpose: This file defines the class LyaAutoInterpolationTable. This class contains the variables necessary to interpolate 1D Lyman-alpha autocorrelation from a paralel separation values.
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 01/21/2015
 
 */

#ifndef _LyaAutoInterpolationMap_h
#define _LyaAutoInterpolationMap_h

// libraries needed
#include <iostream>
#include <fstream>
#include <map>
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

class LyaAutoInterpolationMap: public InterpolationMap{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs object and initializes its variables
    LyaAutoInterpolationMap(const Input& input);
    
private:
    
    
};


#endif
