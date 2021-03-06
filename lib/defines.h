/**
 defines.h
 Purpose: This file contains function definitions and defined constants 
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 10/09/2014
 
 */


#ifndef _defines_h
#define _defines_h

#include <cfloat>
#include <climits>

// classes needed
#include "input.h"


// functions
void ComputePlateNeighbours(const Input& input);

// defined constants
#define _NORM_ -1
#define _BAD_DATA_ DBL_MAX
#define _BAD_DATA_INT_ INT_MAX



#endif
