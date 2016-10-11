/**
 function_to_str.cpp
 Purpose: This files contains the function ToStr()
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 06/17/2014
 */

#ifndef _toStr_h
#define _toStr_h


// libraries used
#include <sstream>
#include <string>
////////

// classes used
////////

// functions used
////////

//-----------------------------------------------------
template <typename T>
std::string ToStr(T i){
    /*
     EXPLANATION:
     transforms variable into string
     INPUTS:
     i - variable to be transformed
     OUTPUTS:
     string version of i
     */
    
    std::stringstream convert;
    convert << i;
    return convert.str();
    
}

#endif