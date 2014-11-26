/**
 baofit_setup.h
 Purpose: This file defines the class BaofitSetup. This class contains the variables necessary setup the baofit ini file
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 11/25/2014
 
 */


#ifndef _BaofitSetup_h
#define _BaofitSetup_h

// libraries needed
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
////////

// classes needed
#include "input.h"
////////

// functions needed
////////


class BaofitSetup{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs empty object
    BaofitSetup(){};
    
    
    // -------------------------------------------------------------
    // other methods
    
    // Runs baofit
    void Run (const Input& input, const bool boostrap = false);

    // Sets ini file
    void Set(const Input& input, const bool boostrap = false);

private:

    // writes ini file required to run baofit
    void WriteIniFile(const Input& input, const bool boostrap = false);
    
};





#endif
