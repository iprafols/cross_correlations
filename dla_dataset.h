/**
 dla_dataset.h
 Purpose: This file defines the class AstroObjectDataset. This class contains the variables necessary to store a object dataset. This class is a specialization of the cDataset class
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 11/25/2014
 
 */


#ifndef _DLADataset_h
#define _DLADdataset_h

// libraries needed
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
////////

// classes needed
#include "astro_object_dataset.h"
#include "input.h"
////////

// functions needed
#include "function_split.hpp"
////////

#include "typedefs.h"


class DLADataset: public AstroObjectDataset{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs empty object
    DLADataset(){};
    
    // constructs object and initializes its variables
    DLADataset(const Input& input);
    
    // load the dataset
    void Load(const double& z_min, const double& z_max, const std::string& dataset1);

};





#endif
