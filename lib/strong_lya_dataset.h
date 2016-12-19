/**
 strong_lya_dataset.h
 Purpose: This file defines the class StrongLyaDataset. This class contains the variables necessary to store a object dataset. This class is a specialization of the AstroObjectDataset class
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 12/16/2016
 
 */


#ifndef _StrongLyaDataset_h
#define _StrongLyaDataset_h

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


class StrongLyaDataset: public AstroObjectDataset{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs empty object
    StrongLyaDataset(){};
    
    // constructs object and initializes its variables
    StrongLyaDataset(const Input& input);
    
    // -------------------------------------------------------------
    // access methods
    
    // access function for flag_verbose_dla_dataset_
    size_t flag_verbose_strong_lya_dataset() const {return flag_verbose_strong_lya_dataset_;}

    // -------------------------------------------------------------
    // other methods
    
    // load the dataset
    void Load(const double& z_min, const double& z_max, const std::string& dataset1);
    
private:
    
    // dla_dataset verbose flag
    size_t flag_verbose_strong_lya_dataset_;
    
    // limits on the lyman alpha normalized flux
    double lya_flux_min_;
    double lya_flux_max_;

};





#endif
