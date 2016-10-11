/**
 dla_dataset.h
 Purpose: This file defines the class AstroObjectDataset. This class contains the variables necessary to store a object dataset. This class is a specialization of the cDataset class
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 11/25/2014
 
 */


#ifndef _QuasarDataset_h
#define _QuasarDataset_h

// libraries needed
#include <iostream>
#include <string>
#include <vector>

#include <CCfits>
////////

// classes needed
#include "astro_object_dataset.h"
#include "input.h"
////////

// functions needed
#include "function_split.hpp"
////////

#include "typedefs.h"


class QuasarDataset: public AstroObjectDataset{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs empty object
    QuasarDataset(){};
    
    // constructs object and initializes its variables
    QuasarDataset(const Input& input);
    
    // -------------------------------------------------------------
    // access methods
    
    // access function for flag_verbose_quasar_dataset_
    size_t flag_verbose_quasar_dataset() const {return flag_verbose_quasar_dataset_;}
    
    // -------------------------------------------------------------
    // other methods

    
    // load the dataset
    void Load(const double& z_min, const double& z_max, const std::string& dataset1);
    
private:
     // quasar_dataset verbose flag
    size_t flag_verbose_quasar_dataset_;
    
};





#endif
