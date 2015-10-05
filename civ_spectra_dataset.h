/**
 civ_spectra_dataset.h
 Purpose: This file defines the class CIVSpectraDataset. This class contains the variables necessary to store a spectra dataset. This class is a specialization of the SpectraDataset class
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 10/02/2014
 
 */

#ifndef _LyaSpectraDataset_h
#define _LyaSpectraDataset_h

// libraries needed
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <stdio.h>
#include <vector>
////////

// classes needed
#include "dataset.h"
#include "input.h"
#include "z_dist_interpolation_map.h"
#include "lya_spectrum.h"
////////

// functions needed
////////

#include "typedefs.h"


class CIVSpectraDataset: public SpectraDataset{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs object and initializes its variables
    CIVSpectraDataset(const Input& input);
    
    // -------------------------------------------------------------
    // access methods
    
    // access function for flag_verbose_lya_spectra_dataset_
    size_t flag_verbose_civ_spectra_dataset() const {return flag_verbose_civ_spectra_dataset_;}
    
    // -------------------------------------------------------------
    // other methods
    
    // load the dataset
    void Load(const std::string& object_list, const std::string& lya_spectra_dir, const double& lya_wl);

    
    
    
    
private:
    // lya_spectra_dataset verbose flag
    size_t flag_verbose_civ_spectra_dataset_;
    
};

#endif
