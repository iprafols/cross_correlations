/**
 lya_spectra_dataset.h
 Purpose: This file defines the class LyaSpectraDataset. This class contains the variables necessary to store a spectra dataset. This class is a specialization of the cDataset class
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 09/17/2014
 
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
#include "interpolation_map.h"
#include "lya_spectrum.h"
////////

// functions needed
////////

#include "typedefs.h"


class LyaSpectraDataset: public Dataset{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs object and initializes its variables
    LyaSpectraDataset(const Input& input);
    
    // -------------------------------------------------------------
    // access methods
    
    // access functions for list_
    PlatesMapVector<LyaSpectrum>::map list() const {return list_;}
    std::vector<LyaSpectrum> list(int plate_number) const {return (*list_.find(plate_number)).second;}
    LyaSpectrum list(int plate_number, int pos) const {return (*list_.find(plate_number)).second[pos];}
    
    // -------------------------------------------------------------
    // other methods
            
    // find the catalog length
    int FindCatalogLength(const std::string& lya_spectra_catalog);
        
    // adds the objects' RA-DEC values to ostream
    void GiveRADEC(std::ostream& out) const;
    
    // adds the objects' redshift values to ostream
    void GiveZ(std::ostream& out) const;
    
    // load the dataset
    void Load(const std::string& object_list, const std::string& lya_spectra_dir, const double& lya_wl);

    // set the distance to every object in the dataset
    void SetDistances(const InterpolationMap& redshift_distance_map);
    
    
    
    
private:
    // map with the plates information
    PlatesMapVector<LyaSpectrum>::map list_;
    
    
};

#endif
