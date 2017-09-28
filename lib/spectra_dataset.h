/**
 spectra_dataset.h
 Purpose: This file defines the class SpectraDataset. This class contains the variables necessary to store a generic spectra dataset. This class is a specialization of the Dataset class
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 10/02/2014
 
 */

#ifndef _SpectraDataset_h
#define _SpectraDataset_h

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


class SpectraDataset: public Dataset{
    
public:
    
    // -------------------------------------------------------------
    // access methods
    
    // access functions for list_
    PlatesMapVector<LyaSpectrum>::map list() const {return list_;}
    std::vector<LyaSpectrum> list(int plate_number) const;
    LyaSpectrum list(int plate_number, size_t pos) const;
    
    // -------------------------------------------------------------
    // other methods
    
    // adds the objects' RA-DEC values to ostream
    void GiveRADEC(std::ostream& out) const;
    
    // adds the objects' redshift values to ostream
    void GiveZ(std::ostream& out) const;
    
    // load the dataset
    virtual void Load(const std::string& object_list, const std::string& lya_spectra_dir, const double& lya_wl, const std::vector<double>& alt_wl) =0;
    
    // set the distance to every object in the dataset
    void SetDistances(const ZDistInterpolationMap& redshift_distance_map);
    
protected:
    // map with the plates information
    PlatesMapVector<LyaSpectrum>::map list_;
    
    
};

#endif
