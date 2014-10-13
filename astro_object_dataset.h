/**
 dataset.h
 Purpose: This file defines the class AstroObjectDataset. This class contains the variables necessary to store a object dataset. This class is a specialization of the cDataset class
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 07/25/2014
 
 */

#ifndef _AstroObjectDataset_h
#define _AstroObjectDataset_h

// libraries needed
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <CCfits>
////////

// classes needed
#include "astro_object.h"
#include "dataset.h"
#include "global_variables.h"
#include "interpolation_map.h"
////////

// functions needed
////////

#include "typedefs.h"


class AstroObjectDataset: public Dataset{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs empty object
    AstroObjectDataset(){};
    
    // constructs object and initializes its variables
    AstroObjectDataset(const GlobalVariables& kGlobalVariables);
    
    // -------------------------------------------------------------
    // access methods
    
    // access functions for list_
    PlatesMapVector<AstroObject>::map list() const {return list_;}
    std::vector<AstroObject> list(int plate_number) const {return (*list_.find(plate_number)).second;}
    AstroObject list(int plate_number, size_t pos) const {return (*list_.find(plate_number)).second[pos];}
    
    // -------------------------------------------------------------
    // other methods
    
    // adds the objects' RA-DEC values to ostream
    void GiveRADEC(std::ostream& out) const;
    
    // adds the objects' redshift values to ostream
    void GiveZ(std::ostream& out) const;
    
    // load the dataset
    void Load(const double& z_min, const double& z_max, const std::string& object_list);
    
    // set the distance to every object in the dataset
    void SetDistances(const InterpolationMap& redshift_distance_map);

    
    
    
    
    
private:        
    // map with the plates information
    PlatesMapVector<AstroObject>::map list_;
    

    
};

#endif
