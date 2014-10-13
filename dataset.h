/**
 dataset.h
 Purpose: This file defines the class Dataset. This class contains the variables necessary to store a dataset
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 07/25/2014
 
 */

#ifndef _Dataset_h
#define _Dataset_h

// libraries needed
#include <sstream>
#include <string>
#include <vector>
////////

// classes needed
#include "astro_object.h"
#include "interpolation_map.h"
////////

// functions needed
////////

#include "typedefs.h"

class Dataset{
    
public:
    
    // -------------------------------------------------------------
    // access methods
        
    // access function for name_
    std::string name() const {return name_;}
    
    // access function for num_objects_in_plate_
    PlatesMapSimple<size_t>::map num_objects_in_plate() const {return num_objects_in_plate_;}
    size_t num_objects_in_plate(int plate_number) const {return (*num_objects_in_plate_.find(plate_number)).second;}
    
    // access function for size_
    size_t size() const {return size_;}
    
    // -------------------------------------------------------------
    // other methods
    
    // adds the objects' RA-DEC values to ostream
    virtual void GiveRADEC(std::ostream& out) const = 0;
    
    // adds the objects' redshift values to ostream
    virtual void GiveZ(std::ostream& out) const = 0;
    
    // set the distance to every object in the dataset
    virtual void SetDistances(const InterpolationMap& redshift_distance_map) = 0;
        

    
protected:
    // name of the catalog
    std::string name_;
    
    // number of objects in plate
    PlatesMapSimple<size_t>::map num_objects_in_plate_;
    
    // size of dataset
    size_t size_;
    


    
};

#endif
