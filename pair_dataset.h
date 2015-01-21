/**
 piars_dataset.h
 Purpose: This file defines the class PairsDataset. This class contains the variables necessary to store the information of a set of object-pixel pairs.
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 01/20/2015
 
 */

#ifndef _PairsDataset_h
#define _PairsDataset_h

// libraries needed
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
////////

// classes needed
#include "dataset.h"
#include "z_dist_interpolation_map.h"
#include "pair.h"
#include "plate_neighbours.h"
////////

// functions needed
#include "function_split.hpp"
#include "function_to_str.hpp"
////////

#include "defines.h"
#include "typedefs.h"


class PairDataset: public Dataset{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs empty object
     PairDataset(){};
     
     // constructs object and initializes its variables
     PairDataset(const Input& input, const size_t bin, const std::vector<int>& plates);
    
    // -------------------------------------------------------------
    // access methods
    
    // access function for bin_
    size_t bin() const {return bin_;}
    
    // access function for flag_verbose_pair_dataset_
    size_t flag_verbose_pair_dataset() const {return flag_verbose_pair_dataset_;}
    // access functions for list_
    PlatesMapVector<Pair>::map list() const {return list_;}
    std::vector<Pair> list(int plate_number) const;
    Pair list(int plate_number, size_t pos) const;
    
    // access function for pairs_file_name_
    std::string pairs_file_name() const {return pairs_file_name_;}
    
    // -------------------------------------------------------------
    // other methods
    
    // adds the objects' RA-DEC values to ostream (disabled)
    void GiveRADEC(std::ostream& out) const {};
    
    // adds the objects' redshift values to ostream (disabled)
    void GiveZ(std::ostream& out) const {};
    
    // load the dataset
    void Load(const std::vector<int>& plates);
    
    // set the distance to every object in the dataset (disabled)
    void SetDistances(const ZDistInterpolationMap& redshift_distance_map){};
    
    // -------------------------------------------------------------
    // other methods
    
    // return list of plates
    std::vector<int> GetPlatesList() const;
    
    // return number of pairs in plate entry
    size_t GetNumberPairs(int plate) const;
    
    
    
private:
    // bin number
    size_t bin_;
    
    // pair_dataset verbose flag
    size_t flag_verbose_pair_dataset_;
    
    // map with the plates information
    PlatesMapVector<Pair>::map list_;
    
    // basename of the files where the pairs' information is stored
    std::string pairs_file_name_;
    
    

    
    // -------------------------------------------------------------
    // other methods
    
    // find the catalog length
    int FindCatalogLength(const std::string& filename);
    

    
};

#endif
