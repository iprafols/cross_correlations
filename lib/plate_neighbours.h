/**
 plate_neighbours.h
 Purpose: This file defines the class PlateNeighbours. This class contains the plates list and the corresponding neighbouring plates
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 06/26/2014
 
 */

#ifndef _PlateNeighbours_h
#define _PlateNeighbours_h

// libraries needed
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <vector>
////////

// classes needed
#include "input.h"
#include "plate.h"
////////

// functions needed
#include "function_split.hpp"
////////

#include "typedefs.h"
#include "defines.h"

class PlateNeighbours{
    
public:
    // -------------------------------------------------------------
    // constructors

    // constructs empty object
    PlateNeighbours(){};
    
    // constructs object and initializes its variables
    PlateNeighbours(const Input& Input);
    
    // -------------------------------------------------------------
    // access methods
    
    // access function for flag_verbose_plate_neighbours_
    size_t  flag_verbose_plate_neighbours() const {return flag_verbose_plate_neighbours_;}
    
    // access function for plates_
    PlatesMapVector<Plate>::map plates() const {return plates_;}
    
    // return list of neighbours
    std::vector<int> GetNeighboursList(int plate) const;
    std::vector<int> GetNeighboursList(const Plate& plate) const {return GetNeighboursList(plate.plate_number());}
    
    // return list of plates
    std::vector<int> GetPlatesList() const;

    // -------------------------------------------------------------
    // other methods

    // add neighbours
    void AddNeighbours(const Plate& plate1, const Plate& plate2);
    
    // add plate entry
    void AddPlate(const Plate& plate);
    
    // check whether or not two plates are neighbours
    bool AreNeighbours(const int plate1, const int plate2) const;
    bool AreNeighbours(const Plate& plate1, const int plate2) const {return AreNeighbours(plate1.plate_number(), plate2);}
    bool AreNeighbours(const int plate1, const Plate& plate2) const {return AreNeighbours(plate1, plate2.plate_number());}
    bool AreNeighbours(const Plate& plate1, const Plate& plate2) const {return AreNeighbours(plate1.plate_number(), plate2.plate_number());}
   
    // save list of plates into file
    void Save(const std::string& filename);

    
    

    
private:
    // plate_neighours verbose flag
    size_t flag_verbose_plate_neighbours_;
    
    // list of plates with the corresponding neighbours
    PlatesMapVector<Plate>::map plates_;

    // -------------------------------------------------------------
    // other methods
    void ReadPlateNeighbours(std::ifstream& plates_file);
};

#endif
