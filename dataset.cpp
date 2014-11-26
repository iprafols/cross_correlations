/**
 dataset.cpp
 Purpose: This files contains the body for the functions defined in dataset.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 07/25/2014
 */

#include "dataset.h"

size_t Dataset::num_objects_in_plate(int plate_number) const{
    /**
     EXPLANATION:
     Access function for num_objects_in_plate_
     
     INPUTS:
     plate_number - indez of the selected num_objects_in_plate_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Dataset
     
     FUNCITONS USED:
     NONE
     */
    
    PlatesMapSimple<size_t>::map::const_iterator it = num_objects_in_plate_.find(plate_number);
    
    if (it == num_objects_in_plate_.end()){
        return 0;
    }
    else{
        return (*it).second;
    }
    
}