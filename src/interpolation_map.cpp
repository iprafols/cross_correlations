/**
 interpolation_table.cpp
 Purpose: This files contains the body for the functions defined in interpolation_table.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 09/30/2014
 */

#include "interpolation_map.h"

double InterpolationMap::interpolation_map(double first) const{
    /**
     EXPLANATION:
     Access function for interpolation_map_
     
     INPUTS:
     it - pointer to the element to be recovered
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     InterpolationMap
     
     FUNCITONS USED:
     NONE
     */
    std::map<double,double>::const_iterator it = interpolation_map_.find(first);
    if (it == interpolation_map_.end()){
        return _BAD_DATA_;
    }
    else{
        return (*it).second;
    }
}

double InterpolationMap::interpolation_map(std::map<double,double>::iterator it) const{
    /**
     EXPLANATION:
     Access function for interpolation_map_
     
     INPUTS:
     it - pointer to the element to be recovered
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     InterpolationMap
     
     FUNCITONS USED:
     NONE
     */
    if (it == interpolation_map_.end()){
        return _BAD_DATA_;
    }
    else{
        return (*it).second;
    }
}

double InterpolationMap::LinearInterpolation(const double& first) const{
    /**
     EXPLANATION:
     Compute the value corresponding to the given redshift by using linear interpolation
     
     INPUTS:
     first - a double with the "x" value
     
     OUTPUTS:
     y - the interpolated "y" value
     
     CLASSES USED:
     InterpolationMap
     
     FUNCITONS USED:
     NONE
     */
    
    std::map<double,double>::const_iterator it,it2;

    it = interpolation_map_.lower_bound(first);
    
    if (it == interpolation_map_.end()){
        return _BAD_DATA_;
    }
    
    // no interpolation needed if distance is already computed for the given redshift
    if ((*it).first == first){
        return (*it).second;
    } else{
        it2 = interpolation_map_.lower_bound(first);
        it2 --;
        return ((*it).second-(*it2).second)/((*it).first-(*it2).first)*(first-(*it).first)+(*it).second;
    }
}