/**
 astro_object_dataset.cpp
 Purpose: This files contains the body for the functions defined in astro_object_dataset.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 08/08/2014
 */

#include "astro_object_dataset.h"

std::vector<AstroObject> AstroObjectDataset::list(int plate_number) const{
    /**
     EXPLANATION:
     Access function for list_
     
     INPUTS:
     plate_number - index of the selected list_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     AstroObjectDataset
     
     FUNCITONS USED:
     NONE
     */
    
    PlatesMapVector<AstroObject>::map::const_iterator it = list_.find(plate_number);
    if (it == list_.end()){
        std::vector<AstroObject> v;
        return v;
    }
    else{
        return (*it).second;
    }    
}

AstroObject AstroObjectDataset::list(int plate_number, size_t pos) const {
    /**
     EXPLANATION:
     Access function for list_
     
     INPUTS:
     plate_number - index of the selected list_ element
     pos - position in the selected list_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     AstroObjectDataset
     
     FUNCITONS USED:
     NONE
     */
    
    PlatesMapVector<AstroObject>::map::const_iterator it = list_.find(plate_number);
    if (it == list_.end()){
        AstroObject v(_BAD_DATA_);
        return v;
    }
    else{
        if (pos < (*it).second.size()){
            return (*it).second[pos];
        }
        else{
            AstroObject v(_BAD_DATA_);
            return v;
        }
    }  
}

void AstroObjectDataset::GiveRADEC(std::ostream& out) const{
    /**
     EXPLANATION:
     Adds the objects' RA-DEC values to out
     
     INPUTS:
     out - an ostream to add the objects' RA-DEC values to
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     AstroObjectDataset
     
     FUNCITONS USED:
     NONE
     */
    
    for (PlatesMapVector<AstroObject>::map::const_iterator it = list_.begin(); it != list_.end(); it ++){
        for (size_t i = 0; i < (*it).second.size(); i ++){
            out << (*it).second[i].angle() << std::endl;
        }
    }
    
}

void AstroObjectDataset::GiveZ(std::ostream& out) const{
    /**
     EXPLANATION:
     Adds the objects' redshift values to out
     
     INPUTS:
     out - an ostream to add the objects' redshift values to
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     AstroObjectDataset
     
     FUNCITONS USED:
     NONE
     */
    
    for (PlatesMapVector<AstroObject>::map::const_iterator it = list_.begin(); it != list_.end(); it ++){
        for (size_t i = 0; i < (*it).second.size(); i ++){
            out << (*it).second[i].z() << std::endl;
        }
    }
    
}

void AstroObjectDataset::SetDistances(const InterpolationMap& redshift_distance_map){
    /**
     EXPLANATION:
     Sets the distance to every object in the dataset
     
     INPUTS:
     redshif_distance_map - a InterpolationMap instance with the redshift-distance relation
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     AstroObjectDataset
     InterpolationMap
     
     FUNCITONS USED:
     NONE
     */

    for (PlatesMapVector<AstroObject>::map::iterator it = list_.begin(); it != list_.end(); it ++){
        for (size_t i = 0; i < (*it).second.size(); i ++){
            (*it).second[i].SetDistance(redshift_distance_map);
        }
    }
        
}

