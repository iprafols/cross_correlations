/**
 spectra_dataset.cpp
 Purpose: This files contains the body for the functions defined in spectra_dataset.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 10/02/2015
 */

#include "spectra_dataset.h"

std::vector<LyaSpectrum> SpectraDataset::list(int plate_number) const{
    /**
     EXPLANATION:
     Access function for list_
     
     INPUTS:
     plate_number - index of the selected list_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaSpectrum
     LyaSpectraDataset
     
     FUNCITONS USED:
     NONE
     */
    
    PlatesMapVector<LyaSpectrum>::map::const_iterator it = list_.find(plate_number);
    if (it == list_.end()){
        std::vector<LyaSpectrum> v;
        return v;
    }
    else{
        return (*it).second;
    }
}

LyaSpectrum SpectraDataset::list(int plate_number, size_t pos) const {
    /**
     EXPLANATION:
     Access function for list_
     
     INPUTS:
     plate_number - index of the selected list_ element
     pos - position in the selected list_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaSpectrum
     LyaSpectraDataset
     
     FUNCITONS USED:
     NONE
     */
    
    PlatesMapVector<LyaSpectrum>::map::const_iterator it = list_.find(plate_number);
    if (it == list_.end()){
        LyaSpectrum v(_BAD_DATA_);
        return v;
    }
    else{
        if (pos < (*it).second.size()){
            return (*it).second[pos];
        }
        else{
            LyaSpectrum v(_BAD_DATA_);
            return v;
        }
    }
}

void SpectraDataset::GiveRADEC(std::ostream& out) const{
    /**
     EXPLANATION:
     Adds the objects' RA-DEC values to out
     
     INPUTS:
     out - an ostream to add the objects' RA-DEC values to
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaSpectrum
     LyaSpectraDataset
     
     FUNCITONS USED:
     NONE
     */
    
    for (PlatesMapVector<LyaSpectrum>::map::const_iterator it = list_.begin(); it != list_.end(); it ++){
        for (size_t i = 0; i < (*it).second.size(); i ++){
            out << (*it).second[i].angle() << std::endl;
        }
    }
    
}

void SpectraDataset::GiveZ(std::ostream& out) const{
    /**
     EXPLANATION:
     Adds the objects' redshift values to out
     
     INPUTS:
     out - an ostream to add the objects' redshift values to
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaPixel
     LyaSpectrum
     LyaSpectraDataset
     
     FUNCITONS USED:
     NONE
     */
    
    for (PlatesMapVector<LyaSpectrum>::map::const_iterator it = list_.begin(); it != list_.end(); it ++){
        for (size_t i = 0; i < (*it).second.size(); i ++){
            for (size_t j = 0; j < (*it).second[i].SpectrumSize(); j++){
                out << (*it).second[i].spectrum(j).z() << std::endl;
            }
        }
    }
    
}

void SpectraDataset::SetDistances(const ZDistInterpolationMap& redshift_distance_map){
    /**
     EXPLANATION:
     Sets the distance to every object in the dataset
     
     INPUTS:
     redshif_distance_map - a InterpolationMap instance with the redshift-distance relation
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaSpectrum
     LyaSpectraDataset
     InterpolationMap
     
     FUNCITONS USED:
     NONE
     */
    
    for (PlatesMapVector<LyaSpectrum>::map::iterator it = list_.begin(); it != list_.end(); it ++){
        for (size_t i = 0; i < (*it).second.size(); i ++){
            (*it).second[i].SetDistance(redshift_distance_map);
        }
    }
    
}
