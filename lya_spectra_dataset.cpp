/**
 lya_spectra_dataset.cpp
 Purpose: This files contains the body for the functions defined in lya_spectra_dataset.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 08/08/2014
 */

#include "lya_spectra_dataset.h"

LyaSpectraDataset::LyaSpectraDataset(const Input& input){
    /**
     EXPLANATION:
     Cosntructs a LyaSpectraDataset instance and loads a catalog of LyaSpectrum objects into it
     
     INPUTS:
     input - object of type cInput
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaSpectraDataset
     Input
     
     FUNCITONS USED:
     NONE
     */
    
    name_ = input.lya_spectra_catalog_name();
    Load(input.lya_spectra_catalog(), input.lya_spectra_dir(), input.lya_wl());
    
}

int LyaSpectraDataset::FindCatalogLength(const std::string& lya_spectra_catalog){
    /**
     EXPLANATION:
     Finds the number of lines in the provided catalog file list
     
     INPUTS:
     lya_spectra_catalog - name of the list of fits files containing the spectra
     
     OUTPUTS:
     length - an integer containing the number of lines in the catalog file list
     
     CLASSES USED:
     NONE
     
     FUNCITONS USED:
     NONE
     */

    std::string cmd = "wc -l "+lya_spectra_catalog;
    char path[PATH_MAX];
    FILE* pipe = popen(cmd.c_str(),"r");
    
    fgets(path, PATH_MAX, pipe);
    
    int length = atoi(strtok(path," "));

    return length;
}

void LyaSpectraDataset::GiveRADEC(std::ostream& out) const{
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

void LyaSpectraDataset::GiveZ(std::ostream& out) const{
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

void LyaSpectraDataset::Load(const std::string& lya_spectra_catalog, const std::string& lya_spectra_dir, const double& lya_wl){
    /**
     EXPLANATION:
     Loads the object dataset from a catalog file
     
     INPUTS:
     lya_spectra_catalog - a string with the name of the list of fits files containing the spectra
     lya_spectra_dir - a string with the directory where the ly-a spectrum files are stored
     lya_wl - a double with the wavelength of the lyman-alpha line (in Angstroms)
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaSpectrum
     LyaSpectraDataset
     
     FUNCITONS USED:
     NONE
     */
    
    // resizing astro_object_pointer
    size_ = FindCatalogLength(lya_spectra_catalog);

    // open catalog
    std::ifstream catalog(lya_spectra_catalog.c_str());
    if (catalog.is_open()){
        std::string file("");        
        while (getline(catalog,file)){
            
            // create LyaSpectrum
            LyaSpectrum object(lya_spectra_dir + file, lya_wl);
            
            // adding object to list_
            if (list_.find(object.plate()) == list_.end()){
                // if necessary, create new entry
                std::vector<LyaSpectrum> v;
                list_[object.plate()] = v;
                num_objects_in_plate_[object.plate()] = 0;
            }
            (*list_.find(object.plate())).second.push_back(object);
            
            // updating number_of_objects_in_plate
            (*num_objects_in_plate_.find(object.plate())).second ++;
        }
    }
    else{
        std::cout << "Error: could not read spectra catalog" << std::endl;
    }
}

void LyaSpectraDataset::SetDistances(const InterpolationMap& redshift_distance_map){
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
