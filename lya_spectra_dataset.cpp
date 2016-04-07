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
    
    flag_verbose_lya_spectra_dataset_ = input.flag_verbose_lya_spectra_dataset();
    name_ = input.dataset2_name();
    Load(input.dataset2(), input.lya_spectra_dir(), input.lya_wl());
    if (input.flag_project_deltas()){
        ProjectDeltas();
    }
    
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
    
    if (flag_verbose_lya_spectra_dataset_ >= 1){
        std::cout << "Loading lya spectra dataset" << std::endl;
    }
    
    // resizing astro_object_pointer
    size_ = 0;

    // open catalog
    std::ifstream catalog(lya_spectra_catalog.c_str());
    if (catalog.is_open()){
        std::string file("");        
        while (getline(catalog,file)){
            
            size_t extension = 1;
            while (extension != 0){
                try{
                    // create LyaSpectrum
                    LyaSpectrum object(lya_spectra_dir + file, lya_wl, extension);
                    
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
                    
                    size_ ++;
                    if (flag_verbose_lya_spectra_dataset_ >= 3 or (flag_verbose_lya_spectra_dataset_ >= 2 and size_ == size_/1000*1000)){
                        std::cout << "Loaded " << size_ << " lya spectra" << std::endl;
                    }
                    extension ++;
                }
                catch (CCfits::FITS::NoSuchHDU x){
                    extension = 0;
                }
            }
        }
        
        if (flag_verbose_lya_spectra_dataset_ >= 1){
            std::cout << "Loaded " << size_ << " lya spectra" << std::endl;
        }
        catalog.close();
    }
    else{
        std::cout << "Error: in LyaSpectraDataset::Load : Could not read spectra catalog" << std::endl;
    }
}

void LyaSpectraDataset::ProjectDeltas(){
    /**
     EXPLANATION:
     Projects the delta field
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaSpectrum
     LyaSpectraDataset
     
     FUNCITONS USED:
     NONE
     */
    
    if (flag_verbose_lya_spectra_dataset_ >= 1){
        std::cout << "Projecting deltas" << std::endl;
    }
    
    // loop over plates
    size_t plates_computed = 0;
    for (PlatesMapVector<LyaSpectrum>::map::iterator it = list_.begin(); it != list_.end(); it ++){
        
        plates_computed ++;
        if (flag_verbose_lya_spectra_dataset_ >= 2 or (flag_verbose_lya_spectra_dataset_ >= 1 and plates_computed == plates_computed/100*100)){
            #pragma omp critical (cout)
            {
                std::cout << plates_computed << " out of " << list_.size() << " plates computed" << std::endl;
            }
        }

        // loop over spectra
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < (*it).second.size(); i ++){
            (*it).second[i].ProjectDeltas();
        }
    }
    
}


