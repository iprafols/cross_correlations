/**
 strong_lya_dataset.cpp
 Purpose: This files contains the body for the functions defined in strong_lya_dataset.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 12/16/2016
 */

#include "strong_lya_dataset.h"

StrongLyaDataset::StrongLyaDataset(const Input& input){
    /**
     EXPLANATION:
     Cosntructs a StrongLyaDataset instance and loads a catalog of strong lyman alpha absorbers into it
     
     INPUTS:
     input - object of type Input
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     DLADataset
     Input
     
     FUNCITONS USED:
     NONE
     */
    
    flag_verbose_strong_lya_dataset_ = input.flag_verbose_strong_lya_dataset();
    name_ = input.dataset1_name();
    lya_flux_min_ = input.lya_flux_min();
    lya_flux_max_ = input.lya_flux_max();
    Load(input.z_min(), input.z_max(), input.dataset1());
    
    
}

void StrongLyaDataset::Load(const double& z_min, const double& z_max, const std::string& dataset1){
    /**
     EXPLANATION:
     Loads the object dataset from a catalog file
     
     INPUTS:
     z_min - minimum redshift to accept AstroObject into the dataset
     z_max - maximum redshift to accept AstroObject into the dataset
     dataset1 - name of the catalog file
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     AstroObjectDataset
     
     FUNCITONS USED:
     NONE
     */
    std::string line, mpf;
    bool read_column_names;
    size_t ra_index, dec_index, mjd_index, plate_index, fiber_index, z_abs_index, z_qso_index, lya_flux_index;
    int plate, fiber, mjd;
    double ra, dec, z_abs, z_qso, nhi, bi, cnr, rf_wl;
    
    if (flag_verbose_dla_dataset_ >= 1){
        std::cout << "Loading StrongLya dataset" << std::endl;
    }
    
    // setting size to zero
    size_ = 0;
    
    std::ifstream file(dataset1.c_str());
    if (file.is_open()){
        read_column_names = true; 
        // reading catalog
        while (getline(file,line)){
            if (line[0] != '#' and line[0] != '-'){
                std::vector<std::string> cols = Split(line," ");
                
                // get the column number for each variable
                if (read_column_names){
                    for (size_t i = 0; i < cols.size(); i++){
                        if (cols[i] == "RA"){
                            ra_index = i;
                        }
                        else if (cols[i] == "Dec"){
                            dec_index = i;
                        }
                        else if (cols[i] == "MJD"){
                            mjd_index = i;
                        }
                        else if (cols[i] == "plate"){
                            plate_index = i;
                        }
                        else if (cols[i] == "fiber"){
                            fiber_index = i;
                        }
                        else if (cols[i] == "z_abs"){
                            z_abs_index = i;
                        }
                        else if (cols[i] == "zqso"){
                            z_qso_index = i;
                        }
                        else if (cols[i] == "lya_flux"){
                            lya_flux_index = i;
                        }
                    }
                    read_column_names = false;
                }
                // read entries
                else{
                    ra = atof(cols[ra_index].c_str());
                    dec = atof(cols[dec_index].c_str());
                    z_abs = atof(cols[z_abs_index].c_str());
                    z_qso = atof(cols[z_qso_index].c_str());
                    mjd = atoi(cols[mjd_index].c_str());
                    plate = atoi(cols[plate_index].c_str());
                    fiber = atoi(cols[fiber_index].c_str());
                    lya_flux = atof(cols[lya_flux_index].c_str());
                    
                    if (z_abs >= z_min and z_abs < z_max and lya_flux >= lya_flux_min_ and lya_flux < lya_flux_max_){
                        
                        // create AstroObject
                        AstroObject object(ra, dec, plate, fiber, mjd, z_abs, false);
                        
                        // adding object to list_
                        if (list_.find(plate) == list_.end()){
                            // if necessary, create new entry
                            std::vector<AstroObject> v;
                            list_[plate] = v;
                            num_objects_in_plate_[plate] = 0;
                        }
                        (*list_.find(plate)).second.push_back(object);
                        
                        // updating size_
                        size_ ++;
                        
                        // updating number_of_objects_in_plate
                        (*num_objects_in_plate_.find(plate)).second ++;
                        
                        if (flag_verbose_strong_lya_dataset_ >= 3 or (flag_verbose_strong_lya_dataset_ >= 2 and size_ == size_/1000*1000)){
                            std::cout << "Loaded " << size_ << " strong lyman alpha absorbers" << std::endl;
                        }

                    }

                }
            }
        }
        if (flag_verbose_strong_lya_dataset_ >= 1){
            std::cout << "Loaded " << size_ << " strong lyman alpha absorbers" << std::endl;
        }
        file.close();
    }
    else{
        std::cout << "Error : In StrongLyaDataset::Load : Unable to open file: " << std::endl << dataset1 << std::endl;
        std::exit(EXIT_FAILURE);
    }
        
}