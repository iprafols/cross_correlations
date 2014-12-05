/**
 dla_dataset.cpp
 Purpose: This files contains the body for the functions defined in dla_dataset.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 11/25/2014
 */

#include "dla_dataset.h"

DLADataset::DLADataset(const Input& input){
    /**
     EXPLANATION:
     Cosntructs a DLADataset instance and loads a catalog of DLAs into it
     
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
    
    flag_verbose_dla_dataset_ = input.flag_verbose_dla_dataset();    
    name_ = input.dataset1_name();
    Load(input.z_min(), input.z_max(), input.dataset1());
    
}

void DLADataset::Load(const double& z_min, const double& z_max, const std::string& dataset1){
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
    size_t ra_index, dec_index, mpf_index, z_index, pos, pos2;
    int plate, fiber, mjd;
    double ra, dec, z;
    
    if (flag_verbose_dla_dataset_ >= 1){
        std::cout << "Loading DLA dataset" << std::endl;
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
                        else if (cols[i] == "MJD-plate-fiber"){
                            mpf_index = i;
                        }
                        else if (cols[i] == "z_abs"){
                            z_index = i;
                        }
                    }
                    read_column_names = false;
                }
                // read entries
                else{
                    ra = atof(cols[ra_index].c_str());
                    dec = atof(cols[dec_index].c_str());
                    z = atof(cols[z_index].c_str());
                    
                    mpf = cols[mpf_index];
                    pos = mpf.find("-");
                    pos2 = mpf.find("-",pos + 1);
                    
                    mjd = atoi(mpf.substr(0, pos).c_str());
                    plate = atoi(mpf.substr(pos+1, pos2).c_str());
                    fiber = atoi(mpf.substr(pos2+1).c_str());
                    
                    if (z > z_min and z < z_max){
                        // create AstroObject
                        AstroObject object(ra, dec, plate, fiber, mjd, z, false);
                        
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
                        
                        if (flag_verbose_dla_dataset_ >= 3 or (flag_verbose_dla_dataset_ >= 2 and size_ == size_/1000*1000)){
                            std::cout << "Loaded " << size_ << " DLAs" << std::endl;
                        }

                    }

                }
            }
        }
        if (flag_verbose_dla_dataset_ >= 1){
            std::cout << "Loaded " << size_ << " DLAs" << std::endl;
        }
        file.close();
    }
    else{
        std::cout << "Error : In DLADataset::Load : Unable to open file: " << std::endl << dataset1 << std::endl;
        std::exit;
    }
        
}