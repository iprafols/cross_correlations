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
    nhi_min_ = input.nhi_min();
    nhi_max_ = input.nhi_max();
    cnr_min_ = input.cnr_min();
    rf_wl_max_ = input.rf_wl_max();
    rf_wl_min_ = input.rf_wl_min();
    rf_wl_forbidden_interval_ = input.rf_wl_forbidden_interval();
    lya_wl_ = input.lya_wl();
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
    size_t ra_index, dec_index, mpf_index, z_abs_index, z_qso_index, nhi_index, bi_index, cnr_index, pos, pos2;
    int plate, fiber, mjd;
    double ra, dec, z_abs, z_qso, nhi, bi, cnr, rf_wl;
    
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
                            z_abs_index = i;
                        }
                        else if (cols[i] == "zqso"){
                            z_qso_index = i;
                        }
                        else if (cols[i] == "NHI"){
                            nhi_index = i;
                        }
                        else if (cols[i] == "BI"){
                            bi_index = i;
                        }
                        else if (cols[i] == "CNR"){
                            cnr_index = i;
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
                    nhi = atof(cols[nhi_index].c_str());
                    bi = atof(cols[bi_index].c_str());
                    cnr = atof(cols[cnr_index].c_str());
                    rf_wl = lya_wl_*(1.0-(z_qso-z_abs)/(1.0+z_qso));
                    
                    mpf = cols[mpf_index];
                    pos = mpf.find("-");
                    pos2 = mpf.find("-",pos + 1);
                    
                    mjd = atoi(mpf.substr(0, pos).c_str());
                    plate = atoi(mpf.substr(pos+1, pos2).c_str());
                    fiber = atoi(mpf.substr(pos2+1).c_str());
                    
                    if (z_abs >= z_min and z_abs < z_max and nhi >= nhi_min_ and nhi < nhi_max_ and cnr >= cnr_min_ and bi == 0.0 and rf_wl <= rf_wl_max_ and rf_wl >= rf_wl_min_ and (rf_wl <= rf_wl_forbidden_interval_.first or rf_wl >= rf_wl_forbidden_interval_.second)){
                        
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
        std::exit(EXIT_FAILURE);
    }
        
}