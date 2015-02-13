/**
 pairs_dataset.cpp
 Purpose: This files contains the body for the functions defined in pairs_dataset.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 01/20/2015
 */

#include "pair_dataset.h"

PairDataset::PairDataset(const Input& input, const size_t bin, const std::vector<int>& plates){
    /**
     EXPLANATION:
     Cosntructs a PairDataset instance and loads a list of pairs into it
     
     INPUTS:
     input - object of type Input
     bin - an unsigned integer with the bin number
     plates -  a vector containing a list of plates
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Input
     PairDataset
     
     FUNCITONS USED:
     ToStr
     */
    
    flag_verbose_pair_dataset_ = input.flag_verbose_pair_dataset();
    name_ = "pairs in bin " + ToStr(bin);
    bin_ = bin;
    pairs_file_name_ = input.detailed_results() + ToStr(bin) + "/" + input.pairs_file_name();
    
    Load(plates);
    
}

std::vector<Pair> PairDataset::list(int plate_number) const{
    /**
     EXPLANATION:
     Access function for list_
     
     INPUTS:
     plate_number - index of the selected list_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Pair
     PairDataset
     
     FUNCITONS USED:
     NONE
     */
    
    PlatesMapVector<Pair>::map::const_iterator it = list_.find(plate_number);
    if (it == list_.end()){
        std::vector<Pair> v;
        return v;
    }
    else{
        return (*it).second;
    }

}

Pair PairDataset::list(int plate_number, size_t pos) const{
    /**
     EXPLANATION:
     Access function for list_
     
     INPUTS:
     plate_number - index of the selected list_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Pair
     PairDataset
     
     FUNCITONS USED:
     NONE
     */
    
    PlatesMapVector<Pair>::map::const_iterator it = list_.find(plate_number);
    if (it == list_.end()){
        Pair v(_BAD_DATA_);
        return v;
    }
    else{
        if (pos < (*it).second.size()){
            return (*it).second[pos];
        }
        else{
            Pair v(_BAD_DATA_);
            return v;
        }
    }
}

int PairDataset::FindCatalogLength(const std::string& filename){
    /**
     EXPLANATION:
     Finds the number of lines in the provided catalog file list
     
     INPUTS:
     filename - name of the file containing the pairs list
     
     OUTPUTS:
     length - an integer containing the number of lines in the catalog file list
     
     CLASSES USED:
     NONE
     
     FUNCITONS USED:
     NONE
     */
    
    std::string cmd = "wc -l " + filename;
    if (flag_verbose_pair_dataset_ >= 2){
        std::cout << cmd << std::endl;
    }
    std::cout << "pas 0; ";
    char path[10];
    std::cout << "pas 1; ";
    FILE* pipe = popen(cmd.c_str(),"r");
    
    std::cout << "pas 2; ";
    fgets(path, 10, pipe);
    std::cout << "pas 3; ";
    int length = atoi(strtok(path," "));
    std::cout << "pas 4; length = " << length << std::endl;
    
    return length - 1;
}

std::vector<int> PairDataset::GetPlatesList() const {
    /**
     EXPLANATION:
     Returns a vector conatining the plate number of all the used plates
     
     INPUTS:
     NONE
     
     OUTPUTS:
     plates - a vector containg the list of used plates
     
     CLASSES USED:
     PLateNeighbours
     
     FUNCITONS USED:
     NONE
     */
    
    std::vector<int> plates_list;
    plates_list.reserve(list_.size());
    
    for (PlatesMapVector<Pair>::map::const_iterator it = list_.begin(); it != list_.end(); it ++){
        plates_list.push_back((*it).first);
    }
    
    return plates_list;
}


void PairDataset::Load(const std::vector<int>& plates){
    /**
     EXPLANATION:
     Loads the list of pairs from file
     
     INPUTS:
     plates - a vector containing the list of plates used in the analysis
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Pair
     PairDataset
     
     FUNCITONS USED:
     NONE
     */
    bool read_column_names;
    size_t spectrum_ra_index, spectrum_dec_index, pixel_number_index, pixel_dist_index, pixel_weight_index;
    double spectrum_ra, spectrum_dec, pixel_dist, pixel_weight;
    int pixel_number;
    std::string filename, line;
    
    if (flag_verbose_pair_dataset_ >= 2){
        std::cout << "Loading " << name_ << std::endl;
    }
    
    // setting size to zero
    size_ = 0;
    
    // loop over plates
    for (size_t i = 0; i < plates.size(); i++){
        
        // checking if the selected plate has any contributions to this bin
        if (flag_verbose_pair_dataset_ >= 3){
            std::cout << "checking plate " << plates[i] << std::endl;
        }
        filename = pairs_file_name_ + ToStr(plates[i]);
        std::ifstream file(filename.c_str());
        if (file.is_open()){
            
            if (list_.find(plates[i]) == list_.end()){
                std::vector<Pair> v;
                list_[plates[i]] = v;
                /*int num_objects = FindCatalogLength(filename);
                if (flag_verbose_pair_dataset_ >= 2){
                    std::cout << "There are " << num_objects << " pairs in this plate" << std::endl;
                }
                /*if (num_objects > 0){
                    num_objects_in_plate_[plates[i]] = num_objects;
                    v.reserve(num_objects);
                    list_[plates[i]] = v;
                }
                else{
                    continue;
                }*/
            }
            else {
                if (flag_verbose_pair_dataset_ >= 1){
                    std::cout << "Warning : In PairDataset::Load : This plate has already been checked out, ignoring..." << std::endl;
                }
                continue;
            }
            
            
            PlatesMapVector<Pair>::map::iterator list_it = list_.find(plates[i]);
            if (list_it == list_.end()){
                if (flag_verbose_pair_dataset_ >= 1){
                    std::cout << "Warning : In PairDataset::Load : Plate entry was not properly created, ignoring..." << std::endl;
                }
                continue;
            }
            read_column_names = true;
            // reading catalog
            while (getline(file, line)){
                std::vector<std::string> cols = Split(line," ");
                    
                // get the column number for each variable
                if (read_column_names){
                    if (flag_verbose_pair_dataset_ >= 4){
                        std::cout << "Reading header" << std::endl;
                    }
                    for (size_t i = 0; i < cols.size(); i++){
                        if (cols[i] == "spectrum_RA"){
                            spectrum_ra_index = i - 1;
                        }
                        else if (cols[i] == "spectrum_DEC"){
                            spectrum_dec_index = i - 1;
                        }
                        else if (cols[i] == "pixel_number"){
                            pixel_number_index = i - 1;
                        }
                        else if (cols[i] == "pixel_dist"){
                            pixel_dist_index = i - 1;
                        }
                        else if (cols[i] == "pixel_w"){
                            pixel_weight_index = i - 1;
                        }
                    }
                    read_column_names = false;
                }
                
                // read entries
                else{
                    if (flag_verbose_pair_dataset_ >= 4){
                        std::cout << "Reading entry" << std::endl;
                    }
                    spectrum_ra = atof(cols[spectrum_ra_index].c_str());
                    spectrum_dec = atof(cols[spectrum_dec_index].c_str());
                    pixel_number = atoi(cols[pixel_number_index].c_str());
                    pixel_dist = atof(cols[pixel_dist_index].c_str());
                    pixel_weight = atof(cols[pixel_weight_index].c_str());
                    
                    // create Pair
                    Pair object(spectrum_ra, spectrum_dec, pixel_number, pixel_dist, pixel_weight);
                        
                    // adding Pair to list
                    (*list_it).second.push_back(object);
                    
                    // updating size_
                    size_ ++;
                    
                    if (flag_verbose_pair_dataset_ >= 4 or (flag_verbose_pair_dataset_ >= 2 and size_ == size_/100000*100000)){
                        std::cout << "Loaded " << size_ << " pairs" << std::endl;
                    }
                }
            }
            file.close();
        }
        else if (flag_verbose_pair_dataset_ >= 2){
            std::cout << "There are no pairs found in this plate " << std::endl;
        }
    }
    
    if (flag_verbose_pair_dataset_ >= 2){
        std::cout << "Loaded " << size_ << " pairs" << std::endl;
    }
}

size_t PairDataset::GetNumberPairs(int plate) const{
    /**
     EXPLANATION:
     Returns the number of pairs contained in the specified plate
     
     INPUTS:
     plate - an integer representing a plate number
     
     OUTPUTS:
     num - number of pairs contained in the specified plate
     
     CLASSES USED:
     PairDataset
     
     FUNCITONS USED:
     NONE
     */
    PlatesMapVector<Pair>::map::const_iterator it = list_.find(plate);
    
    if (it == list_.end()){
        return 0;
    }
    else{
        return (*it).second.size();
    }
}