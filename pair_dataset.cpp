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
    char path[PATH_MAX];
    FILE* pipe = popen(cmd.c_str(),"r");
    
    fgets(path, PATH_MAX, pipe);
    
    int length = atoi(strtok(path," "));
    
    return length - 1;
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
    
    if (flag_verbose_pair_dataset_ >= 1){
        std::cout << "Loading " << name_ << std::endl;
    }
    
    // setting size to zero
    size_ = 0;
    
    // loop over plates
    for (size_t i = 0; i < plates.size(); i++){
        
        // checking if the selected plate has any contributions to this bin
        filename = pairs_file_name_ + ToStr(plates[i]);
        std::ifstream file(filename.c_str());
        if (file.is_open()){
            
            if (list_.find(plates[i]) == list_.end()){
                std::vector<Pair> v;
                num_objects_in_plate_[plates[i]] = FindCatalogLength(filename);
                v.reserve(num_objects_in_plate_[plates[i]]);
                list_[plates[i]] = v;
            }
            
            read_column_names = true;
            // reading catalog
            while (getline(file, line)){
                if (line[0] != '#'){
                    std::vector<std::string> cols = Split(line," ");
                    
                    // get the column number for each variable
                    if (read_column_names){
                        for (size_t i = 0; i < cols.size(); i++){
                            if (cols[i] == "spectrum_RA"){
                                spectrum_ra_index = i;
                            }
                            else if (cols[i] == "spectrum_DEC"){
                                spectrum_dec_index = i;
                            }
                            else if (cols[i] == "pixel_number"){
                                pixel_number_index = i;
                            }
                            else if (cols[i] == "pixel_dist"){
                                pixel_dist_index = i;
                            }
                            else if (cols[i] == "pixel_w"){
                                pixel_weight_index = i;
                            }
                        }
                        read_column_names = false;
                    }
                    
                    // read entries
                    else{
                        spectrum_ra = atof(cols[spectrum_ra_index].c_str());
                        spectrum_dec = atof(cols[spectrum_dec_index].c_str());
                        pixel_number = atoi(cols[pixel_number_index].c_str());
                        pixel_dist = atof(cols[pixel_dist_index].c_str());
                        pixel_weight = atof(cols[pixel_weight_index].c_str());
                        
                        // create Pair
                        Pair object(spectrum_ra, spectrum_dec, pixel_number, pixel_dist, pixel_weight);
                        
                        // adding Pair to list
                        (*list_.find(plates[i])).second.push_back(object);
                        
                        // updating size_
                        size_ ++;
                        
                        if (flag_verbose_pair_dataset_ >= 3 or (flag_verbose_pair_dataset_ >= 2 and size_ == size_/1000*1000)){
                            std::cout << "Loaded " << size_ << " pairs" << std::endl;
                        }
                    }
                }
            }
        }
    }
    
    if (flag_verbose_pair_dataset_ >= 1){
        std::cout << "Loaded " << size_ << " pairs" << std::endl;
    }
}

