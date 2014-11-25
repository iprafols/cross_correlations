/**
 correlation_results.cpp
 Purpose: This files contains the body for the functions defined in correlation_results.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 06/17/2014
 */

#include "correlation_results.h"


CorrelationResults::CorrelationResults(const Input& input, const PlateNeighbours& kPlateNeighbours){
    /**
     EXPLANATION:
     Cosntructs a CorrelationResults instance and initializes all its variables
     
     INPUTS:
     input - object of type Input
     kPlateNeighbours - a PlateNeighbours instance
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     CorrelationResults
     Input
     
     FUNCITONS USED:
     NONE
     */
    
    // setting the number of bins and plates from input
    num_bins_ = input.num_bins();
    
    // setting the results directory and the pairs file name from input
    results_ = input.results() + "detailed_info_bin_";
    pairs_file_name_ = "detailed_info_plate_";
    output_base_name_ = input.output() + input.output_base_name();
    
    // initialization of the plates map
    plates_list_ = kPlateNeighbours.GetPlatesList();
    for (size_t i = 0; i < plates_list_.size(); i ++){
        correlation_plates_[plates_list_[i]] = CorrelationPlate(plates_list_[i], num_bins_, results_, pairs_file_name_, kPlateNeighbours.GetNeighboursList(plates_list_[i]));
    }
    
    // initialization of the normalized cross-correlation variable
    normalized_correlation_ = CorrelationPlate(_NORM_, num_bins_, results_, "", kPlateNeighbours.GetNeighboursList(_NORM_));
    
    // initialization of the bootstrap variable
    flag_compute_bootstrap_ = input.flag_compute_bootstrap();
    if (flag_compute_bootstrap_){
        bootstrap_.reserve(input.num_bootstrap());
        for (size_t i = 0; i < input.num_bootstrap(); i++){
            bootstrap_.push_back(CorrelationPlate(_NORM_, num_bins_, results_, "error", kPlateNeighbours.GetNeighboursList(_NORM_)));
        }

    }
    
    // creating bin files
    if (!input.flag_load_only()){
        CreateBinFiles();
    }
}

void CorrelationResults::ComputeCrossCorrelation(const AstroObjectDataset& object_list, const LyaSpectraDataset& spectra_list, const Input& input){
    /**
     EXPLANATION:
     Computes the cross-correlation for all plates
     
     INPUTS:
     object_list - an AstroObjectDataset instance
     spectra_list - a LyaSpectraDataset instance
     input - a Input instance to load the bin settings
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObjectDataset
     CorrelationResults
     Input
     LyaSpectraDataset
     
     FUNCITONS USED:
     NONE
     */
    
    for (size_t i = 0; i < plates_list_.size(); i++){
        std::cout << i << " out of " << plates_list_.size() << " computed" << std::endl;
        (*correlation_plates_.find(plates_list_[i])).second.ComputeCrossCorrelation(object_list, spectra_list, input);
    }
    
    NormalizeCrossCorrelation();
    
    if (flag_compute_bootstrap_){
        ComputeBootstrapRealizations();
    }
    SaveCrossCorrelation();
    
}

void CorrelationResults::ComputeBootstrapRealizations(){
    /**
     EXPLANATION:
     Computes and normalizes the bootstrap realizations
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     CorrelationResults
     
     FUNCITONS USED:
     NONE
     */
    int number_of_plates = plates_list_.size();
    int plate_chosen;
    
    for (size_t i = 0; i < bootstrap_.size(); i ++){
        for (size_t j = 0; j < number_of_plates; j ++){
            
            // pick plate
            plate_chosen = rand() % number_of_plates;
            
            // add to average
            bootstrap_[i] += (*correlation_plates_.find(plates_list_[plate_chosen])).second;
            
        }
        
        bootstrap_[i].Normalize();
        
    }
}

void CorrelationResults::CreateBinFiles(){
    /**
     EXPLANATION:
     Creates the bin files in where the pair information will be stored. If files already exist, resets them
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationResults
     
     FUNCITONS USED:
     ToStr
     */
    
    for (int bin = 0; bin < num_bins_; bin ++){
        
        // creating directory for the bin
        std::string directory_name = results_ + ToStr(bin)+"/";
        std::string command = "mkdir -p -v " + directory_name;
        system(command.c_str());
        
        // creating bin files for each of the plates
        for (PlatesMapSimple<CorrelationPlate>::map::iterator it = correlation_plates_.begin(); it != correlation_plates_.end(); it ++){
            
            std::string filename = directory_name + (*it).second.pairs_file_name();
            
            // writing headers in file (open the file erasing the previous content)
            std::ofstream bin_file(filename.c_str(),std::ofstream::trunc); 
            if (bin_file.is_open()){
                bin_file << "# object_RA object_DEC object_Z spectrum_RA spectrum_DEC pixel_Z pixel_dist pixel_number pixel_delta pixel_w pi sigma\n";
                bin_file.close();
            }
            else{
                std::cout << "Unable to open file:" << std::endl << filename << std::endl;
            }            
            
        }
        
    }
}

void CorrelationResults::NormalizeCrossCorrelation(){
    /**
     EXPLANATION:
     Normalizes the cross correlation results
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     CorrelationResults
     
     FUNCITONS USED:
     NONE
     */
    
    for (size_t i = 0; i < plates_list_.size(); i++){
        
        normalized_correlation_ += (*correlation_plates_.find(plates_list_[i])).second;
    
    }
    normalized_correlation_.Normalize();

    
}

void CorrelationResults::SaveCrossCorrelation(){
    /**
     EXPLANATION:
     Save the cross correlation results
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     CorrelationResults
     
     FUNCITONS USED:
     ToStr
     */
    
    std::string filename;
    
    // save plate contribution to cross-correlation in each of the bins
    for (int bin = 0; bin < num_bins_; bin ++){
        
        filename = output_base_name_ + ToStr(bin) + ".data";
        
        // open the file erasing the previous content)
        std::ofstream file(filename.c_str(),std::ofstream::trunc); 
        if (file.is_open()){
            file << "# Note: the following are not normalized" << std::endl;
            file << "# " << CorrelationPlate::InfoHeader() << std::endl;
            
            for (size_t i = 0; i < plates_list_.size(); i++){
                
                file << (*correlation_plates_.find(plates_list_[i])).second.Info(bin) << std::endl;
            }
            
            file.close();
        }
        else{
            std::cout << "Unable to open file:" << std::endl << filename << std::endl;
        }
    }
    
    // save normalized cross-correlation
    filename = output_base_name_ + ".data";
    {
        std::ofstream file(filename.c_str(),std::ofstream::trunc); 
        if (file.is_open()){
            for (size_t i = 0; i < num_bins_; i++){
                
                file << i << " " << normalized_correlation_.xi(i) << std::endl;
            }
            
            file.close();
        }
        else{
            std::cout << "Unable to open file:" << std::endl << filename << std::endl;
        }
    }
    filename = output_base_name_ + ".full.data";
    {
        std::ofstream file(filename.c_str(),std::ofstream::trunc); 
        if (file.is_open()){
            file << "# bin_index " << CorrelationPlate::InfoHeader() << std::endl;
        
            for (size_t i = 0; i < num_bins_; i++){
            
                file << i << " " << normalized_correlation_.Info(i) << std::endl;
            }
        
            file.close();
        }
        else{
            std::cout << "Unable to open file:" << std::endl << filename << std::endl;
        }
    }
    
    // save bootstrap realizations
    if (flag_compute_bootstrap_){
        for (size_t i = 0; i < bootstrap_.size(); i ++){
            filename = output_base_name_ + ".bootstrap" + ToStr(i) + ".data";
            {
                std::ofstream file(filename.c_str(),std::ofstream::trunc); 
                if (file.is_open()){
                    
                    for (size_t j = 0; j < num_bins_; j++){
                        
                        file << j << " " << bootstrap_[i].Info(j) << std::endl;
                    }
                    
                    file.close();
                }
                else{
                    std::cout << "Unable to open file:" << std::endl << filename << std::endl;
                }
            }
            
            filename = output_base_name_ + ".bootstrap" + ToStr(i) + ".full.data";
            {
                std::ofstream file(filename.c_str(),std::ofstream::trunc); 
                if (file.is_open()){
                    file << "# bin index " << CorrelationPlate::InfoHeader() << std::endl;
                    
                    for (size_t j = 0; j < num_bins_; j++){
                        
                        file << bootstrap_[i].Info(j) << " " << j << std::endl;
                    }
                    
                    file.close();
                }
                else{
                    std::cout << "Unable to open file:" << std::endl << filename << std::endl;
                }
            }
            
        }

    }

}


