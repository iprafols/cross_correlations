/**
 correlation_results.cpp
 Purpose: This files contains the body for the functions defined in correlation_results.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 06/17/2014
 */

#include "correlation_results.h"


CorrelationResults::CorrelationResults(const GlobalVariables& kGlobalVariables, const PlateNeighbours& kPlateNeighbours){
    /**
     EXPLANATION:
     Cosntructs a CorrelationResults instance and initializes all its variables
     
     INPUTS:
     kGlobalVarialbes - object of type GlobalVariables
     kPlateNeighbours - a PlateNeighbours instance
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     CorrelationResults
     GlobalVariables
     
     FUNCITONS USED:
     NONE
     */
    
    // setting the number of bins and plates from kGlobalVariables
    num_bins_ = kGlobalVariables.num_bins();
    num_sigma_bins_ = kGlobalVariables.num_sigma_bins();
    num_pi_bins_ = kGlobalVariables.num_pi_bins();
    
    // setting the results directory and the pairs file name from kGlobalVariables
    results_ = kGlobalVariables.results() + "qso_spectrum_pairs_info_bin_";
    pairs_file_name_ = kGlobalVariables.pairs_file_name();
    correlation_file_name_ = kGlobalVariables.correlation_file_name();
    
    // initialization of the plates map
    plates_list_ = kPlateNeighbours.GetPlatesList();
    for (size_t i = 0; i < plates_list_.size(); i ++){
        correlation_plates_[plates_list_[i]] = CorrelationPlate(plates_list_[i], num_bins_, results_, pairs_file_name_, kPlateNeighbours.GetNeighboursList(plates_list_[i]));
    }
    
    // initialization of the normalized cross-correlation variable
    normalized_correlation_ = CorrelationPlate(_NORM_, num_bins_, results_, kGlobalVariables.normalized_correlation(), kPlateNeighbours.GetNeighboursList(_NORM_));
    
    // creating bin files
    CreateBinFiles();
}

void CorrelationResults::ComputeCrossCorrelation(const AstroObjectDataset& object_list, const LyaSpectraDataset& spectra_list, const GlobalVariables& kGlobalVariables){
    /**
     EXPLANATION:
     Computes the cross-correlation for all plates
     
     INPUTS:
     object_list - an AstroObjectDataset instance
     spectra_list - a LyaSpectraDataset instance
     kGlobalVariables - a GlobalVariables instance to load the bin settings
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObjectDataset
     CorrelationResults
     GlobalVariables
     LyaSpectraDataset
     
     FUNCITONS USED:
     NONE
     */
    
    for (size_t i = 0; i < plates_list_.size(); i++){
        std::cout << i << " out of " << plates_list_.size() << " computed" << std::endl;
        (*correlation_plates_.find(plates_list_[i])).second.ComputeCrossCorrelation(object_list, spectra_list, kGlobalVariables);
    }
    
    NormalizeCrossCorrelation();
    SaveCrossCorrelation();
    
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
        std::string command = "mkdir " + directory_name;
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
    
    // save plate contribution cross-correlation in each of the bins
    for (int bin = 0; bin < num_bins_; bin ++){
        
        filename = correlation_file_name_ + ToStr(bin) + ".dat";
        
        // open the file erasing the previous content)
        std::ofstream file(filename.c_str(),std::ofstream::trunc); 
        if (file.is_open()){
            file << "# Note: the following numbers are not normalized" << std::endl;
            file << CorrelationPlate::InfoHeader() << std::endl;
            
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
    filename = normalized_correlation_.pairs_file_name();
    
    // open the file erasing the previous content)
    std::ofstream file(filename.c_str(),std::ofstream::trunc); 
    if (file.is_open()){
        file << CorrelationPlate::InfoHeader() << " bin_index" << std::endl;
        
        for (size_t i = 0; i < num_bins_; i++){
            
            file << normalized_correlation_.Info(i) << " " << i << std::endl;
        }
        
        file.close();
    }
    else{
        std::cout << "Unable to open file:" << std::endl << filename << std::endl;
    }
    
    
    
}


