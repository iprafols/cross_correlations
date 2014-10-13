/**
 correlation_results.h
 Purpose: This file defines the class CorrelationResults. This class contains the variables necessary to store the measurement of the cross-correlation 
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 06/17/2014

 */

#ifndef _CorrelationResults_h
#define _CorrelationResults_h


// libraries needed
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
////////

// classes needed
#include "astro_object_dataset.h"
#include "correlation_plate.h"
#include "global_variables.h"
#include "lya_spectra_dataset.h"
#include "plate_neighbours.h"
////////

// functions needed
#include "function_to_str.hpp"
////////

#include "typedefs.h"
#include "defines.h"

class CorrelationResults{

public:
    // -------------------------------------------------------------
    // constructors

    // constructs empty object
    CorrelationResults(){};

    // constructs object and initializes its variables
    CorrelationResults(const GlobalVariables& kGlobalVariables, const PlateNeighbours& kPlateNeighbours);

    // -------------------------------------------------------------
    // access methods
    
    // access function for correlation_file_name_
    std::string correlation_file_name() const {return correlation_file_name_;}
    
    // access function for correlation_plates_
    PlatesMapSimple<CorrelationPlate>::map correlation_plates() const {return correlation_plates_;}
    CorrelationPlate correlation_plates(int plate_num) const {return (*correlation_plates_.find(plate_num)).second;}
    
    // access function for normalized_correlation_
    CorrelationPlate normalized_correlation() const {return normalized_correlation_;}
    
    // access function for num_bins_
    size_t num_bins() const {return num_bins_;}
    
    // access function for num_pi_bins_
    size_t num_pi_bins() const {return num_pi_bins_;}
    
    // access function for num_sigma_bins_
    size_t num_sigma_bins() const {return num_sigma_bins_;}
    
    // access function for pairs_file_name_
    std::string pairs_file_name() const {return pairs_file_name_;}
    
    // access functions for plates_list_
    std::vector<int> plates_list() const {return plates_list_;}
    int plates_list(int index) const {return plates_list_[index];}
    
    // access function for results_
    std::string results() const {return results_;}
    
    // -------------------------------------------------------------
    // other methods
    
    // compute cross-correlation
    void ComputeCrossCorrelation(const AstroObjectDataset& object_list, const LyaSpectraDataset& spectra_list, const GlobalVariables& kGlobalVariables);
    
    // create bin files
    void CreateBinFiles();
    
    
    
    
    
    
private:
    // correlation file name
    std::string correlation_file_name_;
    
    // map containing the correlation in the different plates
    PlatesMapSimple<CorrelationPlate>::map correlation_plates_;
    
    // number of bins
    size_t num_bins_;
    
    // number of bins in parallel separation
    size_t num_pi_bins_;
    
    // number of bins in perpendicular separation
    size_t num_sigma_bins_;
    
    // pairs file name
    std::string pairs_file_name_;
    
    // normalized cross-correlation
    CorrelationPlate normalized_correlation_;
    
    // list of plates
    std::vector<int> plates_list_;
    
    // results directory
    std::string results_;
    
    
    // -------------------------------------------------------------
    // methods
    
    // Normalizes the cross correlation results
    void NormalizeCrossCorrelation();
    
    // save the cross correlation results
    void SaveCrossCorrelation();

};


#endif
