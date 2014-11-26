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
#include "input.h"
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
    CorrelationResults(const Input& input, const PlateNeighbours& kPlateNeighbours);

    // -------------------------------------------------------------
    // access methods
    
    // access function for bootstrap_
    std::vector<CorrelationPlate> bootstrap() const {return bootstrap_;}
    CorrelationPlate bootstrap(size_t i) const;
    
    // access function for output_base_name_
    std::string output_base_name() const {return output_base_name_;}
    
    // access function for correlation_plates_
    PlatesMapSimple<CorrelationPlate>::map correlation_plates() const {return correlation_plates_;}
    CorrelationPlate correlation_plates(int plate_num) const;
    
    // access function for compute_bootstrap_
    bool flag_compute_bootstrap() const {return flag_compute_bootstrap_;}

    // access function for normalized_correlation_
    CorrelationPlate normalized_correlation() const {return normalized_correlation_;}
    
    // access function for num_bins_
    size_t num_bins() const {return num_bins_;}
        
    // access function for pairs_file_name_
    std::string pairs_file_name() const {return pairs_file_name_;}
    
    // access functions for plates_list_
    std::vector<int> plates_list() const {return plates_list_;}
    int plates_list(int index) const;
    
    // access function for results_
    std::string results() const {return results_;}
    
    // -------------------------------------------------------------
    // other methods
        
    // compute cross-correlation
    void ComputeCrossCorrelation(const AstroObjectDataset& object_list, const LyaSpectraDataset& spectra_list, const Input& input);
    
    // create bin files
    void CreateBinFiles();
    
    
    
    
    
    
private:
    // bootstrap variables
    std::vector<CorrelationPlate> bootstrap_;
    
    // output base name
    std::string output_base_name_;
    
    // map containing the correlation in the different plates
    PlatesMapSimple<CorrelationPlate>::map correlation_plates_;
    
    // boolean to specify whether or not to compute the bootstrap realizations
    bool flag_compute_bootstrap_;

    // number of bins
    size_t num_bins_;
        
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
    
    // compute bootstrap realizations
    void ComputeBootstrapRealizations();

    // Normalizes the cross correlation results
    void NormalizeCrossCorrelation();
    
    // save the cross correlation results
    void SaveCrossCorrelation();

};


#endif
