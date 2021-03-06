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
#include "omp.h"
////////

// classes needed
#include "astro_object_dataset.h"
#include "correlation_plate.h"
#include "input.h"
#include "plate_neighbours.h"
#include "spectra_dataset.h"
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
    
    // access function for bootstrap_results_
    std::string bootstrap_results() const {return bootstrap_results_;}
    
    // access function for detailed_results_
    std::string detailed_results() const {return detailed_results_;}
    
    // access function for flag_verbose_correlation_results_
    size_t flag_verbose_correlation_results() const {return flag_verbose_correlation_results_;}
    
    // access function for flag_write_partial_results_
    size_t flag_write_partial_results() const {return flag_write_partial_results_;}
    
    // access function for output_base_name_
    std::string output_base_name() const {return output_base_name_;}
    
    // access function for compute_bootstrap_
    bool flag_compute_bootstrap() const {return flag_compute_bootstrap_;}
    
    // access function for flag_compute_covariance_
    bool flag_compute_covariance() const {return flag_compute_covariance_;}

    // access function for normalized_correlation_
    CorrelationPlate normalized_correlation() const {return normalized_correlation_;}
    
    // access function for num_bins_
    size_t num_bins() const {return num_bins_;}
            
    // access functions for plates_list_
    std::vector<int> plates_list() const {return plates_list_;}
    int plates_list(int index) const;
    
    // access function for results_
    std::string results() const {return results_;}
    
    
    // -------------------------------------------------------------
    // other methods
        
    // compute cross-correlation
    void ComputeCrossCorrelation(const AstroObjectDataset& object_list, const SpectraDataset& spectra_list, const Input& input, const PlateNeighbours& kPlateNeighbours);
    
    // create bin files
    void CreateBinFiles();
    
    
    
    
    
    
private:
    // bootstrap variables
    std::vector<CorrelationPlate> bootstrap_;
    
    // bootstrap results folder
    std::string bootstrap_results_;
    
    // vector containing the correlation computed with different threads
    std::vector<CorrelationPlate> correlation_threads_;
    
    // detailed results directory (this is where the pairs detailed contribution to the corresponding plate and bin is stored)
    std::string detailed_results_;
    
    // boolean to specify whether or not to compute the bootstrap realizations
    bool flag_compute_bootstrap_;
    
    // boolean to specify whether or not to compute the covariance matrix
    bool flag_compute_covariance_;
    
    // boolean to specify whether or not the correction to the cross-correlation due to the projection of the delta field is computed
    bool flag_projection_correction_;

    // correlation_results verbose flag
    size_t flag_verbose_correlation_results_;
    
    // flag to write partial results
    size_t flag_write_partial_results_;
    
    // output base name
    std::string output_base_name_;
    
    
    
    // mean value of the projected deltas as a function of redshift
    LyaMeanProjectedDeltasInterpolationMap mean_proj_deltas_;
    
    // number of bins
    size_t num_bins_;
    
    // normalized cross-correlation
    CorrelationPlate normalized_correlation_;
    
    // list of plates
    std::vector<int> plates_list_;
    
    // results directory
    std::string results_;
    
    // number of plates that have to be skipped
    int skip_plates_;
    
    
    // -------------------------------------------------------------
    // methods
    
    // add the contribution of a single plate to the corresponding bootstrap realizations
    void AddToBootstrapRealizations(const std::vector<std::vector<size_t> >& picked_plates, const CorrelationPlate& plate);

    // normalize the cross correlation results
    void NormalizeCrossCorrelation();
    
    // save the cross correlation results
    void SaveCrossCorrelation();

};


#endif
