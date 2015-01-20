/**
 covariance_matrix.h
 Purpose: This file defines the class CovarianceMatrix. This class contains the variables necessary to compute and store the covariance matrix of the measurement.
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 11/17/2014
 
 */

#ifndef _CovarianceMatrix_h
#define _CovarianceMatrix_h

// libraries needed
#include <string>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>
////////

// classes needed
#include "correlation_plate.h"
#include "input.h"
#include "plate_neighbours.h"
#include "pair_dataset.h"
////////

// functions needed
////////

#include "typedefs.h"
#include "defines.h"

class CovarianceMatrix{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs object and initializes its variables
    CovarianceMatrix(const Input& input);
    
    
    // -------------------------------------------------------------
    // access methods
    
    // access function for bootstrap_cov_mat_
    CovMat bootstrap_cov_mat() const {return bootstrap_cov_mat_;}
    double bootstrap_cov_mat(size_t i, size_t j) const;
    
    // access function for covmat_
    CovMat cov_mat() const {return cov_mat_;}
    double cov_mat(size_t i, size_t j) const;
    
    // access function for flag_verbose_covariance_matrix_
    size_t flag_verbose_covariance_matrix() const {return flag_verbose_covariance_matrix_;}
    
    // access function for num_bins_
    size_t num_bins() const {return num_bins_;}

    // access function for output_base_name_
    std::string output_base_name() const {return output_base_name_;}

    
    
    // -------------------------------------------------------------
    // other methods
    
    // compute covariance matrix from bootstrap realizations
    void ComputeBootstrapCovMat(const std::vector<CorrelationPlate>& bootstrap);
    
    // compute covariance matrix
    void ComputeCovMat(const Input& input, const PlateNeighbours& kPlateNeighbours);

    
    
    
private:
    // bootstrap covariance matrix
    CovMat bootstrap_cov_mat_;
    
    // covariance matrix
    CovMat cov_mat_;
    
    // covariance_matrix verbose flag
    size_t flag_verbose_covariance_matrix_;
    
    // number of bins
    size_t num_bins_;
    
    // output base name
    std::string output_base_name_;
    
    // -------------------------------------------------------------
    // methods
    
    // save the bootstrap covariance matrix
    void SaveBootstrapCovMat();
    
};



#endif