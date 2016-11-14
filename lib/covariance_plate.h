/**
 correlation_plate.h
 Purpose: This file defines the class CovariancePlate. This class contains the variables necessary to store the measurement of the cross-correlation in a single plate
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 10/02/2014
 
 */

#ifndef _CovariancePlate_h
#define _CovariancePlate_h

// libraries needed
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
// TODO: remove this test
#include "omp.h"
// end of test
////////

// classes needed
#include "astro_object.h"
#include "astro_object_dataset.h"
#include "lya_auto_interpolation_map.h"
#include "lya_pixel.h"
#include "lya_spectrum.h"
#include "pair.h"
#include "plate.h"
#include "spectra_dataset.h"
#include "sphere_point.h"
////////

// functions needed
#include "function_to_str.hpp"
////////

#include "defines.h"
#include "typedefs.h"

class CovariancePlate: public Plate{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs empty object
    CovariancePlate(){};
    
    // constructs "bad data" CovariancePlate
    CovariancePlate(int bad_data);
    
    // constructs object and initializes its variables
    CovariancePlate(const Input& input, const int plate_number, const std::vector<int>& plate_neighbours);
    
    // constructs object and initializes its variables
    CovariancePlate(const int plate_number, const int num_bins, const std::vector<int>& plate_neighbours, size_t flag_verbose_covariance_plate);
    
    // -------------------------------------------------------------
    // access methods
    
    
    // access function for cov_mat_
    CovMat cov_mat() const {return cov_mat_;}
    
    // access functions for weight_
    //CovMat weight() const {return weight_;}
    std::vector<double> weight() const {return weight_;}
    double weight(size_t i) const;
    
    // access function for flag_verbose_covariance_plate_
    size_t flag_verbose_covariance_plate() const {return flag_verbose_covariance_plate_;}
    
    // access function for num_bins_
    size_t num_bins() const {return num_bins_;}
    
    // access function for plate_neigbours_
    std::vector<int> plate_neighbours() const {return plate_neighbours_;}

    
    // -------------------------------------------------------------
    // set methods

    // set cov_mat_
    void set_cov_mat(size_t i, size_t j, double value);
    
    // set flag_verbose_covariance_plate_
    void set_flag_verbose_covariance_plate(size_t value) {flag_verbose_covariance_plate_ = value;}
    
    // set num_averaged_pairs_
    void set_num_averaged_pairs(size_t index, int value);
    
    // set weight_
    //void set_weight(size_t i, size_t j, double value);
    void set_weight(size_t i, double value);
    
    // -------------------------------------------------------------
    // other methods
    
    // compute covariance matrix
    void ComputeCovMat(const AstroObjectDataset& object_list, const SpectraDataset& spectra_list, const Input& input, const std::vector<LyaAutoInterpolationMap>& lya_auto);
    
    // Normalizes the cross correlation results
    void Normalize();
    
    // Saves the covariance matrix in this plate
    void SaveCovMat(const Input& input);
    
    
    // -------------------------------------------------------------
    // operator overload
    
    // adition
    void operator+= (const CovariancePlate& other);
    
    // subtraction
    CovariancePlate operator- (const CovariancePlate& other);
    
    // multiplication
    CovariancePlate operator* (const CovariancePlate& other);


    // -------------------------------------------------------------
    // static access methods
    
    // value of gamma/2
    static double half_gamma() {return CovariancePlate::half_gamma_;}
    
    // value of (1+z_0)^(gamma/2)
    static double one_plus_z0_to_the_half_gamma() {return CovariancePlate::one_plus_z0_to_the_half_gamma_;}

    
    
private:
    // covariance matrix
    CovMat cov_mat_;
    
    // verbose flag
    size_t flag_verbose_covariance_plate_;

    // maximum number of pairs stored in each bin
    size_t max_pairs_;
    
    // number of bins
    size_t num_bins_;
    
    // pairs file name
    std::string pairs_file_name_;
    
    // plate number for neighbouring plates
    std::vector<int> plate_neighbours_;
    
    // weight
    //CovMat weight_;
    std::vector<double> weight_;

    // -------------------------------------------------------------
    // other methods
    
    // adding contribution to covariance matrix in the specified bin
    void AddPair(const LyaPixel& pixel1, const LyaPixel& pixel2, const size_t& i, const size_t& j, const LyaAutoInterpolationMap& lya_auto);
    
    // -------------------------------------------------------------
    // static variables
    
    // value of gamma/2
    static double half_gamma_;
    
    // value of (1+z_0)^(gamma/2)
    static double one_plus_z0_to_the_half_gamma_;

    
};


#endif
