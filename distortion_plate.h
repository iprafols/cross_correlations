/**
 distortion_plate.h
 Purpose: This file defines the class DistortionPlate. This class contains the variables necessary to store the measurement of the distortion matrix in a single plate
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 04/04/2016
 
 */

#ifndef _DistortionPlate_h
#define _DistortionPlate_h

// libraries needed
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
////////

// classes needed
#include "astro_object.h"
#include "astro_object_dataset.h"
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

class DistortionPlate: public Plate{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs empty object
    DistortionPlate(){};
    
    // constructs "bad data" DistortionPlate
    DistortionPlate(int bad_data);
    
    // constructs object and initializes its variables
    DistortionPlate(const Input& input, const int plate_number, const std::vector<int>& plate_neighbours);
    
    // constructs object and initializes its variables
    DistortionPlate(const int plate_number, const int num_bins, const std::vector<int>& plate_neighbours, size_t flag_verbose_distortion_plate);
    
    // -------------------------------------------------------------
    // access methods
    
    
    // access function for dist_mat_
    CovMat dist_mat() const {return dist_mat_;}
    
    // access function for flag_verbose_distortion_plate_
    size_t flag_verbose_distortion_plate() const {return flag_verbose_distortion_plate_;}
    
    // access function for num_bins_
    size_t num_bins() const {return num_bins_;}
    
    // access function for plate_neigbours_
    std::vector<int> plate_neighbours() const {return plate_neighbours_;}
    
    // access function for results_
    std::string results() const {return results_;}

    // access functions for weight_
    std::vector<double> weight() const {return weight_;}
    double weight(size_t index) const;
    
    // -------------------------------------------------------------
    // set methods

    // set dist_mat_
    void set_dist_mat(size_t i, size_t j, double value);
    
    // set flag_verbose_distortion_plate_
    void set_flag_verbose_distortion_plate(size_t value) {flag_verbose_distortion_plate_ = value;}
    
    // set num_averaged_pairs_
    void set_num_averaged_pairs(size_t index, int value);
    
    // set weight_
    void set_weight(size_t index, double value);
    
    // -------------------------------------------------------------
    // other methods
    
    // compute covariance matrix
    void ComputeDistMat(const AstroObjectDataset& object_list, const SpectraDataset& spectra_list, const Input& input);
    
    // Normalizes the cross correlation results
    void Normalize();    
    
    
    // -------------------------------------------------------------
    // operator overload
    
    // adition
    void operator+= (const DistortionPlate& other);
    
    // subtraction
    DistortionPlate operator- (const DistortionPlate& other);
    
    // multiplication
    DistortionPlate operator* (const DistortionPlate& other);
    
    
    
private:
    // distortion matrix
    CovMat dist_mat_;
    
    // verbose flag
    size_t flag_verbose_distortion_plate_;

    // number of bins
    size_t num_bins_;
    
    // plate number for neighbouring plates
    std::vector<int> plate_neighbours_;
    
    // results directory (missing the bin number)
    std::string results_;
    
    // weight
    std::vector<double> weight_;

    // -------------------------------------------------------------
    // other methods
    
    // adding contribution to distortion matrix in the specified bin
    void AddPair(const LyaPixel& pixel, const LyaPixel& pixel2, const size_t& i, const size_t& j, const double& forest_total_weight, const double& forest_mean_loglam, const double& forest_aux);
    
};


#endif
