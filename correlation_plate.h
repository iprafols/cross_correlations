/**
 correlation_plate.h
 Purpose: This file defines the class CorrelationPlate. This class contains the variables necessary to store the measurement of the cross-correlation in a single plate
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 10/02/2014
 
 */

#ifndef _CorrelationPlate_h
#define _CorrelationPlate_h

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
#include "lya_spectra_dataset.h"
#include "lya_spectrum.h"
#include "pair.h"
#include "plate.h"
#include "sphere_point.h"
////////

// functions needed
#include "function_to_str.hpp"
////////

#include "defines.h"
#include "typedefs.h"

class CorrelationPlate: public Plate{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs empty object
    CorrelationPlate(){};
    
    // constructs "bad data" CorrelationPlate
    CorrelationPlate(int bad_data);
    
    // constructs object and initializes its variables
    CorrelationPlate(const Input& input, const int plate_number, const std::vector<int>& plate_neighbours, const bool flag_covariance);
    
    // constructs object and initializes its variables
    CorrelationPlate(const int plate_number, const int num_bins, const std::string& results, const std::string& pairs_file_name, const std::vector<int>& plate_neighbours, size_t flag_verbose_correlation_plate, size_t flag_write_partial_results, const bool flag_covariance);
    
    // -------------------------------------------------------------
    // access methods
    
    // access function for flag_compute_covariance_
    bool flag_covariance() const {return flag_covariance_;}
    
    // access function for cov_mat_
    CovMat cov_mat() const {return cov_mat_;}
    
    // access function for flag_verbose_correlation_plate_
    size_t flag_verbose_correlation_plate() const {return flag_verbose_correlation_plate_;}
    
    // access function for flag_write_partial_results_
    size_t flag_write_partial_results() const {return flag_write_partial_results_;}
    
    // access functions for mean_pi_
    std::vector<double> mean_pi() const {return mean_pi_;}
    double mean_pi(size_t index) const;
    
    // access functions for mean_sigma_
    std::vector<double> mean_sigma() const {return mean_sigma_;}
    double mean_sigma(size_t index) const;
    
    // access function for num_bins_
    size_t num_bins() const {return num_bins_;}
    
    // access funtion for num_averaged_pairs_
    std::vector<int> num_averaged_pairs() const {return num_averaged_pairs_;}
    int num_averaged_pairs(size_t index) const;
    
    // access function for pairs_file_name_
    std::string pairs_file_name() const {return pairs_file_name_;}
    
    // access function for plate_neigbours_
    std::vector<int> plate_neighbours() const {return plate_neighbours_;}
    
    // access function for results_
    std::string results() const {return results_;}

    // access functions for weight_
    std::vector<double> weight() const {return weight_;}
    double weight(size_t index) const;
    
    // access functions for xi_
    std::vector<double> xi() const {return xi_;}
    double xi(size_t index) const;
    
    
    // -------------------------------------------------------------
    // set methods

    // set cov_mat_
    void set_cov_mat(size_t i, size_t j, double value);
    
    // set flag_verbose_correlation_plate_
    void set_flag_verbose_correlation_plate(size_t value) {flag_verbose_correlation_plate_ = value;}
    
    // set mean_pi_
    void set_mean_pi(size_t index, double value);
    
    // set mean_sigma_
    void set_mean_sigma(size_t index, double value);
    
    // set num_averaged_pairs_
    void set_num_averaged_pairs(size_t index, int value);
    
    // set weight_
    void set_weight(size_t index, double value);
    
    // set xi_
    void set_xi(size_t index, double value);
    
    
    // -------------------------------------------------------------
    // other methods
    
    // compute covariance matrix
    void ComputeCovMat(const AstroObjectDataset& object_list, const LyaSpectraDataset& spectra_list, const Input& input);
    
    // compute cross-correlation
    void ComputeCrossCorrelation(const AstroObjectDataset& object_list, const LyaSpectraDataset& spectra_list, const Input& input);
    
    //returns a string with the information for the selected bin
    std::string Info(size_t bin);
    
    //returns a string with the column information
    static std::string InfoHeader();
    
    // Normalizes the cross correlation results
    void Normalize();    
    
    
    // -------------------------------------------------------------
    // operator overload
    
    // adition
    void operator+= (const CorrelationPlate& other);
    
    // subtraction
    CorrelationPlate operator- (const CorrelationPlate& other);
    
    // multiplication
    CorrelationPlate operator* (const CorrelationPlate& other);


    // -------------------------------------------------------------
    // static access methods
    
    // value of gamma/2
    static double half_gamma() {return CorrelationPlate::half_gamma_;}
    
    // value of (1+z_0)^(gamma/2)
    static double one_plus_z0_to_the_half_gamma() {return CorrelationPlate::one_plus_z0_to_the_half_gamma_;}

    
    
private:
    // covariance matrix
    CovMat cov_mat_;
    
    // boolean to specify whether or not to compute the covariance matrix
    bool flag_covariance_;
    
    // verbose flag
    size_t flag_verbose_correlation_plate_;

    // flag to write partial results
    size_t flag_write_partial_results_;

    // maximum number of pairs stored in each bin
    size_t max_pairs_;
    
    // mean value of parallel separation in bin
    std::vector<double> mean_pi_;
    
    // mean value of perpendicular separation in bin
    std::vector<double> mean_sigma_;
    
    // number of bins
    size_t num_bins_;
    
    // number of pairs averaged
    std::vector<int> num_averaged_pairs_;
    
    // pairs file name
    std::string pairs_file_name_;
    
    // pairs information
    std::vector<std::vector<Pair> > pairs_information_;
    
    // plate number for neighbouring plates
    std::vector<int> plate_neighbours_;
    
    // position inside pairs_information_ of the next pair information
    std::vector<size_t> position_;
    
    // results directory (missing the bin number)
    std::string results_;
    
    // weight
    std::vector<double> weight_;
    
    // cross correlation in bin
    std::vector<double> xi_;
    
    
    // -------------------------------------------------------------
    // other methods
    
    // adding contribution to xi in the specified bin
    void AddPair(const int& k_index, const LyaPixel& pixel, const double& pi, const double& sigma);

    // adding contribution to covariance matrix in the specified bin
    void AddPair(const LyaPixel& pixel1, const LyaPixel& pixel2, const size_t& i, const size_t& j);
    
    // keeps the pair information for latter storage
    void KeepPair(const int& k_index, const LyaSpectrum& lya_spectrum, const size_t& pixel_number);
    
    // write down pair information in bin file
    void SavePairs(const int& k_index);
    
    
    // -------------------------------------------------------------
    // static variables
    
    // value of gamma/2
    static double half_gamma_;
    
    // value of (1+z_0)^(gamma/2)
    static double one_plus_z0_to_the_half_gamma_;

    
};


#endif
