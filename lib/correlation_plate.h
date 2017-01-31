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
#include "lya_auto_interpolation_map.h"
#include "lya_mean_projected_deltas_interpolation_map.h"
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

class CorrelationPlate: public Plate{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs empty object
    CorrelationPlate(){};
    
    // constructs "bad data" CorrelationPlate
    CorrelationPlate(int bad_data);
    
    // constructs object and initializes its variables
    CorrelationPlate(const Input& input, const int plate_number, const std::vector<int>& plate_neighbours);
    
    // constructs object and initializes its variables
    CorrelationPlate(const int plate_number, const int num_bins, const std::string& results, const std::string& pairs_file_name, const std::vector<int>& plate_neighbours, size_t flag_verbose_correlation_plate, size_t flag_write_partial_results);
    
    // -------------------------------------------------------------
    // access methods
    
    // access function for flag_verbose_correlation_plate_
    size_t flag_verbose_correlation_plate() const {return flag_verbose_correlation_plate_;}
    
    // access function for flag_write_partial_results_
    size_t flag_write_partial_results() const {return flag_write_partial_results_;}
    
    // access function for flag_projection_correction_
    bool flag_projection_correction() const {return flag_projection_correction_;}
    
    // access functions for mean_pi_
    std::vector<double> mean_pi() const {return mean_pi_;}
    double mean_pi(size_t index) const;
    
    // access functions for mean_sigma_
    std::vector<double> mean_sigma() const {return mean_sigma_;}
    double mean_sigma(size_t index) const;
    
    // access function for mean_redshift_
    double mean_z() const {return mean_z_;}
    
    // access function for mean_z_in_bin_
    std::vector<double> mean_z_in_bin() const {return mean_z_in_bin_;}
    double mean_z_in_bin(size_t index) const;
    
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
    
    // access function for weight_z_
    double weight_z() const {return weight_z_;}
    
    // access functions for xi_
    std::vector<double> xi() const {return xi_;}
    double xi(size_t index) const;
    
    // access functions for xi_correction_
    std::vector<double> xi_correction() const {return xi_correction_;}
    double xi_correction(size_t index) const;

    
    // -------------------------------------------------------------
    // set methods
    
    // set flag_verbose_correlation_plate_
    void set_flag_verbose_correlation_plate(size_t value) {flag_verbose_correlation_plate_ = value;}
    
    // set mean_pi_
    void set_mean_pi(size_t index, double value);
    
    // set mean_sigma_
    void set_mean_sigma(size_t index, double value);
    
    // set mean_z_
    void set_mean_z(double mean_z) {mean_z_ = mean_z;}
    
    // set mean_z_in_bin_
    void set_mean_z_in_bin(size_t index, double value);
    
    // set num_averaged_pairs_
    void set_num_averaged_pairs(size_t index, int value);
    
    // set weight_
    void set_weight(size_t index, double value);
    
    // set weight_z_
    void set_weight_z(double weight_z) {weight_z_ = weight_z;}
    
    // set xi_
    void set_xi(size_t index, double value);
    
    // set xi_correction_
    void set_xi_correction(size_t index, double value);
    
    
    // -------------------------------------------------------------
    // other methods
    
    // compute cross-correlation
    void ComputeCrossCorrelation(const AstroObjectDataset& object_list, const SpectraDataset& spectra_list, const Input& input);
    
    //returns a string with the information for the selected bin
    std::string Info(size_t bin);
    
    //returns a string with the column information
    static std::string InfoHeader();
    
    // Normalizes the cross correlation results
    void Normalize();
    
    // Saves the cross correlation resutls
    void SaveCrossCorrelation(const Input& input);
    
    
    // -------------------------------------------------------------
    // operator overload
    
    // adition
    void operator+= (const CorrelationPlate& other);
    
    // subtraction
    CorrelationPlate operator- (const CorrelationPlate& other);
    
    // multiplication
    CorrelationPlate operator* (const CorrelationPlate& other);


private:
    // verbose flag
    size_t flag_verbose_correlation_plate_;

    // flag to write partial results
    size_t flag_write_partial_results_;
    
    // flag to compute the projection correction
    bool flag_projection_correction_;

    // maximum number of pairs stored in each bin
    size_t max_pairs_;
    
    // mean value of the projected deltas as a function of redshift
    LyaMeanProjectedDeltasInterpolationMap mean_proj_deltas_;
    
    // mean value of parallel separation in bin
    std::vector<double> mean_pi_;
    
    // mean value of perpendicular separation in bin
    std::vector<double> mean_sigma_;
    
    // mean redshift
    double mean_z_;
    
    // mean redshift in bin
    std::vector<double> mean_z_in_bin_;
    
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

    // mean redshift weight
    double weight_z_;
    
    // cross correlation in bin
    std::vector<double> xi_;
    
    // correction to the cross correlation in bin
    std::vector<double> xi_correction;

    
    // -------------------------------------------------------------
    // other methods
    
    // adding contribution to xi in the specified bin
    void AddPair(const int& k_index, const LyaPixel& pixel, const double& pi, const double& sigma);

    // keeps the pair information for latter storage
    void KeepPair(const int& k_index, const LyaSpectrum& lya_spectrum, const size_t& pixel_number, const size_t obj_plate, const size_t obj_num);
    
    // write down pair information in bin file
    void SavePairs(const int& k_index);
    
    
};


#endif
