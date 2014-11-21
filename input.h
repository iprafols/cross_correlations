/**
 input.h
 Purpose: This file defines the class Input. This class contains the constant variables which are used in a "global" sense
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 on 17/06/14
 
 */

#ifndef _Input_h
#define _Input_h

// libraries needed
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
////////

// classes needed
////////

// functions needed
#include "function_remove_blank_spaces.hpp"
////////

#include "typedefs.h"

class Input{
    
public:
    // -------------------------------------------------------------
    // constuctors
    
    // constructs object and initializes its variables
    Input(const std::string& filename = "");
    
    // -------------------------------------------------------------
    // access methods
    
    // access function for c_
    double c() const {return c_;}
        
    // access function for flag_compute_bootstrap_
    bool flag_compute_bootstrap() const {return flag_compute_bootstrap_;}
    
    // access function for flag_compute_plate_neighbours_
    bool flag_compute_plate_neighbours() const {return flag_compute_plate_neighbours_;}

    // access function for h
    double h() const {return h_;}
    
    // access function for H0_
    double h0() const {return h0_;}
    
    // access function for dataset2_
    std::string dataset2() const {return dataset2_;}
    
    // access function for dataset2_name_
    std::string dataset2_name() const {return dataset2_name_;}
    
    // access function for lya_spectra_dir_
    std::string lya_spectra_dir() const {return lya_spectra_dir_;}
    
    // access function for lya_wl_
    double lya_wl() const {return lya_wl_;}
    
    // access function for max_pi_
    double max_pi() const {return max_pi_;}
    
    // access function for max_sigma_
    double max_sigma() const {return max_sigma_;}
    
    // access function for neighbours_max_distance_
    double neighbours_max_distance() const {return neighbours_max_distance_;}
        
    // access function for num_bootstrap_
    size_t num_bootstrap() const {return num_bootstrap_;}
    
    // access function for num_bins_
    int num_bins() const {return num_bins_;}
    
    // access function for num_pi_bins_
    int num_pi_bins() const {return num_pi_bins_;}
    
    // access function for num_plates_
    int num_plates() const {return num_plates_;}
    
    // access function for num_points_interpolation_
    int num_points_interpolation() const {return num_points_interpolation_;}
    
    // access function for num_sigma_bins_
    int num_sigma_bins() const {return num_sigma_bins_;}
    
    // access function for object_catalog_
    std::string dataset1() const {return dataset1_;}
    
    // access function for object_catalog_name_
    std::string dataset1_name() const {return dataset1_name_;}
    
    // access function for output_
    std::string output() const {return output_;}
    
    // access function for output_base_name_
    std::string output_base_name() const {return output_base_name_;}
        
    // acces function for plate_neighbours_
    std::string plate_neighbours() const {return plate_neighbours_;}
    
    // access function for plots_
    std::string plots() const {return plots_;}
    
    // access function for input_
    std::string input() const {return input_;}
    
    // access function for results_
    std::string results() const {return results_;}
    
    // access function for step_pi_
    double step_pi() const {return step_pi_;}
    
    // access function for step_sigma_
    double step_sigma() const {return step_sigma_;}
    
    // access function for wm_
    double wm() const {return wm_;}
    
    // access function for z_max_
    double z_max() const {return z_max_;}
    
    // access function for z_max_interpolation_
    double z_max_interpolation() const {return z_max_interpolation_;}
    
    // access function for z_min_
    double z_min() const {return z_min_;}
    
    // access function for z_min_interpolation_
    double z_min_interpolation() const {return z_min_interpolation_;}
    
    
    
    
    
private:
    // -------------------------------------------------------------
    // flags
    
    // flag to determine whether to compute the bootstrap realizations or not
    bool flag_compute_bootstrap_;
    
    // flag to determine whether to compute the plate neighbours list or not
    bool flag_compute_plate_neighbours_;
    
    
    
    // -------------------------------------------------------------
    // input settings
    
    // input directory
    std::string input_;
    
    // number of plates 
    int num_plates_;
    
    // objects catalog filename
    std::string dataset1_;
    
    // objects catalog name
    std::string dataset1_name_;
    
    // name of the file containing the plate neighbours
    std::string plate_neighbours_;
    
    // spectra directory
    std::string lya_spectra_dir_;
    
    // spectra catalog filename
    std::string dataset2_;
    
    // spectra catalog name
    std::string dataset2_name_;
    
    
    
    // -------------------------------------------------------------
    // output settings
    
    // output directory
    std::string output_;
    
    // output base name
    std::string output_base_name_;
    
    // plots directory
    std::string plots_;
    
    // partial results directory
    std::string results_;
    
    
    
    // -------------------------------------------------------------
    // bin settings 
    
    // maximum value of parallel separation (in Mpc/h)
    double max_pi_;
    
    // maximum value of perpendicular separation (in Mpc/h)
    double max_sigma_;
    
    // maximum angular distance at which two plates are considered as neighbours (in radians)
    double neighbours_max_distance_;
    
    // number of bins in parallel separation
    int num_pi_bins_;
    
    // number of bins in perpendicular separation 
    int num_sigma_bins_;
    
    // Total number of bins
    int num_bins_;
    
    // step value of parallel separation (in Mpc/h)
    double step_pi_;
    
    // step value of perpendicular separation (in Mpc/h)
    double step_sigma_;

    
    
    
    
    
    
    
    // -------------------------------------------------------------
    // bootstrap settings
        
    // number of bootstrap realizations
    size_t num_bootstrap_;
    
    
    
    // -------------------------------------------------------------
    // Fidutial model
    
    // Hubble constatnt at present time (in km/s/Mpc)
    double h0_;
    
    // Reduced Hubble constant: H0/100
    double h_; 
    
    // Omega matter
    double wm_; 
    
    
    // -------------------------------------------------------------
    // line and redshift settings 
    
    // lyman-alpha wavelength (in Angstroms)
    double lya_wl_;
    
    // minimum redshift values for accepting a quasar
    double z_min_;
    
    // maximum redshift values for accepting a quasar
    double z_max_;
    
    // minimum redshift considered for the interpolation grid
    double z_min_interpolation_; 
    
    // maximum redshift considered for the interpolation grid
    double z_max_interpolation_; 
    
    // number of points the interpolation grid will have
    int num_points_interpolation_; 
    
    // -------------------------------------------------------------
    // Some mathematical and physical constants
    double c_; // speed of light (in km/s)
    
    // -------------------------------------------------------------
    // control variables

    // used parameters
    std::string used_params_;
    
    // unused parameters
    std::string unused_params_;
    
    
    // -------------------------------------------------------------
    // methods
    
    // Reads the input files and sets the specified parameters
    void ReadInputValues(const std::string& filename);
    
    // Set all parameters to default values
    void SetDefaultValues();
    
    // Set a specified variable to the specified value
    void SetValue(const std::string& name, const std::string& value, InputFlag& input_flag);
    
    // Update the composed parameters whenever necessary
    void UpdateComposedParams(const InputFlag& input_flag);
    
    // Write this run's configuration
    void WriteLog();
    
    // Write the used and the unused parameters
    void WriteParams();
    
};


#endif
