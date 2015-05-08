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
#include <stdio.h>
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
    
    // access function for best_fit_
    std::string best_fit() const {return best_fit_;}
    
    // access function for baofit_model_root_
    std::string baofit_model_root() const {return baofit_model_root_;}
    
    // access function for bootstrap_results_
    std::string bootstrap_results() const {return bootstrap_results_;}
    
    // access function for c_
    double c() const {return c_;}
            
    // access function for object_catalog_
    std::string dataset1() const {return dataset1_;}
    
    // access function for object_catalog_name_
    std::string dataset1_name() const {return dataset1_name_;}
    
    // access function for dataset1_type_
    std::string dataset1_type() const {return dataset1_type_;}
    
    // access function for dataset1_type_options_
    std::string dataset1_type_options() const {return dataset1_type_options_;}

    // access function for dataset2_
    std::string dataset2() const {return dataset2_;}
    
    // access function for dataset2_name_
    std::string dataset2_name() const {return dataset2_name_;}
    
    // access function for fit_
    std::string fit() const {return fit_;}
    
    // access function for detailed_results_
    std::string detailed_results() const {return detailed_results_;}
    
    // access function for flag_compute_bootstrap_
    bool flag_compute_bootstrap() const {return flag_compute_bootstrap_;}
    
    // access function for flag_compute_covariance_
    bool flag_compute_covariance() const {return flag_compute_covariance_;}
    
    // access function for flag_compute_cross_correlation_
    bool flag_compute_cross_correlation() const {return flag_compute_cross_correlation_;}

    // access function for flag_compute_plate_neighbours_
    bool flag_compute_plate_neighbours() const {return flag_compute_plate_neighbours_;}
        
    // access function for flag_load_only_
    bool flag_load_only() const {return flag_load_only_;}
    
    // access function for flag_covariance_matrix_from_file_
    bool flag_covariance_matrix_from_file() const {return flag_covariance_matrix_from_file_;}
    
    // access function for flag_plot_
    bool flag_plot() const {return flag_plot_;}
    
    // access function for flag_plot_catalog_info_
    bool flag_plot_catalog_info() const {return flag_plot_catalog_info_;}
    
    // access function for flag_run_baofit_
    bool flag_run_baofit() const {return flag_run_baofit_;}
    
    // access function for flag_run_baofit_best_fit_
    bool flag_run_baofit_best_fit() const {return flag_run_baofit_best_fit_;}
    
    // access function for flag_set_baofit_
    bool flag_set_baofit() const {return flag_set_baofit_;}
    
    // access function for flag_set_baofit_best_fit_
    bool flag_set_baofit_best_fit() const {return flag_set_baofit_best_fit_;}
    
    // access function for flag_write_partial_results_
    size_t flag_write_partial_results() const {return flag_write_partial_results_;}

    // access function for flag_verbose_
    size_t flag_verbose() const {return flag_verbose_;}
    
    // access function for flag_verbose_baofit_setup_
    size_t flag_verbose_baofit_setup() const {return flag_verbose_baofit_setup_;}
    
    // access function for flag_verbose_compute_plate_neighbours_
    size_t flag_verbose_compute_plate_neighbours() const {return flag_verbose_compute_plate_neighbours_;}
    
    // access function for flag_verbose_correlation_plate_
    size_t flag_verbose_correlation_plate() const {return flag_verbose_correlation_plate_;}
    
    // access function for flag_verbose_correlation_results_
    size_t flag_verbose_correlation_results() const {return flag_verbose_correlation_results_;}
    
    // access function for flag_verbose_covariance_matrix_
    size_t flag_verbose_covariance_matrix() const {return flag_verbose_covariance_matrix_;}
    
    // access function for flag_verbose_dla_dataset_
    size_t  flag_verbose_dla_dataset() const {return flag_verbose_dla_dataset_;}
    
    // access function for flag_verbose_lya_spectra_dataset_
    size_t flag_verbose_lya_spectra_dataset() const {return flag_verbose_lya_spectra_dataset_;}
    
    // access function for flag_verbose_main_
    size_t flag_verbose_main() const {return flag_verbose_main_;}
    
    // access function for flag_verbose_pair_dataset_
    size_t flag_verbose_pair_dataset() const {return flag_verbose_pair_dataset_;}
    
    // access function for flag_verbose_plate_neighbours_
    size_t flag_verbose_plate_neighbours() const {return flag_verbose_plate_neighbours_;}
    
    // access function for flag_verbose_quasar_dataset_
    size_t flag_verbose_quasar_dataset() const {return flag_verbose_quasar_dataset_;}
    
    // access function for h_
    double h() const {return h_;}    
    
    // access function for H0_
    double h0() const {return h0_;}
            
    // access function for include_distorsions_
    bool include_distorsions() const {return include_distorsions_;}
    
    // access function for lya_auto_correlation_
    std::string lya_auto_correlation() const {return lya_auto_correlation_;}
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
    
    // access function for num_bins_lya_auto_correlation_
    int num_bins_lya_auto_correlation() const {return num_bins_lya_auto_correlation_;}
    
    // access function for num_plates_
    int num_plates() const {return num_plates_;}
    
    // access function for num_points_interpolation_
    int num_points_interpolation() const {return num_points_interpolation_;}
    
    // access function for num_sigma_bins_
    int num_sigma_bins() const {return num_sigma_bins_;}
    
    // access function for output_
    std::string output() const {return output_;}
    
    // access function for output_base_name_
    std::string output_base_name() const {return output_base_name_;}
    
    // access function for pairs_file_name_
    std::string pairs_file_name() const {return pairs_file_name_;}
    
    // acces function for plate_neighbours_
    std::string plate_neighbours() const {return plate_neighbours_;}
    
    // access function for plots_
    std::string plots() const {return plots_;}
    
    // access function for input_
    std::string input() const {return input_;}
    
    // access function for results_
    std::string results() const {return results_;}
    
    // access function for running_pwd_
    std::string running_pwd() const {return running_pwd_;}
    
    // access function for skip_plates_
    int skip_plates() const {return skip_plates_;}
    
    // access function for step_lya_auto_correlation_
    double step_lya_auto_correlation() const {return step_lya_auto_correlation_;}
    
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
        
    // flag to compute the bootstrap realizations
    bool flag_compute_bootstrap_;
    
    // flag to compute the covariane matrix
    bool flag_compute_covariance_;
    
    // flag to compute the cross correlation
    bool flag_compute_cross_correlation_;
    
    // flag to compute the plate neighbours list
    bool flag_compute_plate_neighbours_;
    
    // flag to compute the covariance matrix from files containing the pairs information
    bool flag_covariance_matrix_from_file_;
    
    // flag to end program after loading the catalogs and plotting their information
    bool flag_load_only_;
    
    // flag to plot results
    bool flag_plot_;
    
    // flag to plot catalogues info
    bool flag_plot_catalog_info_;
    
    // flag to run baofit
    bool flag_run_baofit_;
    
    // flag to run baofit with with all the parameters fixed to the best-fit model
    bool flag_run_baofit_best_fit_;    
    
    // flag to set baofit ini file
    bool flag_set_baofit_;
    
    // flag to set baofit with with all the parameters fixed to the best-fit model
    bool flag_set_baofit_best_fit_;    
    
    // verbose flag
    size_t flag_verbose_;
        
    // baofit_setup verbose flag
    size_t flag_verbose_baofit_setup_;

    // compute_plate_neighbours verbose flag
    size_t flag_verbose_compute_plate_neighbours_;
    
    // correlation_plate verbose flag (only works if flag_verbose_correlation_results_ >= 1)
    size_t flag_verbose_correlation_plate_;
    
    // correlation_results verbose flag
    size_t flag_verbose_correlation_results_;
    
    // covariance_matrix verbose flag
    size_t flag_verbose_covariance_matrix_;

    // dla_dataset verbose flag
    size_t flag_verbose_dla_dataset_;
    
    // lya_spectra_dataset flag
    size_t flag_verbose_lya_spectra_dataset_;
    
    // main verbose flag
    size_t flag_verbose_main_;

    // pairs dataset verbose flag
    size_t flag_verbose_pair_dataset_;

    // plate_neighbours verbose flag
    size_t flag_verbose_plate_neighbours_;

    // quasar_dataset verbose flag
    size_t flag_verbose_quasar_dataset_;
    
    // flag to write partial results (ignored if flag_compute_covariance_ is set)
    size_t flag_write_partial_results_;


    
    // -------------------------------------------------------------
    // input settings
    
    // objects catalog filename
    std::string dataset1_;
    
    // objects catalog name
    std::string dataset1_name_;
    
    // objects type
    std::string dataset1_type_;
    
    // enabled objects type
    std::string dataset1_type_options_;
    
    // spectra catalog filename
    std::string dataset2_;
    
    // spectra catalog name
    std::string dataset2_name_;
    
    // input directory
    std::string input_;
    
    // name of the file containing the lyman-alpha auto-correlation
    std::string lya_auto_correlation_;
    
    // spectra directory
    std::string lya_spectra_dir_;
    
    // number of plates
    int num_plates_;
    
    // number of bins in the lyman-alpha auto-correlation
    int num_bins_lya_auto_correlation_;
    
    // name of the file containing the plate neighbours
    std::string plate_neighbours_;
    
    // number of plates that have to be skipped
    int skip_plates_;
    
    // step value of parallel separation in the lyman-alpha auto-correlation (in Mpc/h)
    double step_lya_auto_correlation_;

    
    
    
    // -------------------------------------------------------------
    // output settings
    
    // basename of the directory where the pairs' detailed information will be stored
    std::string detailed_results_;

    // output directory
    std::string output_;
    
    // output base name
    std::string output_base_name_;
    
    // plots directory
    std::string plots_;
    
    // partial results directory
    std::string results_;
    
    // basename of the files where the pairs' detailed information will be stored
    std::string pairs_file_name_;
    
    
    
    // -------------------------------------------------------------
    // fit settings
    
    // flag to include distorsions into the fitting model
    bool include_distorsions_;
    
    // baofit modelroot
    std::string baofit_model_root_;
    
    // directory where the fit results will be saved
    std::string fit_;
    
    // directory where the best-fit predictions will be saved
    std::string best_fit_;

    
    
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
    
    // bootstrap results folder
    std::string bootstrap_results_;
    
    
    
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
    
    // running directory
    std::string running_pwd_;
    
    
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
