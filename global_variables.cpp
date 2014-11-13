/**
 global_variables.cpp
 Purpose: This files contains the body for the functions defined in global_variables.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 06/17/2014
 */

#include "global_variables.h"

GlobalVariables::GlobalVariables(){
    /**
     EXPLANATION:
     Cosntructs a GlobalVariables instance and initializes all its variables
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     GlobalVariables
     
     FUNCITONS USED:
     NONE
     */
    
    //
    // general settings
    //
    pwd_ = "/Users/iprafols/cross_correlations/";
    //pwd_ = "/triforce/iprafols/cross_correlations/";
    results_ = pwd_ + "results2/";
    plots_ = pwd_ + "plots2/";
    objects_catalog_ = pwd_ + "DR11Q_alpha_v0.fits";
    objects_catalog_name_ = "DR11Q_alpha_v0";
    pairs_file_name_ = "qso_spectrum_pairs_plate_";
    correlation_file_name_ = results_ + "correlation_bin_";
    normalized_correlation_ = results_ + "normalized_correlation.dat";
    plate_neighbours_ = pwd_ + "plate_neighbours.dat";
    lya_spectra_dir_ = pwd_ + "spectrum_fits_files/";
    lya_spectra_catalog_ = pwd_ + "DR11Q_spectra_forest_one_spectrum.ls";// versió per fer proves
    //lya_spectra_catalog_ = pwd_ + "DR11Q_spectra_forest_some_spectrum.ls";// versió per fer proves
    //lya_spectra_catalog_ = pwd_ + "DR11Q_spectra_forest_list.ls"; // versió definitiva
    lya_spectra_catalog_name_ = "DR11Q_spectra_forest";
    num_plates_ = 2044; // DR11

    // bootstrap settings
    bootstrap_flag_ = true;
    num_bootstrap_ = 50000;
    bootstrap_ = results_ + "bootstrap_realization_";
    bootstrap_dispersion_squared_ = results_ + "bootstrap_dispersion_squared.dat";
    
    //
    // Fidutial model
    //
    h0_ = 68.0;
    h_ = h0_/100.0;
    wm_ = 0.3;
    
    //
    // bin setting
    //
    neighbours_max_distance_ = 4.0*acos(-1.0)/180.0; // (in radians)
    max_pi_ = 150.0; // (in Mpc/h)
    max_sigma_ = 150.0; // (in Mpc/h)
    step_pi_ = 5.0; // (in Mpc/h)
    step_sigma_ = 5.0; // (in Mpc/h)
    num_pi_bins_ = int(2.0*max_pi_/step_pi_);
    num_sigma_bins_ = int(max_sigma_/step_sigma_);
    num_bins_ = num_pi_bins_*num_sigma_bins_;
    
    //
    // line and redshift settings
    //
    lya_wl_ = 1215.67;
    z_min_ = 2.0;
    z_max_ = 3.5;
    z_min_interpolation_ = 1.5; 
    z_max_interpolation_ = 4.0; 
    num_points_interpolation_ = 30000; 
    
    //
    // Some mathematical and physical constants
    //
    c_ = 299792.458;
}