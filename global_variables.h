/**
 global _variables.h
 Purpose: This file defines the class GlobalVariables. This class contains the constant variables which are used in a "global" sense

 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 on 17/06/14

 */

#ifndef _GlobalVariables_h
#define _GlobalVariables_h

// libraries needed
#include <cmath>
#include <string>
////////

// classes needed
////////

// functions needed
////////


class GlobalVariables{
    
public:
    // -------------------------------------------------------------
    // constuctors

    // constructs object and initializes its variables
    GlobalVariables();

    // -------------------------------------------------------------
    // access methods
    
    // access function for c_
    double c() const {return c_;}
    
    // access function for correlation_file_name_
    std::string correlation_file_name() const {return correlation_file_name_;}
    
    // access function for h
    double h() const {return h_;}
    
    // access function for H0_
    double h0() const {return h0_;}
    
    // access function for lya_spectra_catalog_
    std::string lya_spectra_catalog() const {return lya_spectra_catalog_;}
    
    // access function for lya_spectra_catalog_name_
    std::string lya_spectra_catalog_name() const {return lya_spectra_catalog_name_;}
    
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
    
    // access function for normalized_correlation_
    std::string normalized_correlation() const {return normalized_correlation_;}
    
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
    std::string objects_catalog() const {return objects_catalog_;}
    
    // access function for object_catalog_name_
    std::string objects_catalog_name() const {return objects_catalog_name_;}
    
    // access function for pairs_file_name_
    std::string pairs_file_name() const {return pairs_file_name_;}
    
    // acces function for plate_neighbours_
    std::string plate_neighbours() const {return plate_neighbours_;}
    
    // access function for plots_
    std::string plots() const {return plots_;}
    
    // access function for pwd_
    std::string pwd() const {return pwd_;}
    
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
    // general settings
    
    // base name of the files containing the cross-correlation measurements (filename = correlation_file_name_ + 'bin number' + ".dat"
    std::string correlation_file_name_;
    
    // name of the file where the normalized cross-correlation values will be saved
    std::string normalized_correlation_;

    // number of plates 
    int num_plates_;
    
    // objects catalog filename
    std::string objects_catalog_;
    
    // objects catalog name
    std::string objects_catalog_name_;
    
    // base name of the files in which pairs information is saved (filename = results_ + 'some_directory' + pairs_file_name_)
    std::string pairs_file_name_;
    
    // name of the file containing the plate neighbours
    std::string plate_neighbours_;

    // plots directory
    std::string plots_;
    
    // main directory
    std::string pwd_;
    
    // results directory
    std::string results_;
    
    // spectra directory
    std::string lya_spectra_dir_;
    
    // spectra catalog filename
    std::string lya_spectra_catalog_;
    
    // spectra catalog name
    std::string lya_spectra_catalog_name_;
    
    
    // -------------------------------------------------------------
    // Fidutial model
    
    // Hubble constatnt at present time (in km/s/Mpc)
    double h0_;
    
    // Reduced Hubble constant: H0/100
    double h_; 
    
    // Omega matter
    double wm_; 
    
    
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
    
};


#endif
