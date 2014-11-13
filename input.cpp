/**
 input.cpp
 Purpose: This files contains the body for the functions defined in input.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 10/20/2014
 */

#include "input.h"

Input::Input(const std::string& filename){
    /**
     EXPLANATION:
     Constructs and initializes an Input instance
     
     INPUTS:
     filename[optional] - parameters file. Must end with ".ini". If it is not supplied, then default values for all parameters are found
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Input
     
     FUNCITONS USED:
     NONE
     */
    
    // load default values
    std::cout << "Loading default parameters" << std::endl;
    SetDefaultValues();
    
    // load input values
    std::cout << "checking whether or not there is a file to read" << std::endl;
    if (filename != ""){
        std::cout << "there is" << std::endl;
        std::cout << "Loading input parameters" << std::endl;
        ReadInputValues(filename);
    }
    else{
        std::cout << "there isn't" << std::
        endl;
    }

}

void Input::ReadInputValues(const std::string& filename){
    /**
     EXPLANATION:
     Reads the input files and sets the configured parameters
     
     INPUTS:
     filename - a string containing the name of the input file
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Input
     
     FUNCITONS USED:
     RemoveBlankSpaces
     */
    
    if (filename.substr(filename.size()-4,4) != ".ini"){
        std::cout << "Warning: the input filename shuold have .ini extension. Parameters set to default values, ignoring file..." << std::endl;
        return;
    }
    
    // declaring useful variables
    std::string line, name, value;
    size_t equal_sign_position, comment_sign_position;
    
    // open input filename
    std::ifstream input(filename.c_str());
    if (input.is_open()){
        std::string line("");        
        while (getline(input,line)){
            
            equal_sign_position = line.find('=');
            comment_sign_position = line.find('#');
            // if the line is a comment, ignore it
            if (comment_sign_position <= equal_sign_position){
                continue;
            }
            else{
                name = RemoveBlankSpaces(line.substr(0, equal_sign_position));
                value = RemoveBlankSpaces(line.substr(equal_sign_position+1, comment_sign_position-equal_sign_position+1));
                std::cout << "parameter read: " << name << " = " << value << std::endl;
                SetValue(name, value);
            }
            
        }
        WriteParams();
    }
    
    else{
        std::cout << "Error: could not read input file. Using default parameters instead" << std::endl;
    }
}
    
void Input::SetDefaultValues(){
    /**
     EXPLANATION:
     Sets all parameters to default values
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Input
     
     FUNCITONS USED:
     NONE
     */
    
    //
    // general settings
    //
    pwd_ = "../";
    output_ = "output/";
    results_ = output_ + "results/";
    plots_ = output_ + "plots/";
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
    
    //
    // internal settings
    //
    used_params_ = "";
    unused_params_ = "";
}

void Input::SetValue(const std::string& name, const std::string& value){
    /**
     EXPLANATION:
     Sets a specified variable to the specified value 
     
     INPUTS:
     name - a string containing the name of the variable to set
     value - a string containing the value the variable has to contain
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Input
     
     FUNCITONS USED:
     NONE
     */
    
    if (name == "bootstrap"){
        bootstrap_ = value;
    }
    else if (name == "bootstrap_dispersion_squared"){
        bootstrap_dispersion_squared_ = value;
    }
    else if (name == "bootstrap_flag"){
        if (value == "true" or value == "TRUE" or value == "True"){
            bootstrap_flag_ = true;
        }
        else{
            bootstrap_flag_ = false;
        }
    }
    else if (name == "c"){
        c_ = double(atof(value.c_str()));
    }
    else if (name == "correlation_file_name"){
        correlation_file_name_ = value;
    }
    else if (name == "H0"){
        h0_ = double(atof(value.c_str()));
        h_ = h0_/100.0;
    }
    else if (name == "lya_spectra_catalog"){
        lya_spectra_catalog_ = value;
    }
    else if (name == "lya_spectra_catalog_name"){
        lya_spectra_catalog_name_ = value;
    }
    else if (name == "lya_spectra_dir"){
        lya_spectra_dir_ = value;
    }
    else if (name == "lya_wl"){
        lya_wl_ = double(atof(value.c_str()));
    }
    else if (name == "max_pi"){
        max_pi_ = double(atof(value.c_str()));
    }
    else if (name == "max_sigma"){
        max_sigma_ = double(atof(value.c_str()));
    }
    else if (name == "neighbours_max_distance"){
        neighbours_max_distance_ = double(atof(value.c_str()));
    }
    else if (name == "normalized_correlation"){
        normalized_correlation_ = value;
    }
    else if (name == "num_bootstrap"){
        num_bootstrap_ = size_t(atoi(value.c_str()));
    }
    else if (name == "num_bins"){
        num_bins_ = atoi(value.c_str());
    }
    else if (name == "num_pi_bins"){
        num_pi_bins_ = atoi(value.c_str());
    }
    else if (name == "num_plates"){
        num_plates_ = atoi(value.c_str());
    }
    else if (name == "num_points_interpolation"){
        num_points_interpolation_ = atoi(value.c_str());
    }
    else if (name == "num_sigma_bins"){
        num_sigma_bins_ = atoi(value.c_str());
    }
    else if (name == "objects_catalog"){
        objects_catalog_ = value;
    }
    else if (name == "objects_catalog_name"){
        objects_catalog_name_ = value;
    }
    else if (name == "output"){
        output_ = value;
    }
    else if (name == "pairs_file_name"){
        pairs_file_name_ = value;
    }
    else if (name == "plate_neighbours"){
        plate_neighbours_ = value;
    }
    else if (name == "plots"){
        plots_ = value;
    }
    else if (name == "pwd"){
        pwd_ = value;
    }
    else if (name == "results"){
        results_ = value;
    }
    else if (name == "step_pi"){
        step_pi_ = double(atof(value.c_str()));
    }
    else if (name == "step_sigma"){
        step_sigma_ = double(atof(value.c_str()));
    }
    else if (name == "wm"){
        wm_ = double(atof(value.c_str()));
    }
    else if (name == "z_max"){
        z_max_ = double(atof(value.c_str()));
    }
    else if (name == "z_max_interpolation"){
        z_max_interpolation_ = double(atof(value.c_str()));
    }
    else if (name == "z_min"){
        z_min_ = double(atof(value.c_str()));
    }
    else if (name == "z_min_interpolation"){
        z_min_interpolation_ = double(atof(value.c_str()));
    }
    else{
        unused_params_ += name + " = " + value + "\n";
        return;
    }
    used_params_ += name + " = " + value + "\n";
};

void Input::WriteParams(){
    /**
     EXPLANATION:
     Writes the used and the unused parameters
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Input
     
     FUNCITONS USED:
     NONE
     */

    std::ofstream used_params_file;
    used_params_file.open((output_ + "used.param").c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
    if (used_params_file.is_open()){
        
        std::cout << std::endl << "used_params : " << used_params_ << std::endl;
        used_params_file << used_params_;
        used_params_file.close();
    }
    else{
        std::cout << "Unable to open file:" << std::endl << output_ << "used.param" << std::endl;
    }
    
    std::ofstream unused_params_file;
    unused_params_file.open((output_ + "unused.param").c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
    if (unused_params_file.is_open()){
        
        unused_params_file << unused_params_;
        unused_params_file.close();
    }
    else{
        std::cout << "Unable to open file:" << std::endl << output_ << "unused.param" << std::endl;
    }
}
