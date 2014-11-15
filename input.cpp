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
    std::cout << "Checking whether or not there is a file to read" << std::endl;
    if (filename != ""){
        std::cout << "There is: " << filename << std::endl;

        std::cout << "Loading input parameters" << std::endl;
        ReadInputValues(filename);
    }
    else{
        std::cout << "There isn't" << std::
        endl;
    }
    
    // creating folders if required
    std::string command;
    command = "mkdir -p -v " + results_;
    system(command.c_str());
    command = "mkdir -p -v " + plots_;
    system(command.c_str());
    
    // write the used and unused parameters
    if (filename != ""){
        WriteParams();
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
    InputFlag input_flag; 
    
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
                SetValue(name, value, input_flag);
            }
            
        }
        input.close();
        
        UpdateComposedParams(input_flag);
        
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
    
    // -------------------------------------------------------------
    // flags
    flag_compute_bootstrap_ = true;
    flag_compute_plate_neighbours_ = false;
    
    
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
    //lya_spectra_catalog_ = pwd_ + "DR11Q_spectra_forest_one_spectrum.ls";// versió per fer proves
    //lya_spectra_catalog_ = pwd_ + "DR11Q_spectra_forest_some_spectrum.ls";// versió per fer proves
    lya_spectra_catalog_ = pwd_ + "DR11Q_spectra_forest_list.ls";
    lya_spectra_catalog_name_ = "DR11Q_spectra_forest";
    num_plates_ = 2044; // DR11
    
    // bootstrap settings
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

void Input::SetValue(const std::string& name, const std::string& value, InputFlag& input_flag){
    /**
     EXPLANATION:
     Sets a specified variable to the specified value 
     
     INPUTS:
     name - a string containing the name of the variable to set
     value - a string containing the value the variable has to contain
     input_flag - a InputFlag (a flag map) instance
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Input
     
     FUNCITONS USED:
     NONE
     */
    
    if (name == "bootstrap"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            bootstrap_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "bootstrap_dispersion_squared"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            bootstrap_dispersion_squared_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "c"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            c_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "correlation_file_name"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            correlation_file_name_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "flag_compute_bootstrap"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                flag_compute_bootstrap_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                flag_compute_bootstrap_ = false;
            }
            else{
                unused_params_ += name + " = " + value + "\n";
                return;
            }
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "flag_compute_plate_neighbours"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                flag_compute_plate_neighbours_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                flag_compute_plate_neighbours_ = false;
            }
            else{
                unused_params_ += name + " = " + value + "\n";
                return;
            }
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "H0"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            it = input_flag.find("h");
            if (it == input_flag.end()){
                h0_ = double(atof(value.c_str()));
                h_ = h0_/100.0;
                input_flag[name] = true;
            }
            else{
                std::cout << "Input file contains an entry for both H0 and h, pick one" << std::endl << "quiting..." << std::endl;
                std::exit;
            }
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit;
        }
    }
    else if (name == "h"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            it = input_flag.find("H0");
            if (it == input_flag.end()){
                h_ = double(atof(value.c_str()));
                h0_ = h_*100.0;
                input_flag[name] = true;
            }
            else{
                std::cout << "Input file contains an entry for both H0 and h, pick one" << std::endl << "quiting..." << std::endl;
                std::exit;
            }
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "lya_spectra_catalog"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            lya_spectra_catalog_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "lya_spectra_catalog_name"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            lya_spectra_catalog_name_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "lya_spectra_dir"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            lya_spectra_dir_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "lya_wl"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            lya_wl_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "max_pi"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            max_pi_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "max_sigma"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            max_sigma_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "neighbours_max_distance"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            neighbours_max_distance_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "normalized_correlation"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            normalized_correlation_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "num_bootstrap"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            num_bootstrap_ = size_t(atoi(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "num_bins"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            num_bins_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "num_pi_bins"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            num_pi_bins_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "num_plates"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            num_plates_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "num_points_interpolation"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            num_points_interpolation_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "num_sigma_bins"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            num_sigma_bins_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "objects_catalog"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            objects_catalog_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "objects_catalog_name"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            objects_catalog_name_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "output"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            output_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "pairs_file_name"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            pairs_file_name_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "plate_neighbours"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            plate_neighbours_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "plots"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            plots_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "pwd"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            pwd_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "results"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            results_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "step_pi"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            step_pi_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "step_sigma"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            step_sigma_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "wm"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            wm_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "z_max"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            z_max_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "z_max_interpolation"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            z_max_interpolation_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "z_min"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            z_min_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else if (name == "z_min_interpolation"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            z_min_interpolation_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::exit;
        }
    }
    else{
        unused_params_ += name + " = " + value + "\n";
        return;
    }
    used_params_ += name + " = " + value + "\n";
}

void Input::UpdateComposedParams(const InputFlag& input_flag){
    /**
     EXPLANATION:
     Updates the composed parameters whenever necessary
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Input
     
     FUNCITONS USED:
     NONE
     */
    
    InputFlag::const_iterator it, it2, it3, it4;
    
    // updating results_ if necessary
    it = input_flag.find("output");
    it2 = input_flag.find("results");    
    if (it != input_flag.end() and it2 == input_flag.end()){
        results_ = output_ + "results/";
    }
    
    //updating plots_ if necessary
    it = input_flag.find("output");
    it2 = input_flag.find("plots");
    if (it != input_flag.end() and it2 == input_flag.end()){
        plots_ = output_ + "plots/";
    }
    
    // updating objects_catalog_ if necessary
    it = input_flag.find("pwd");
    it2 = input_flag.find("objects_catalog");
    if (it != input_flag.end() and it2 == input_flag.end()){
        objects_catalog_ = pwd_ + "DR11Q_alpha_v0.fits";
    }
    
    // updating correlation_file_name_ if necessary
    it = input_flag.find("output");
    it2 = input_flag.find("results");
    it3 = input_flag.find("correlation_file_name");
    if ((it != input_flag.end() or it2 != input_flag.end()) and it3 == input_flag.end()){
        if (it2 == input_flag.end()){
            correlation_file_name_ = output_ + "results/correlation_bin_";
        }
        else{
            correlation_file_name_ = results_ + "correlation_bin_";
        }
    }
    
    // updating normalized_correlation_ if necessary
    it = input_flag.find("output");
    it2 = input_flag.find("results");
    it3 = input_flag.find("normalized_correlation");
    if ((it != input_flag.end() or it2 != input_flag.end()) and it3 == input_flag.end()){
        if (it2 == input_flag.end()){
            normalized_correlation_ = output_ + "results/normalized_correlation.dat";
        }
        else{
            normalized_correlation_ = results_ + "normalized_correlation.dat";
        }
    }
    
    // updating plate_neighbours_ if necessary
    it = input_flag.find("pwd");
    it2 = input_flag.find("plate_neighbours");
    if (it != input_flag.end() and it2 == input_flag.end()){
        plate_neighbours_ = pwd_ + "plate_neighbours.dat";
    }
    
    // updating lya_spectra_dir_ if necessary
    it = input_flag.find("pwd");
    it2 = input_flag.find("lya_spectra_dir");
    if (it != input_flag.end() and it2 == input_flag.end()){
        lya_spectra_dir_ = pwd_ + "spectrum_fits_files/";
    }
    
    // updating lya_spectra_catalog_ if necessary
    it = input_flag.find("pwd");
    it2 = input_flag.find("lya_spectra_catalog");
    if (it != input_flag.end() and it2 == input_flag.end()){
        lya_spectra_catalog_ = pwd_ + "DR11Q_spectra_forest_list.ls";
    }
    
    // updating normalized_correlation_ if necessary
    it = input_flag.find("output");
    it2 = input_flag.find("results");
    it3 = input_flag.find("bootstrap");
    if ((it != input_flag.end() or it2 != input_flag.end()) and it3 == input_flag.end()){
        if (it2 == input_flag.end()){
            bootstrap_ = output_ + "results/bootstrap_realization_";
        }
        else{
            bootstrap_ = results_ + "bootstrap_realization_";
        }
    }
    
    // updating bootstrap_dispersion_squared_ if necessary
    it = input_flag.find("output");
    it2 = input_flag.find("results");
    it3 = input_flag.find("bootstrap_dispersion_squared");
    if ((it != input_flag.end() or it2 != input_flag.end()) and it3 == input_flag.end()){
        if (it2 == input_flag.end()){
            bootstrap_dispersion_squared_ = output_ + "results/bootstrap_dispersion_squared.dat";
        }
        else{
            bootstrap_dispersion_squared_ = results_ + "bootstrap_dispersion_squared.dat";
        }
    }
    
    // updating num_pi_bins_ if necessary
    it = input_flag.find("max_pi");
    it2 = input_flag.find("step_pi");
    if (it != input_flag.end() or it2 != input_flag.end()){
        num_pi_bins_ = int(2.0*max_pi_/step_pi_);
    }
    
    // updating num_sigma_bins_ if necessary
    it = input_flag.find("max_sigma");
    it2 = input_flag.find("step_sigma");
    if (it != input_flag.end() or it2 != input_flag.end()){
        num_sigma_bins_ = int(max_sigma_/step_sigma_);
    }
    
    // updating num_bins_ if necessary
    it = input_flag.find("max_pi");
    it2 = input_flag.find("step_pi");
    it3 = input_flag.find("max_sigma");
    it4 = input_flag.find("step_sigma");
    if (it != input_flag.end() or it2 != input_flag.end() or it3 != input_flag.end() or it4 != input_flag.end()){
        num_bins_ = int(2.0*max_pi_/step_pi_)*int(max_sigma_/step_sigma_);
    }
    
}

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
