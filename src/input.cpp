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
    
    // look for the path to the folder in whicht the program is being executed
    running_pwd_ = "";
    FILE* pipe = popen("pwd", "r");
    char buffer[128];
    while(!feof(pipe)) {
        if(fgets(buffer, 128, pipe) != NULL){
            running_pwd_ += buffer;
        }
    }
    pclose(pipe);
    std::cout << "Running from : " << running_pwd_;

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
    command = "mkdir -p -v " + output_;
    system(command.c_str());
    if (flag_write_partial_results_ >= 1){
        command = "mkdir -p -v " + results_;
        system(command.c_str());
    }
    command = "mkdir -p -v " + output_ + plots_;
    system(command.c_str());
    if (flag_compute_bootstrap_){
        command = "mkdir -p -v " + bootstrap_results_;
        system(command.c_str());
    }
    
    // write the used and unused parameters from the .ini file
    if (filename != ""){
        WriteParams();
    }
    
    // write down this run's configuration
    WriteLog();
    
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
        std::cout << "Warning: in Input::ReadInputValues(filename): The input filename shuold have .ini extension. Parameters set to default values, ignoring file..." << std::endl;
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
        std::cout << "Error : In Input::ReadInputValues : Could not read input file. Using default parameters instead" << std::endl;
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
    flag_compute_covariance_ = true;
    flag_compute_cross_correlation_ = true;
    flag_compute_distortion_ = true;
    flag_compute_plate_neighbours_ = false;
    flag_load_only_ = false;
    flag_plot_ = true;
    flag_plot_catalog_info_ = flag_load_only_;
    flag_project_deltas_ = true;
    flag_projection_correction_ = flag_project_deltas_
    flag_verbose_ = 1;
    flag_verbose_civ_spectra_dataset_ = flag_verbose_;
    flag_verbose_compute_plate_neighbours_ = flag_verbose_;
    flag_verbose_correlation_plate_ = flag_verbose_;
    flag_verbose_correlation_results_ = flag_verbose_;
    flag_verbose_covariance_matrix_ = flag_verbose_;
    flag_verbose_covariance_plate_ = flag_verbose_;
    flag_verbose_dla_dataset_ = flag_verbose_;
    flag_verbose_distortion_matrix_ = flag_verbose_;
    flag_verbose_distortion_plate_ = flag_verbose_;
    flag_verbose_lya_spectra_dataset_ = flag_verbose_;
    flag_verbose_main_ = flag_verbose_;
    flag_verbose_pair_dataset_ = flag_verbose_;
    flag_verbose_plate_neighbours_ = flag_verbose_;
    flag_verbose_quasar_dataset_ = flag_verbose_;
    flag_verbose_strong_lya_dataset_ = flag_verbose_;
    flag_write_partial_results_ = 0;
    
    
    // -------------------------------------------------------------
    // input settings
    input_ = "/triforce/catalogues/"; // Note: default value must always start with "/"
    dataset1_ = input_ + "DR12Q.fits";
    dataset1_name_ = "DR11Q";
    dataset1_type_ = "quasar";
    dataset1_type_options_ = "quasar, dla, strong_lya";
    plate_neighbours_ = input_ + "plate_neighbours.dat";
    skip_plates_ = 0;
    lya_spectra_dir_ = input_ + "spectrum_fits_files/";
    dataset2_ = input_ + "DR11Q_spectra_forest_list.ls";
    dataset2_name_ = "DR11LyaF";
    dataset2_type_ = "lya";
    dataset2_type_options_ = "lya, civ";
    num_plates_ = 2044; // DR11
    
    
    // -------------------------------------------------------------
    // output settings
    output_ = "output/";
    output_base_name_ = dataset1_name_ + "-" + dataset2_name_;
    results_ = output_ + "partial_results/";
    detailed_results_ = results_ + "detailed_info_bin_";
    //pairs_file_name_ = "detailed_info_plate_";
    plots_ = "plots/";
    
    
    // -------------------------------------------------------------
    // bin setting
    neighbours_max_distance_ = 4.0*acos(-1.0)/180.0; // (in radians)
    max_pi_ = 80.0; // (in Mpc/h)
    max_sigma_ = 80.0; // (in Mpc/h)
    step_pi_ = 2.0; // (in Mpc/h)
    step_sigma_ = 2.0; // (in Mpc/h)
    num_pi_bins_ = int(2.0*max_pi_/step_pi_);
    num_sigma_bins_ = int(max_sigma_/step_sigma_);
    num_bins_ = num_pi_bins_*num_sigma_bins_;

    
    // -------------------------------------------------------------
    // bootstrap settings
    num_bootstrap_ = 100;
    bootstrap_results_ = output_ + "bootstrap_realizations/";
    
    
    // -------------------------------------------------------------
    // lya autocorrelation and projection correction settings
    lya_projection_correction_ = output_ + "projection_correction.dat";
    if (flag_project_deltas_){
        lya_auto_correlation_1d_ = output_ + "lya_auto_correlation_1d_projected.dat";
        lya_auto_correlation_3d_ = output_ + "lya_auto_correlation_3d_projected.dat";
    }
    else{
        lya_auto_correlation_1d_ = output_ + "lya_auto_correlation_1d.dat";
        lya_auto_correlation_3d_ = output_ + "lya_auto_correlation_3d.dat";
    }
    pixels_separation_ = 5; // (in number of pixels)
    z_min_interpolation_ = 1.96;
    z_max_interpolation_ = 3.44;
    num_points_interpolation_ = 400;
    max_pi_auto_ = 50.0; // (in Mpc/h)
    max_sigma_auto_ = 50.0; // (in Mpc/h)
    step_pi_auto_ = 5.0; // (in Mpc/h)
    step_sigma_auto_ = 5.0; // (in Mpc/h)
    num_pi_bins_auto_ = int(max_pi_auto_/step_pi_auto_);
    num_sigma_bins_auto_ = int(max_sigma_auto_/step_sigma_auto_);
    num_bins_auto_ = num_pi_bins_auto_*num_sigma_bins_auto_;
    
    
    // -------------------------------------------------------------
    // Fidutial model
    h0_ = 67.74;
    h_ = h0_/100.0;
    wm_ = 0.3089;
      
    
    // -------------------------------------------------------------
    // line and redshift settings
    lya_wl_ = 1215.67; // in angs
    z_min_ = 2.0;
    z_max_ = 3.5;
    nhi_min_ = 20.0;
    nhi_max_ = 22.0;
    cnr_min_ = 3.0;
    rf_wl_min_ = 1026.0; // in angs
    rf_wl_max_ = 1195.395; // in angs
    rf_wl_forbidden_interval_.first = 1005.0; // in angs
    rf_wl_forbidden_interval_.second = 1037.0; // in angs
    lya_flux_min_ = -0.05;
    lya_flux_min_ = 0.25;
    
    // -------------------------------------------------------------
    // Some mathematical and physical constants
    c_ = 299792.458;
    
    
    // internal settings
    // -------------------------------------------------------------
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
    
    if (name == "bootstrap_results"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            bootstrap_results_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "c"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            c_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "cnr_min"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            cnr_min_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "dataset1"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            dataset1_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "dataset1_name"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            dataset1_name_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "dataset1_type"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            dataset1_type_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "dataset2"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            dataset2_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." <<        std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "dataset2_name"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            dataset2_name_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "dataset2_type"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            dataset2_type_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
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
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_compute_covariance"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                flag_compute_covariance_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                flag_compute_covariance_ = false;
            }
            else{
                unused_params_ += name + " = " + value + "\n";
                return;
            }
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_compute_cross_correlation"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                flag_compute_cross_correlation_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                flag_compute_cross_correlation_ = false;
            }
            else{
                unused_params_ += name + " = " + value + "\n";
                return;
            }
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_compute_distortion"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                flag_compute_distortion_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                flag_compute_distortion_ = false;
            }
            else{
                unused_params_ += name + " = " + value + "\n";
                return;
            }
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
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
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_load_only"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                flag_load_only_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                flag_load_only_ = false;
            }
            else{
                unused_params_ += name + " = " + value + "\n";
                return;
            }
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_plot"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                flag_plot_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                flag_plot_ = false;
            }
            else{
                unused_params_ += name + " = " + value + "\n";
                return;
            }
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_plot_catalog_info"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                flag_plot_catalog_info_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                flag_plot_catalog_info_ = false;
            }
            else{
                unused_params_ += name + " = " + value + "\n";
                return;
            }
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_project_deltas"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                flag_project_deltas_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                flag_project_deltas_ = false;
            }
            else{
                unused_params_ += name + " = " + value + "\n";
                return;
            }
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_projection_correction"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                flag_projection_correction_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                flag_projection_correction_ = false;
            }
            else{
                unused_params_ += name + " = " + value + "\n";
                return;
            }
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_write_partial_results"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_write_partial_results_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose_civ_spectra_dataset"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_civ_spectra_dataset_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose_compute_plate_neighbours"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_compute_plate_neighbours_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose_correlation_plate"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_correlation_plate_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose_correlation_results"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_correlation_results_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose_covariance_matrix"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_covariance_matrix_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." <<    std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose_covariance_plate"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_covariance_plate_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose_dla_dataset"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_dla_dataset_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose_distortion_matrix"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_distortion_matrix_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose_distortion_plate"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_distortion_plate_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose_lya_spectra_dataset"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_lya_spectra_dataset_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose_main"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_main_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose_pair_dataset"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_pair_dataset_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose_plate_neighbours"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_plate_neighbours_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose_quasar_dataset"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_quasar_dataset_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "flag_verbose_strong_lya_dataset"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_strong_lya_dataset_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
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
                std::exit(EXIT_FAILURE);
            }
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
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
                std::exit(EXIT_FAILURE);
            }
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "input"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            input_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "lya_spectra_dir"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            lya_spectra_dir_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "lya_wl"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            lya_wl_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "max_pi"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            max_pi_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "max_pi_auto"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            max_pi_auto_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "max_sigma"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            max_sigma_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "max_sigma_auto"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            max_sigma_auto_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "neighbours_max_distance"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            neighbours_max_distance_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "nhi_max"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            nhi_max_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "nhi_min"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            nhi_min_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "num_bootstrap"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            num_bootstrap_ = size_t(atoi(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "num_pi_bins"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            it = input_flag.find("step_pi");
            if (it == input_flag.end()){
                num_pi_bins_ = atoi(value.c_str());
                input_flag[name] = true;
            }
            else{
                std::cout << "Input file contains an entry for both num_pi_bins and step_pi, pick one" << std::endl << "quiting..." << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." <<    std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "num_pi_bins_auto"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            it = input_flag.find("step_pi_auto");
            if (it == input_flag.end()){
                num_pi_bins_auto_ = atoi(value.c_str());
                input_flag[name] = true;
            }
            else{
                std::cout << "Input file contains an entry for both num_pi_bins and step_pi, pick one" << std::endl << "quiting..." << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." <<    std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "num_plates"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            num_plates_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "num_points_interpolation"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            num_points_interpolation_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "num_sigma_bins"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            it = input_flag.find("step_sigma");
            if (it == input_flag.end()){
                num_sigma_bins_ = atoi(value.c_str());
                input_flag[name] = true;
            }
            else{
                std::cout << "Input file contains an entry for both num_sigma_bins and step_sigma, pick one" << std::endl << "quiting..." << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "num_sigma_bins_auto"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            it = input_flag.find("step_sigma_auto");
            if (it == input_flag.end()){
                num_sigma_bins_auto_ = atoi(value.c_str());
                input_flag[name] = true;
            }
            else{
                std::cout << "Input file contains an entry for both num_sigma_bins and step_sigma, pick one" << std::endl << "quiting..." << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "output"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            output_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "output_base_name"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            output_base_name_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "pixels_separation"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            pixels_separation_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "plate_neighbours"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            plate_neighbours_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "results"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            results_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "rf_wl_min"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            rf_wl_min_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "rf_wl_max"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            rf_wl_max_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "rf_wl_forbidden_interval"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            size_t open_parenthesis_position, close_parenthesis_position, comma_position;
            open_parenthesis_position = value.find('(');
            close_parenthesis_position = value.find(')');
            comma_position = value.find(',');
            
            rf_wl_forbidden_interval_.first = double(atof(value.substr(open_parenthesis_position + 1, comma_position).c_str()));
            rf_wl_forbidden_interval_.second = double(atof(value.substr(comma_position + 1, close_parenthesis_position).c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "skip_plates"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            skip_plates_ = atoi(value.c_str());
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "step_pi"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            it = input_flag.find("num_pi_bins");
            if (it == input_flag.end()){
                step_pi_ = double(atof(value.c_str()));
                input_flag[name] = true;
            }
            else{
                std::cout << "Input file contains an entry for both num_pi_bins and step_pi, pick one" << std::endl << "quiting..." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "step_pi_auto"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            it = input_flag.find("num_pi_bins_auto");
            if (it == input_flag.end()){
                step_pi_auto_ = double(atof(value.c_str()));
                input_flag[name] = true;
            }
            else{
                std::cout << "Input file contains an entry for both num_pi_bins and step_pi, pick one" << std::endl << "quiting..." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "step_sigma"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            it = input_flag.find("num_sigma_bins");
            if (it == input_flag.end()){
                step_sigma_ = double(atof(value.c_str()));
                input_flag[name] = true;
            }
            else{
                std::cout << "Input file contains an entry for both num_sigma_bins and step_sigma, pick one" << std::endl << "quiting..." << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "step_sigma_auto"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            it = input_flag.find("num_sigma_bins_auto");
            if (it == input_flag.end()){
                step_sigma_auto_ = double(atof(value.c_str()));
                input_flag[name] = true;
            }
            else{
                std::cout << "Input file contains an entry for both num_sigma_bins and step_sigma, pick one" << std::endl << "quiting..." << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "wm"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            wm_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "z_max"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            z_max_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "z_max_interpolation"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            z_max_interpolation_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "z_min"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            z_min_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "z_min_interpolation"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            z_min_interpolation_ = double(atof(value.c_str()));
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl; 
            std::exit(EXIT_FAILURE);
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
    
    InputFlag::const_iterator it, it2, it3, it4, it5, it6;
    
    // updating bootstrap_results_ if necessary
    it = input_flag.find("output");
    it2 = input_flag.find("bootstrap_results");
    if (it != input_flag.end() and it2 == input_flag.end()){
        bootstrap_results_ = output_ + "bootstrap_realizations/";
    }
    
    // updating dataset1_ if necessary
    it = input_flag.find("input");
    it2 = input_flag.find("dataset1");
    if (it != input_flag.end() and it2 == input_flag.end()){
        dataset1_ = input_ + "DR11Q_alpha_v0.fits";
    }
    else if (it2 != input_flag.end() and dataset1_[0] != '/'){
        dataset1_ = input_ + dataset1_;
    }
    
    // updating dataset2_ if necessary
    it = input_flag.find("input");
    it2 = input_flag.find("dataset2");
    if (it != input_flag.end() and it2 == input_flag.end()){
        dataset2_ = input_ + "DR11Q_spectra_forest_list.ls";
    }
    else if (it2 != input_flag.end() and dataset2_[0] != '/'){
        dataset2_ = input_ + dataset2_;
    }

    // updating detailed_results_ if necessary
    it = input_flag.find("output");
    it2 = input_flag.find("results");
    if (it != input_flag.end() or it2 != input_flag.end()){
        if (it2 != input_flag.end()){
            detailed_results_ = results_ + "detailed_info_bin_";
        }
        else{
            detailed_results_ =  output_ + "partial_results/detailed_info_bin_";
        }
    }

    // updating flag_plot_catalog_info_ if necessary
    it = input_flag.find("flag_load_only");
    it2 = input_flag.find("flag_plot_catalog_info");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_plot_catalog_info_ = flag_load_only_;
    }
    
    // updating flag_projection_correction if necessary
    it = input_flag.find("flag_project_deltas");
    it2 = input_flag.find("flag_compute_distortion");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_projection_correction_ = flag_project_deltas_;
    }
    flag_projection_correction_ = flag_projection_correction_ and flag_project_deltas_;
    
    // updating flag_verbose_compute_plate_neighbours_ if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_compute_plate_neighbours");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_compute_plate_neighbours_ = flag_verbose_;
    }
    
    // updating flag_verbose_civ_spectra_dataset_ if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_civ_spectra_dataset");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_civ_spectra_dataset_ = flag_verbose_;
    }

    // updating flag_verbose_correlation_plate_ if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_correlation_plate");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_correlation_plate_ = flag_verbose_;
    }
    
    // updating flag_verbose_correlation_results_ if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_correlation_results");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_correlation_results_ = flag_verbose_;
    }
    
    // updating flag_verbose_covariance_matrix_ if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_covariance_matrix");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_covariance_matrix_ = flag_verbose_;
    }
    
    // updating flag_verbose_covariance_plate if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_covariance_plate");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_covariance_plate_ = flag_verbose_;
    }
    
    // updating flag_verbose_dla_dataset_ if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_dla_dataset");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_dla_dataset_ = flag_verbose_;
    }
    
    // updating flag_verbose_distortion_matrix if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_distortion_matrix");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_distortion_matrix_ = flag_verbose_;
    }
    
    // updating flag_verbose_distortion_plate if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_distortion_plate");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_distortion_plate_ = flag_verbose_;
    }
    
    // updating flag_verbose_lya_spectra_dataset_ if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_lya_spectra_dataset");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_lya_spectra_dataset_ = flag_verbose_;
    }
    
    // updating flag_verbose_main_ if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_main");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_main_ = flag_verbose_;
    }
    
    // updating flag_verbose_pair_dataset_ if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_pair_dataset");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_pair_dataset_ = flag_verbose_;
    }

    // updating flag_verbose_plate_neighbours_ if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_plate_neighbours");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_plate_neighbours_ = flag_verbose_;
    }
    
    // updating flag_verbose_plate_neighbours_ if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_quasar_dataset");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_quasar_dataset_ = flag_verbose_;
    }
    
    // updating flag_verbose_strong_lya_dataset_ if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_strong_lya_dataset");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_strong_lya_dataset_ = flag_verbose_;
    }

    // updating lya_auto_correlation_ if necessary
    it = input_flag.find("output");
    it2 = input_flag.find("flag_project_deltas");
    it3 = input_flag.find("flag_compute_distortion");
    if (it != input_flag.end() or it2 != input_flag.end() or it3 != input_flag.end()){
        if (flag_project_deltas_){
            lya_auto_correlation_1d_ = output_ + "lya_auto_correlation_1d_projected.dat";
            lya_auto_correlation_3d_ = output_ + "lya_auto_correlation_3d_projected.dat";
        }
        else{
            lya_auto_correlation_1d_ = output_ + "lya_auto_correlation_1d.dat";
            lya_auto_correlation_3d_ = output_ + "lya_auto_correlation_3d.dat";
        }
    }
    
    // updating lya_projection_correction_ if necessary
    it = input_flag.find("output");
    if (it != input_flag.end()){
        lya_projection_correction_ = output_ + "projection_correction.dat";
    }
    
    // updating lya_spectra_dir_ if necessary
    it = input_flag.find("input");
    it2 = input_flag.find("lya_spectra_dir");
    if (it != input_flag.end() and it2 == input_flag.end()){
        lya_spectra_dir_ = input_ + "spectrum_fits_files/";
    }
    else if (it2 != input_flag.end() and lya_spectra_dir_[0] != '/'){
        lya_spectra_dir_ = input_ + lya_spectra_dir_;
    }

    // updating num_bins_ if necessary
    it = input_flag.find("max_pi");
    it2 = input_flag.find("step_pi");
    it3 = input_flag.find("num_pi_bins");
    it4 = input_flag.find("max_sigma");
    it5 = input_flag.find("step_sigma");
    it6 = input_flag.find("num_sigma_bins");
    if (it != input_flag.end() or it2 != input_flag.end() or it3 != input_flag.end() or it4 != input_flag.end() or it5 != input_flag.end() or it6 != input_flag.end()){
        if (it3 != input_flag.end() and it6 != input_flag.end()){
            num_bins_ = num_pi_bins_*num_sigma_bins_;
        }
        else if (it3 != input_flag.end()){
            num_bins_ = num_pi_bins_*int(max_sigma_/step_sigma_);
        }
        else if (it6 != input_flag.end()){
            num_bins_ = int(2.0*max_pi_/step_pi_)*num_sigma_bins_;
        }
        else{
            num_bins_ = int(2.0*max_pi_/step_pi_)*int(max_sigma_/step_sigma_);
        }
    }
    
    // updating num_bins_auto_ if necessary
    it = input_flag.find("max_pi_auto");
    it2 = input_flag.find("step_pi_auto");
    it3 = input_flag.find("num_pi_bins_auto");
    it4 = input_flag.find("max_sigma_auto");
    it5 = input_flag.find("step_sigma_auto");
    it6 = input_flag.find("num_sigma_bins_auto");
    if (it != input_flag.end() or it2 != input_flag.end() or it3 != input_flag.end() or it4 != input_flag.end() or it5 != input_flag.end() or it6 != input_flag.end()){
        if (it3 != input_flag.end() and it6 != input_flag.end()){
            num_bins_auto_ = num_pi_bins_auto_*num_sigma_bins_auto_;
        }
        else if (it3 != input_flag.end()){
            num_bins_auto_ = num_pi_bins_auto_*int(max_sigma_auto_/step_sigma_auto_);
        }
        else if (it6 != input_flag.end()){
            num_bins_auto_ = int(2.0*max_pi_auto_/step_pi_auto_)*num_sigma_bins_auto_;
        }
        else{
            num_bins_auto_ = int(2.0*max_pi_auto_/step_pi_auto_)*int(max_sigma_auto_/step_sigma_auto_);
        }
    }

    // updating num_pi_bins_ and step_pi_ if necessary
    it = input_flag.find("max_pi");
    it2 = input_flag.find("step_pi");
    it3 = input_flag.find("num_pi_bins");
    if (it != input_flag.end() or it2 != input_flag.end() or it3 != input_flag.end()){
        if (it3 != input_flag.end()){
            step_pi_ = 2.0*max_pi_/double(num_pi_bins_);
        }
        else{
            num_pi_bins_ = int(2.0*max_pi_/step_pi_);
        }
    }
    
    // updating num_pi_bins_auto and step_pi_auto if necessary
    it = input_flag.find("max_pi_auto");
    it2 = input_flag.find("step_pi_auto");
    it3 = input_flag.find("num_pi_bins_auto");
    if (it != input_flag.end() or it2 != input_flag.end() or it3 != input_flag.end()){
        if (it3 != input_flag.end()){
            step_pi_auto_ = max_pi_auto_/double(num_pi_bins_auto_);
        }
        else{
            num_pi_bins_auto_ = int(max_pi_auto_/step_pi_auto_);
        }
    }
    
    // updating num_sigma_bins_ and step_sigma_ if necessary
    it = input_flag.find("max_sigma");
    it2 = input_flag.find("step_sigma");
    it3 = input_flag.find("num_sigma_bins");
    if (it != input_flag.end() or it2 != input_flag.end() or it3 != input_flag.end()){
        if (it3 != input_flag.end()){
            step_sigma_ = max_sigma_/double(num_sigma_bins_);
        }
        else{
            num_sigma_bins_ = int(max_sigma_/step_sigma_);
        }
    }
    
    // updating num_sigma_bins_auto and step_sigma_auto if necessary
    it = input_flag.find("max_sigma_auto");
    it2 = input_flag.find("step_sigma_auto");
    it3 = input_flag.find("num_sigma_bins_auto");
    if (it != input_flag.end() or it2 != input_flag.end() or it3 != input_flag.end()){
        if (it3 != input_flag.end()){
            step_sigma_auto_ = max_sigma_auto_/double(num_sigma_bins_auto_);
        }
        else{
            num_sigma_bins_auto_ = int(max_sigma_auto_/step_sigma_auto_);
        }
    }
    
    // updating output_base_name_ if necessary
    it = input_flag.find("dataset1_name");
    it2 = input_flag.find("dataset2_name");
    it3 = input_flag.find("output_dataset_name");
    if ((it != input_flag.end() or it2 != input_flag.end()) and it3 == input_flag.end()){
        output_base_name_ = dataset1_name_ + "-" + dataset2_name_;
    }
    if (flag_project_deltas_){
        output_base_name_ += "_projected";
    }
    
    // updating plate_neighbours_ if necessary
    it = input_flag.find("input");
    it2 = input_flag.find("plate_neighbours");
    if (it != input_flag.end() and it2 == input_flag.end()){
        plate_neighbours_ = input_ + "plate_neighbours.dat";
    }
    else if (it2 != input_flag.end() and plate_neighbours_[0] != '/'){
        plate_neighbours_ = input_ + plate_neighbours_;
    }
    
    // updating results_ if necessary
    it = input_flag.find("output");
    it2 = input_flag.find("results");    
    if (it != input_flag.end() and it2 == input_flag.end()){
        results_ = output_ + "partial_results/";
    }

}

void Input::WriteLog(){
    /**
     EXPLANATION:
     Writes down this run's configuration
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Input
     
     FUNCITONS USED:
     NONE
     */
    std::ofstream log;
    if (flag_project_deltas_){
        log.open((output_ + "log_projected.param").c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
    }
    else{
        log.open((output_ + "log.param").c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
    }
    if (log.is_open()){
        
        log << std::endl;
        log << "// -------------------------------------------------------------" << std::endl;
        log << "// flags" << std::endl;
        if (flag_compute_bootstrap_){
            log << "flag_compute_bootstrap = true" << std::endl;
        }
        else{
            log << "flag_compute_bootstrap = false" << std::endl;
        }
        if (flag_compute_covariance_){
            log << "flag_compute_covariance = true" << std::endl;
        }
        else{
            log << "flag_compute_covariance = false" << std::endl;
        }
        if (flag_compute_cross_correlation_){
            log << "flag_compute_cross_correlation = true" << std::endl;
        }
        else{
            log << "flag_compute_cross_correlation = false" << std::endl;
        }
        if (flag_compute_distortion_){
            log << "flag_compute_distortion = true" << std::endl;
        }
        else{
            log << "flag_compute_distortion = false" << std::endl;
        }
        if (flag_compute_plate_neighbours_){
            log << "flag_compute_plate_neighbours = true" << std::endl;
        }
        else{
            log << "flag_compute_plate_neighbours = false" << std::endl;
        }
        if (flag_load_only_){
            log << "flag_load_only = true" << std::endl;
        }
        else{
            log << "flag_load_only = false" << std::endl;
        }
        if (flag_plot_){
            log << "flag_plot = true" << std::endl;
        }
        else{
            log << "flag_plot = false" << std::endl;
        }
        if (flag_plot_catalog_info_){
            log << "flag_plot_catalog_info = true" << std::endl;
        }
        else{
            log << "flag_plot_catalog_info = false" << std::endl;
        }
        if (flag_project_deltas_){
            log << "flag_project_deltas = true" << std::endl;
        }
        else{
            log << "flag_project_deltas = false" << std::endl;
        }
        if (flag_projection_correction_){
            log << "flag_projection_correction = true" << std::endl;
        }
        else{
            log << "flag_projection_correction = false" << std::endl;
        }
        log << "flag_verbose = " << flag_verbose_ << std::endl;
        log << "flag_verbose_civ_spectra_dataset = " << flag_verbose_civ_spectra_dataset_ << std::endl;
        log << "flag_verbose_compute_plate_neighbours = " << flag_verbose_compute_plate_neighbours_ << std::endl;
        log << "flag_verbose_correlation_plate = " << flag_verbose_correlation_plate_ << std::endl;
        log << "flag_verbose_correlation_results = " << flag_verbose_correlation_results_ << std::endl;
        log << "flag_verbose_covariance_matrix_ = " << flag_verbose_covariance_matrix_ << std::endl;
        log << "flag_verbose_covariance_plate = " << flag_verbose_covariance_plate_ << std::endl;
        log << "flag_verbose_dla_dataset = " << flag_verbose_dla_dataset_ << std::endl;
        log << "flag_verbose_distortion_matrix = " << flag_verbose_distortion_matrix_ << std::endl;
        log << "flag_verbose_distortion_plate_ = " << flag_verbose_distortion_plate_ << std::endl;
        log << "flag_verbose_lya_spectra_dataset = " << flag_verbose_lya_spectra_dataset_ << std::endl;
        log << "flag_verbose_main = " << flag_verbose_main_ << std::endl;
        log << "flag_verbose_pair_dataset = " << flag_verbose_pair_dataset_ << std::endl;
        log << "flag_verbose_plate_neighbours = " << flag_verbose_plate_neighbours_ << std::endl;
        log << "flag_verbose_quasar_dataset = " << flag_verbose_quasar_dataset_ << std::endl;
        log << "flag_verbose_strong_lya_dataset = " << flag_verbose_strong_lya_dataset_ << std::endl;
        log << "flag_write_partial_results = " << flag_write_partial_results_ << std::endl;
        log << std::endl;
        
        
        log << std::endl;
        log << "// -------------------------------------------------------------" << std::endl;
        log << "// input settings" << std::endl;
        log << "input = " << input_ << std::endl;
        log << "dataset1 = " << dataset1_ << std::endl;
        log << "dataset1_name = " << dataset1_name_ << std::endl;
        log << "dataset1_type = " << dataset1_type_ << std::endl;
        log << "plate_neighbours = " << plate_neighbours_ << std::endl;
        log << "skip_plates = " << skip_plates_ << std::endl;
        log << "lya_spectra_dir = " << lya_spectra_dir_ << std::endl;
        log << "dataset2 = " << dataset2_ << std::endl;
        log << "dataset2_name = " << dataset2_name_ << std::endl;
        log << "dataset2_type = " << dataset2_type_ << std::endl;
        log << "num_plates = " << num_plates_ << std::endl;
        log << std::endl;
        
        
        log << std::endl;
        log << "// -------------------------------------------------------------" << std::endl;
        log << "// output settings" << std::endl;
        log << "output = " << output_ << std::endl;
        log << "output_base_name = " << output_base_name_ << std::endl;
        log << "results = " << results_ << std::endl;
        log << std::endl;
        
        
        log << std::endl;
        log << "// -------------------------------------------------------------" << std::endl;
        log << "// bin setting" << std::endl;
        log << "neighbours_max_distance_ = " << neighbours_max_distance_ << " # (in radians)" << std::endl;
        log << "max_pi = " << max_pi_ << " # (in Mpc/h)" << std::endl;
        log << "max_sigma = " << max_sigma_ << " # (in Mpc/h)" << std::endl;
        log << "step_pi = " << step_pi_ << " # (in Mpc/h)" << std::endl;
        log << "step_sigma = " << step_sigma_ << " # (in Mpc/h)" << std::endl;
        log << std::endl;
        
        
        log << std::endl;
        log << "// -------------------------------------------------------------" << std::endl;
        log << "// bootstrap settings" << std::endl;
        log << "num_bootstrap = " << num_bootstrap_ << std::endl;
        log << "bootstrap_results = " << bootstrap_results_ << std::endl;
        log << std::endl;
        
        
        log << std::endl;
        log << "// -------------------------------------------------------------" << std::endl;
        log << "// lya autocorrelation and projection correction settings" << std::endl;
        log << "pixels_separation = " << pixels_separation_ << std::endl;
        log << "max_pi_auto = " << max_pi_auto_ << std::endl;
        log << "max_sigma_auto = " << max_sigma_auto_ << std::endl;
        log << "step_pi_auto = " << step_pi_auto_ << std::endl;
        log << "step_sigma_auto = " << step_sigma_auto_ << std::endl;
        log << std::endl;
        
        log << std::endl;
        log << "// -------------------------------------------------------------" << std::endl;
        log << "// Fidutial model" << std::endl;
        log << "h = " << h_ << std::endl;
        log << "wm = " << wm_ << std::endl;
        log << std::endl;
        
        
        log << std::endl;
        log << "// -------------------------------------------------------------" << std::endl;
        log << "// line and redshift settings" << std::endl;
        log << "lya_wl = " << lya_wl_ << std::endl;
        log << "z_min = " << z_min_ << std::endl;
        log << "z_max = " << z_max_ << std::endl;
        log << "z_min_interpolation = " << z_min_interpolation_ << std::endl;
        log << "z_max_interpolation = " << z_max_interpolation_ << std::endl;
        log << "num_points_interpolation = " << num_points_interpolation_ << std::endl;
        log << "nhi_min = " << nhi_min_ << std::endl;
        log << "nhi_max = " << nhi_max_ << std::endl;
        log << "cnr_min = " << cnr_min_ << std::endl;
        log << "rf_wl_min = " << rf_wl_min_ << std::endl;
        log << "rf_wl_max = " << rf_wl_max_ << std::endl;
        log << "rf_wl_forbidden_interval = (" << rf_wl_forbidden_interval_.first << ", " << rf_wl_forbidden_interval_.second << ")" << std::endl;
        log << std::endl;
        
        
        log << std::endl;
        log << "// -------------------------------------------------------------" << std::endl;
        log << "// Some mathematical and physical constants" << std::endl;
        log << "c = " << c_ << std::endl;
        log << std::endl;
        
        
        log.close();
    }
    else{
        std::cout << "Error : In Input::WriteLog : Unable to open file:" << std::endl << output_ << "log.param" << std::endl;
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
        std::cout << "Error : In Input::WriteParams : Unable to open file:" << std::endl << output_ << "used.param" << std::endl;
    }
    
    std::ofstream unused_params_file;
    unused_params_file.open((output_ + "unused.param").c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
    if (unused_params_file.is_open()){
        
        unused_params_file << unused_params_;
        unused_params_file.close();
    }
    else{
        std::cout << "Error : In Input::WriteParams : Unable to open file:" << std::endl << output_ << "unused.param" << std::endl;
    }

}
