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
    if (flag_set_baofit_){
        command = "mkdir -p -v " + output_ + fit_;
        system(command.c_str());
    }
    if (flag_set_baofit_best_fit_){
        command = "mkdir -p -v " + output_ + best_fit_;
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
    flag_compute_plate_neighbours_ = false;
    flag_covariance_matrix_from_file_ = false;
    flag_load_only_ = false;
    flag_plot_ = true;
    flag_plot_catalog_info_ = flag_load_only_;
    flag_run_baofit_ = true;
    flag_run_baofit_best_fit_ = true;
    flag_set_baofit_ = true;
    flag_set_baofit_best_fit_ = true;
    flag_verbose_ = 1;
    flag_verbose_baofit_setup_ = flag_verbose_;
    flag_verbose_compute_plate_neighbours_ = flag_verbose_;
    flag_verbose_correlation_plate_ = flag_verbose_;
    flag_verbose_correlation_results_ = flag_verbose_;
    flag_verbose_covariance_matrix_ = flag_verbose_;
    flag_verbose_dla_dataset_ = flag_verbose_;
    flag_verbose_lya_spectra_dataset_ = flag_verbose_;
    flag_verbose_main_ = flag_verbose_;
    flag_verbose_pair_dataset_ = flag_verbose_;
    flag_verbose_plate_neighbours_ = flag_verbose_;
    flag_verbose_quasar_dataset_ = flag_verbose_;
    flag_write_partial_results_ = 0;
    
    
    // -------------------------------------------------------------
    // input settings
    input_ = "/triforce/catalogues/"; // Note: default value must always start with "/"
    dataset1_ = input_ + "DR11Q_alpha_v0.fits";
    dataset1_name_ = "DR11Q";
    dataset1_type_ = "quasar";
    dataset1_type_options_ = "quasar, dla";
    plate_neighbours_ = input_ + "plate_neighbours.dat";
    skip_plates_ = 0;
    lya_spectra_dir_ = input_ + "spectrum_fits_files/";
    dataset2_ = input_ + "DR11Q_spectra_forest_list.ls";
    dataset2_name_ = "DR11LyaF";
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
    // fit settings
    include_distorsions_ = false;
    baofit_model_root_ = "/Users/iprafols/Downloads/programes/baofit/models/";
    fit_ = "fit/";
    best_fit_ = "best_fit/";
    
    
    // -------------------------------------------------------------
    // bin setting
    neighbours_max_distance_ = 4.0*acos(-1.0)/180.0; // (in radians)
    max_pi_ = 50.0; // (in Mpc/h)
    max_sigma_ = 50.0; // (in Mpc/h)
    step_pi_ = 5.0; // (in Mpc/h)
    step_sigma_ = 5.0; // (in Mpc/h)
    num_pi_bins_ = int(2.0*max_pi_/step_pi_);
    num_sigma_bins_ = int(max_sigma_/step_sigma_);
    num_bins_ = num_pi_bins_*num_sigma_bins_;

    
    // -------------------------------------------------------------
    // bootstrap settings
    num_bootstrap_ = 10000;
    bootstrap_results_ = results_ + "bootstrap_realizations/";
    
    
    // -------------------------------------------------------------
    // lya autocorrelation settings
    lya_auto_correlation_ = input_ + "PalanqueDelabrouille_1DPk-Fft.out";
    lya_pixel_width_ = 210.0; // (in km/s)
    pixels_separation_ = 5; // (in number of pixels)
    
    
    // -------------------------------------------------------------
    // Fidutial model
    h0_ = 70;
    h_ = h0_/100.0;
    wm_ = 0.27;
      
    
    // -------------------------------------------------------------
    // line and redshift settings
    lya_wl_ = 1215.67;
    z_min_ = 2.0;
    z_max_ = 3.5;
    z_min_interpolation_ = 1.5; 
    z_max_interpolation_ = 4.0; 
    num_points_interpolation_ = 30000; 
    
    
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
    
    if (name == "baofit_model_root"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            baofit_model_root_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "bootstrap_results"){
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
    else if (name == "flag_covariance_matrix_from_file"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                flag_covariance_matrix_from_file_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                flag_covariance_matrix_from_file_ = false;
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
    else if (name == "flag_run_baofit"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                flag_run_baofit_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                flag_run_baofit_ = false;
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
    else if (name == "flag_run_baofit_best_fit"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                flag_run_baofit_best_fit_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                flag_run_baofit_best_fit_ = false;
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
    else if (name == "flag_set_baofit"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                flag_set_baofit_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                flag_set_baofit_ = false;
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
    else if (name == "flag_set_baofit_best_fit"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                flag_set_baofit_best_fit_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                flag_set_baofit_best_fit_ = false;
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
    else if (name == "flag_verbose_baofit_setup"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            flag_verbose_baofit_setup_ = atoi(value.c_str());
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
    else if (name == "include_distorsions"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            if (value == "true" or value == "TRUE" or value == "True"){
                include_distorsions_ = true;
            }
            else if (value == "false" or value == "FALSE" or value == "False"){
                include_distorsions_ = false;
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
    else if (name == "lya_auto_correlation"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            lya_auto_correlation_ = value;
            input_flag[name] = true;
        }
        else{
            std::cout << "Repeated line in input file: " << name << std::endl << "quiting..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else if (name == "lya_pixel_width"){
        InputFlag::iterator it = input_flag.find(name);
        if (it == input_flag.end()){
            lya_pixel_width_ = double(atof(value.c_str()));
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
    it2 = input_flag.find("results");  
    it3 = input_flag.find("bootstrap_results");
    if ((it != input_flag.end() or it2 != input_flag.end()) and it3 == input_flag.end()){
        if (it2 != input_flag.end()){
            bootstrap_results_ = results_ + "bootstrap_realizations/";
        }
        else{
            bootstrap_results_ = output_ + "partial_results/bootstrap_realizations/";
        }
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
    
    // updating flag_verbose_baofit_setup_ if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_baofit_setup");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_baofit_setup_ = flag_verbose_;
    }
    
    // updating flag_verbose_baofit_setup_ if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_compute_plate_neighbours");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_compute_plate_neighbours_ = flag_verbose_;
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
    
    // updating flag_verbose_dla_dataset_ if necessary
    it = input_flag.find("flag_verbose");
    it2 = input_flag.find("flag_verbose_dla_dataset");
    if (it != input_flag.end() and it2 == input_flag.end()){
        flag_verbose_dla_dataset_ = flag_verbose_;
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

    // updating lya_auto_correlation_ if necessary
    it = input_flag.find("input");
    it2 = input_flag.find("lya_auto_correlation");
    if (it != input_flag.end() and it2 == input_flag.end()){
        lya_auto_correlation_ = input_ + "PalanqueDelabrouille_1DPk-Fft.out";
    }
    else if (it2 != input_flag.end() and lya_auto_correlation_[0] != '/'){
        lya_auto_correlation_ = input_ + lya_auto_correlation_;
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

    // updating num_pi_bins_ and step_pi_ if necessary
    it = input_flag.find("max_pi");
    it2 = input_flag.find("step_pi");
    it3 = input_flag.find("num_pi_bins");
    if (it != input_flag.end() or it2 != input_flag.end() or it3 != input_flag.end()){
        if (it3 != input_flag.end()){
            step_pi_ = max_pi_/double(num_pi_bins_);
        }
        else{
            num_pi_bins_ = int(2.0*max_pi_/step_pi_);
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
    
    // updating output_base_name_ if necessary
    it = input_flag.find("dataset1_name");
    it2 = input_flag.find("dataset2_name");
    it3 = input_flag.find("output_dataset_name");
    if ((it != input_flag.end() or it2 != input_flag.end()) and it3 == input_flag.end()){
        output_base_name_ = dataset1_name_ + "-" + dataset2_name_;
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
    log.open((output_ + "log.param").c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
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
        if (flag_compute_plate_neighbours_){
            log << "flag_compute_plate_neighbours = true" << std::endl;
        }
        else{
            log << "flag_compute_plate_neighbours = false" << std::endl;
        }
        if (flag_covariance_matrix_from_file_){
            log << "flag_covariance_matrix_from_file = true" << std::endl;
        }
        else{
            log << "flag_covariance_matrix_from_file = false" << std::endl;
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
        if (flag_run_baofit_){
            log << "flag_run_baofit = true" << std::endl;
        }
        else{
            log << "flag_run_baofit = false" << std::endl;
        }
        if (flag_run_baofit_best_fit_){
            log << "flag_run_baofit_best_fit = true" << std::endl;
        }
        else{
            log << "flag_run_baofit_best_fit = false" << std::endl;
        }
        if (flag_set_baofit_){
            log << "flag_set_baofit = true" << std::endl;
        }
        else{
            log << "flag_set_baofit = false" << std::endl;
        }
        if (flag_set_baofit_best_fit_){
            log << "flag_set_baofit = true" << std::endl;
        }
        else{
            log << "flag_set_baofit = false" << std::endl;
        }
        log << "flag_verbose = " << flag_verbose_ << std::endl;
        log << "flag_verbose_baofit_setup = " << flag_verbose_baofit_setup_ << std::endl;
        log << "flag_verbose_compute_plate_neighbours = " << flag_verbose_compute_plate_neighbours_ << std::endl;
        log << "flag_verbose_correlation_plate = " << flag_verbose_correlation_plate_ << std::endl;
        log << "flag_verbose_correlation_results = " << flag_verbose_correlation_results_ << std::endl;
        log << "flag_verbose_covariance_matrix_ = " << flag_verbose_covariance_matrix_ << std::endl;
        log << "flag_verbose_dla_dataset = " << flag_verbose_dla_dataset_ << std::endl;
        log << "flag_verbose_lya_spectra_dataset = " << flag_verbose_lya_spectra_dataset_ << std::endl;
        log << "flag_verbose_main = " << flag_verbose_main_ << std::endl;
        log << "flag_verbose_pair_dataset = " << flag_verbose_pair_dataset_ << std::endl;
        log << "flag_verbose_plate_neighbours = " << flag_verbose_plate_neighbours_ << std::endl;
        log << "flag_verbose_quasar_dataset = " << flag_verbose_quasar_dataset_ << std::endl;
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
        log << "// fit settings" << std::endl;
        if (include_distorsions_){
            log << "include_distorsions = true" << std::endl;
        }
        else{
            log << "include_distorsions = false" << std::endl;
        }
        log << "baofit_model_root = " << baofit_model_root_ << std::endl;
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
        log << "// lya autocorrelation settings" << std::endl;
        log << "lya_auto_correlation = " << lya_auto_correlation_ << std::endl;
        log << "lya_pixel_width = " << lya_pixel_width_ << std::endl;
        log << "pixels_separation = " << pixels_separation_ << std::endl;
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
