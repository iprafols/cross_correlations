/**
 baofit_setup.cpp
 Purpose: This files contains the body for the functions defined in baofit_setup.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 11/25/2014
 */

#include "baofit_setup.h"

BaofitSetup::BaofitSetup(const Input& input){
    /**
     EXPLANATION:
     Runs baofit
     
     INPUTS:
     input - object of type Input
     bootstrap - a boolean specifying if the bootstrap covariance matrix is to be used
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     BaofitSetup
     Input
     
     FUNCITONS USED:
     NONE
     */
    flag_verbose_baofit_setup_ = input.flag_verbose_baofit_setup();
}

void BaofitSetup::Run(const Input& input, const bool bootstrap){
    /**
     EXPLANATION:
     Runs baofit
     
     INPUTS:
     input - object of type Input
     bootstrap - a boolean specifying if the bootstrap covariance matrix is to be used
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     BaofitSetup
     Input
     
     FUNCITONS USED:
     NONE
     */
    std::string command;
    
    if (flag_verbose_baofit_setup_ >= 1){
        std::cout << "Running BAOFIT" << std::endl;
    }
    
    // changing to fit directory
    command = "cd " + input.output() + input.fit();
    if (flag_verbose_baofit_setup_ >= 1){
        std::cout << command << std::endl;
    }
    system(command.c_str());
    
    if (bootstrap){
        command = "baofit -i " + input.output_base_name() + "_baofit.bootstrap.diag.ini > " + input.output_base_name() + "_baofit.bootstrap.diag.log";
        if (flag_verbose_baofit_setup_ >= 1){
            std::cout << command << std::endl;
        }
        system(command.c_str());
    }
    else{
        command = "baofit -i " + input.output_base_name() + "_baofit.ini > " + input.output_base_name() + "_baofit.log";
        if (flag_verbose_baofit_setup_ >= 1){
            std::cout << command << std::endl;
        }
        system(command.c_str());
    }
    
    // changing back to call directory
    command = "cd " + input.running_pwd();
    if (flag_verbose_baofit_setup_ >= 1){
        std::cout << command;
    }
    system(command.c_str());
    
}

void BaofitSetup::RunBestFit(const Input& input, const bool bootstrap){
    /**
     EXPLANATION:
     Runs baofit with all the parameters fixed to the best-fit model
     
     INPUTS:
     input - object of type Input
     bootstrap - a boolean specifying if the bootstrap covariance matrix is to be used
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     BaofitSetup
     Input
     
     FUNCITONS USED:
     NONE
     */
    std::string command;
    
    if (flag_verbose_baofit_setup_ >= 1){
        std::cout << "Running BAOFIT" << std::endl;
    }
    
    // changing to fit directory
    command = "cd " + input.output() + input.best_fit();
    if (flag_verbose_baofit_setup_ >= 1){
        std::cout << command << std::endl;
    }
    system(command.c_str());
    
    if (bootstrap){
        command = "baofit -i " + input.output_base_name() + "_baofit.bootstrap.diag.ini > " + input.output_base_name() + "_baofit.bootstrap.diag.log";
        if (flag_verbose_baofit_setup_ >= 1){
            std::cout << command << std::endl;
        }
        system(command.c_str());
    }
    else{
        command = "baofit -i " + input.output_base_name() + "_baofit.ini > " + input.output_base_name() + "_baofit.log";
        if (flag_verbose_baofit_setup_ >= 1){
            std::cout << command << std::endl;
        }
        system(command.c_str());
    }
    
    // changing back to call directory
    command = "cd " + input.running_pwd();
    if (flag_verbose_baofit_setup_ >= 1){
        std::cout << command;
    }
    system(command.c_str());
    
}


void BaofitSetup::Set(const Input& input, const bool bootstrap){
    /**
     EXPLANATION:
     Sets baofit ini file
     
     INPUTS:
     input - object of type Input
     bootstrap - a boolean specifying if the bootstrap covariance matrix is to be used
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     BaofitSetup
     Input
     
     FUNCITONS USED:
     NONE
     */
    if (flag_verbose_baofit_setup_ >= 1){
        std::cout << "Setting BAOFIT ini file" << std::endl;
    }
    
    if (bootstrap){
        WriteIniFile(input, true);
        std::string command;
        command = "cp " + input.output() + input.output_base_name() + ".data " + input.output() + input.output_base_name() + ".bootstrap.diag.data";
        if (flag_verbose_baofit_setup_ >= 1){
            std::cout << command << std::endl;
        }
        system(command.c_str());
    }
    else{
        WriteIniFile(input);
    }
}

void BaofitSetup::SetBestFit(const Input& input, const bool bootstrap){
    /**
     EXPLANATION:
     Sets baofit ini file with all the parameters fixed to the best-fit model
     
     INPUTS:
     input - object of type Input
     bootstrap - a boolean specifying if the bootstrap covariance matrix is to be used
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     BaofitSetup
     Input
     
     FUNCITONS USED:
     NONE
     */
    if (flag_verbose_baofit_setup_ >= 1){
        std::cout << "Setting BAOFIT best-fit ini file" << std::endl;
    }
    
    if (bootstrap){
        WriteBestFitIniFile(input, true);
        std::string command;
        command = "cp " + input.output() + input.output_base_name() + ".data " + input.output() + input.output_base_name() + ".bootstrap.diag.data";
        if (flag_verbose_baofit_setup_ >= 1){
            std::cout << command << std::endl;
        }
        system(command.c_str());
    }
    else{
        WriteBestFitIniFile(input);
    }
}

void BaofitSetup::WriteBestFitIniFile(const Input& input, const bool bootstrap){
    /**
     EXPLANATION:
     Writes baofit ini file with all the parameters fixed to the best-fit model
     
     INPUTS:
     input - object of type Input
     bootstrap - a boolean specifying if the bootstrap covariance matrix is to be used
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     BaofitSetup
     Input
     
     FUNCITONS USED:
     NONE
     */
    
    std::string filename, fit_config_filename, line;
    size_t ini, end;
    
    if (bootstrap){
        filename = input.output() + input.best_fit() + input.output_base_name() + "_baofit.bootstrap.diag.ini";
        fit_config_filename = input.output() + input.fit() + input.output_base_name() + + "_baofit.bootstrap.diag_fit.config";
    }
    else{
        filename = input.output() + input.best_fit() + input.output_base_name() + "_baofit.ini";
        fit_config_filename = input.output() + input.fit() + input.output_base_name() + + "_baofit_fit.config";
    }
    
    std::ofstream file(filename.c_str(),std::ofstream::trunc);
    std::ifstream parameters(fit_config_filename.c_str());
    if (file.is_open() and parameters.is_open()){
        
        file << "#################################################################################" << std::endl;
        file << "## Best-Fit " << input.output_base_name() << " cross-correlation" << std::endl;
        file << "## See https://github.com/iprafols/cross_correlations.git to donwload the code used to compute the cross correlation" << std::endl;
        file << "#################################################################################" << std::endl;
        file << std::endl;
        file << "### Model Options ############################################" << std::endl;
        file << std::endl;
        file << "# Cosmology templates" << std::endl;
        if (input.include_distorsions()){
            file << "#kspace = true" << std::endl;
            file << "kspace-fft = true # recommended when including distorsions into the model" << std::endl;
        }
        else{
            file << "kspace = true" << std::endl;
            file << "#kspace-fft = true # recommended when including distorsions into the model" << std::endl;
        }
        file << "modelroot = " << input.baofit_model_root() << std::endl;
        file << "fiducial = DR9LyaMocksLCDM" << std::endl;
        file << "nowiggles = DR9LyaMocksLCDMSB" << std::endl;
        file << std::endl;
        file << "# Model configuration" << std::endl;
        file << "cross-correlation = yes" << std::endl;
        file << "anisotropic = yes" << std::endl;
        file << "decoupled = yes" << std::endl;
        file << std::endl;
        if (input.include_distorsions()){
            file << "# Broadband distortion model" << std::endl;
            file << "#dist-add = rP,rT=0:2,-3:1" << std::endl;
            file << std::endl;            
        }
        else{
            file << "# Broadband distortion model" << std::endl;
            file << "dist-add = rP,rT=0:2,-3:1" << std::endl;
            file << std::endl;
        }
        file << "# Parameter setup" << std::endl;
        while (getline(parameters,line)){
            ini = line.find_first_not_of("value");
            end = line.find(";");
            file << "model-config = fix" << line.substr(ini,end-ini) << std::endl;
        }
        file <<std::endl;
        file << "### Data Options #############################################" << std::endl;
        file << std::endl;
        if (input.flag_compute_bootstrap()){
            file << "data = " << input.output() << input.output_base_name() << ".bootstrap.diag" << std::endl;
        }
        else{
            file << "data = " << input.output() << input.output_base_name() << std::endl;
        }
        file << std::endl;
        file << "data-format = comoving-cartesian" << std::endl;
        file << std::endl;
        file << "axis1-bins = [-" << input.max_pi() << ":" << input.max_pi() << "]*" << input.num_pi_bins() << std::endl;
        file << "axis2-bins = [0:" << input.max_sigma() << "]*" << input.num_sigma_bins() << std::endl;
        file << "axis3-bins = {2.35717}" << std::endl;
        file << std::endl;
        file << "### Analysis Options #########################################" << std::endl;
        file << std::endl;
        file << "# Cuts to apply before fitting" << std::endl;                                                      
        file << "rmin = 1" << std::endl;
        file << "rmax = 220" << std::endl;
        file << std::endl;
        if (input.flag_compute_bootstrap()){
            file << "output-prefix = " << input.output_base_name() << ".bootstrap.diag_" << std::endl;   
        }
        else{
            file << "output-prefix = " << input.output_base_name() << "_" << std::endl;
        }
        file << std::endl;
        file << "# Generate a second set of outputs with the additive distortion turned off" << std::endl;
        file << "alt-config = fix[dist*]=0" << std::endl;
        file << std::endl;
        file << "# Do not dump multipoles (since the distortion model multipole integrals are singular)" << std::endl;                                                                                 
        file << "ndump = 0" << std::endl;
        
        file.close();
    }
    else{
        if (parameters.is_open()){
            std::cout << "Error : In BaofitSetup::WriteBestFitIniFile : Unable to open file:" << std::endl << filename << std::endl;
        }
        else{
            std::cout << "Error : In BaofitSetup::WriteBestFitIniFile : Unable to open file:" << std::endl << parameters << std::endl;
        }
    }
}

void BaofitSetup::WriteIniFile(const Input& input, const bool bootstrap){
    /**
     EXPLANATION:
     Writes baofit ini file
     
     INPUTS:
     input - object of type Input
     bootstrap - a boolean specifying if the bootstrap covariance matrix is to be used
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     BaofitSetup
     Input
     
     FUNCITONS USED:
     NONE
     */
    
    std::string filename;
    if (bootstrap){
        filename = input.output() + input.fit() + input.output_base_name() + "_baofit.bootstrap.diag.ini";
    }
    else{
        filename = input.output() + input.fit() + input.output_base_name() + "_baofit.ini";
    }
    
    std::ofstream file(filename.c_str(),std::ofstream::trunc); 
    if (file.is_open()){
            
        file << "#################################################################################" << std::endl;
        file << "## Fit " << input.output_base_name() << " cross-correlation" << std::endl;
        file << "## See https://github.com/iprafols/cross_correlations.git to donwload the code used to compute the cross correlation" << std::endl;
        file << "#################################################################################" << std::endl;
        file << std::endl;
        file << "### Model Options ############################################" << std::endl;
        file << std::endl;
        file << "# Cosmology templates" << std::endl;
        if (input.include_distorsions()){
            file << "#kspace = true" << std::endl;
            file << "kspace-fft = true # recommended when including distorsions into the model" << std::endl;
        }
        else{
            file << "kspace = true" << std::endl;
            file << "#kspace-fft = true # recommended when including distorsions into the model" << std::endl;
        }
        file << "modelroot = " << input.baofit_model_root() << std::endl;
        file << "fiducial = DR9LyaMocksLCDM" << std::endl;
        file << "nowiggles = DR9LyaMocksLCDMSB" << std::endl;
        file << std::endl;
        file << "# Model configuration" << std::endl;
        file << "cross-correlation = yes" << std::endl;
        file << "anisotropic = yes" << std::endl;
        file << "decoupled = yes" << std::endl;
        file << std::endl;
        file << "# Parameter setup" << std::endl;
        file << "model-config = value[beta]=1.1;" << std::endl;
        file << "model-config = fix[(1+beta)*bias]=-0.336;" << std::endl;
        file << "model-config = fix[gamma-bias]=0.9; fix[gamma-beta]=0;" << std::endl;
        file << "model-config = value[bias2]=3.64;" << std::endl;
        file << "model-config = fix[beta2*bias2]=0.962524;" << std::endl;
        file << "model-config = value[delta-v]=0;" << std::endl;
        file << "model-config = fix[BAO amplitude]=1;" << std::endl;
        file << "model-config = fix[BAO alpha-iso]; value[BAO alpha-p*]=1;" << std::endl;
        file << "model-config = fix[gamma-scale]=0;" << std::endl;
        file << std::endl;
        file << "## 2D chisq scan in BAO parameters" << std::endl;
        file << "#model-config = binning[BAO alpha-perp]={0.6:1.3}*50" << std::endl;
        file << "#model-config = binning[BAO alpha-parallel]={0.7:1.4}*50" << std::endl;
        file << std::endl;
        file << "## 2D chisq scan in linear bias & RSD parameters" << std::endl;
        file << "#model-config = binning[beta]={0.5:2.0}*40" << std::endl;
        file << "#model-config = binning[bias2]={2.3:3.8}*40" << std::endl;
        file << std::endl;
        if (input.include_distorsions()){
            file << "## Distorsion parameters" << std::endl;
            file << "model-config = value[cont-kc]=0.02" << std::endl;
            file << "model-config = value[cont-pc]=1" << std::endl;
            file << std::endl;
            file << "## Distorsion priors" << std::endl;
            file << "model-config = boxprior[cont-kc] @ (0,0.1)" << std::endl;
            file << "model-config = gaussprior[cont-kc] @ (0,0.02)" << std::endl;
            file << "model-config = boxprior[cont-pc] @ (0,3)" << std::endl;
            file << "model-config = gaussprior[cont-pc] @ (0,4)" << std::endl; 
            file << std::endl;
            file << "# Broadband distortion model" << std::endl;
            file << "#dist-add = rP,rT=0:2,-3:1" << std::endl;
            file << std::endl;

        }
        else{
            file << "## Distorsion parameters" << std::endl;
            file << "#model-config = value[cont-kc]=0.02" << std::endl;
            file << "#model-config = value[cont-pc]=1" << std::endl;
            file << std::endl;
            file << "## Distorsion priors" << std::endl;
            file << "#model-config = boxprior[cont-kc] @ (0,0.1)" << std::endl;
            file << "#model-config = gaussprior[cont-kc] @ (0,0.02)" << std::endl;
            file << "#model-config = boxprior[cont-pc] @ (0,3)" << std::endl;
            file << "#model-config = gaussprior[cont-pc] @ (0,4)" << std::endl;
            file << std::endl;
            file << "# Broadband distortion model" << std::endl;
            file << "dist-add = rP,rT=0:2,-3:1" << std::endl;
            file << std::endl;
        }
        file << "### Data Options #############################################" << std::endl;
        file << std::endl;
        if (input.flag_compute_bootstrap()){
            file << "data = " << input.output() << input.output_base_name() << ".bootstrap.diag" << std::endl;
        }
        else{
            file << "data = " << input.output() << input.output_base_name() << std::endl;
        }
        file << std::endl;
        file << "data-format = comoving-cartesian" << std::endl;
        file << std::endl;
        file << "axis1-bins = [-" << input.max_pi() << ":" << input.max_pi() << "]*" << input.num_pi_bins() << std::endl;
        file << "axis2-bins = [0:" << input.max_sigma() << "]*" << input.num_sigma_bins() << std::endl;
        file << "axis3-bins = {2.35717}" << std::endl;
        file << std::endl;
        file << "### Analysis Options #########################################" << std::endl;
        file << std::endl;
        file << "# Cuts to apply before fitting" << std::endl;                                                      
        file << "rmin = 40" << std::endl;
        file << "rmax = 180" << std::endl;
        file << std::endl;
        if (input.flag_compute_bootstrap()){
            file << "output-prefix = " << input.output_base_name() << ".bootstrap.diag_" << std::endl;   
        }
        else{
            file << "output-prefix = " << input.output_base_name() << "_" << std::endl;
        }
        file << std::endl;
        file << "# Generate a second set of outputs with the additive distortion turned off" << std::endl;
        file << "alt-config = fix[dist*]=0" << std::endl;
        file << std::endl;
        file << "# Do not dump multipoles (since the distortion model multipole integrals are singular)" << std::endl;                                                                                 
        file << "ndump = 0" << std::endl;

        file.close();
    }
    else{
        std::cout << "Error : In BaofitSetup::WriteIniFile : Unable to open file:" << std::endl << filename << std::endl;
    }
}