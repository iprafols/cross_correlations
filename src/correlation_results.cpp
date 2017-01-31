/**
 correlation_results.cpp
 Purpose: This files contains the body for the functions defined in correlation_results.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 06/17/2014
 */

#include "correlation_results.h"


CorrelationResults::CorrelationResults(const Input& input, const PlateNeighbours& kPlateNeighbours){
    /**
     EXPLANATION:
     Cosntructs a CorrelationResults instance and initializes all its variables
     
     INPUTS:
     input - object of type Input
     kPlateNeighbours - a PlateNeighbours instance
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     CorrelationResults
     Input
     
     FUNCITONS USED:
     NONE
     */
    
    // set flags from input
    flag_compute_bootstrap_ = input.flag_compute_bootstrap();
    flag_verbose_correlation_results_ = input.flag_verbose_correlation_results();
    flag_write_partial_results_ = input.flag_write_partial_results();
    flag_compute_covariance_ =  input.flag_compute_covariance();
    flag_projection_correction_ = input.flag_projection_correction();
    
    // setting the number of bins from input
    if (flag_verbose_correlation_results_ >= 3){
        std::cout << "CorrelationResults: setting the number of bins" << std::endl;
    }
    num_bins_ = input.num_bins();
    
    // setting the results directory and the pairs file name from input
    if (flag_verbose_correlation_results_ >= 3){
        std::cout << "CorrelationResults: setting the results directory" << std::endl;
    }
    results_ = input.results();
    detailed_results_ = input.detailed_results();
    output_base_name_ = input.output() + input.output_base_name();
    
    // initialization of the plates map
    if (flag_verbose_correlation_results_ >= 3){
        std::cout << "CorrelationResults: initialize plates map" << std::endl;
    }
    plates_list_ = kPlateNeighbours.GetPlatesList();
    int num_threads = atoi(std::getenv("OMP_NUM_THREADS"));
    correlation_threads_.reserve(num_threads);
    for (size_t i = 0; i < num_threads; i ++){
        correlation_threads_.push_back(CorrelationPlate(input, _NORM_, kPlateNeighbours.GetNeighboursList(_NORM_)));
    }
    skip_plates_ = input.skip_plates();

    
    // initialization of the normalized cross-correlation variable
    if (flag_verbose_correlation_results_ >= 3){
        std::cout << "CorrelationResults: initialize the normalized cross-correlation variable" << std::endl;
    }
    normalized_correlation_ = CorrelationPlate(input, _NORM_, kPlateNeighbours.GetNeighboursList(_NORM_));
    
    // initialization of the bootstrap variable
    if (flag_compute_bootstrap_){
        if (flag_verbose_correlation_results_ >= 3){
            std::cout << "CorrelationResults: initialize bootstrap varaibles" << std::endl;
        }
        bootstrap_results_ = input.bootstrap_results() + input.output_base_name();
        bootstrap_.reserve(input.num_bootstrap());
        for (size_t i = 0; i < input.num_bootstrap(); i++){
            bootstrap_.push_back(CorrelationPlate(input, _NORM_, kPlateNeighbours.GetNeighboursList(_NORM_)));
        }
    }
    
    // creating bin files
    if (flag_write_partial_results_ >= 2 and skip_plates_ == 0){
        if (flag_verbose_correlation_results_ >= 2){
            std::cout << "Creating detailed info files" << std::endl;
        }
        CreateBinFiles();
    }
}

CorrelationPlate CorrelationResults::bootstrap(size_t i) const {
    /**
     EXPLANATION:
     Access function for xi_
     
     INPUTS:
     i - index of the selected bootstrap_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationResults
     
     FUNCITONS USED:
     NONE
     */
    if (i < bootstrap_.size()){
        return bootstrap_[i];
    }
    else{
        CorrelationPlate cp(_BAD_DATA_INT_);
        return cp;
    }
}

int CorrelationResults::plates_list(int index) const {
    /**
     EXPLANATION:
     Access function for plates_list_
     
     INPUTS:
     index - index of the selected correlation_plates_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationResults
     
     FUNCITONS USED:
     NONE
     */
    if (index < plates_list_.size()){
        return plates_list_[index];
    }
    else{
        return _BAD_DATA_INT_;
    }
}

void CorrelationResults::ComputeCrossCorrelation(const AstroObjectDataset& object_list, const SpectraDataset& spectra_list, const Input& input, const PlateNeighbours& kPlateNeighbours){
    /**
     EXPLANATION:
     Computes the cross-correlation for all plates
     
     INPUTS:
     object_list - an AstroObjectDataset instance
     spectra_list - a SpectraDataset instance
     input - a Input instance to load the bin settings
     kPlateNeighbours - a PlateNeighbours instance
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObjectDataset
     CorrelationResults
     Input
     SpectraDataset
     
     FUNCITONS USED:
     NONE
     */
    
    if (input.flag_verbose() >= 1){
        std::cout << "Computing the cross-correlation" << std::endl;
    }
    
    // prepare bootstrap variables if necessary
    std::vector<std::vector<size_t> > picked_plates;
    if (flag_compute_bootstrap_){
        size_t number_of_plates = plates_list_.size();
        std::vector<size_t> aux;
        aux.resize(number_of_plates);
        picked_plates.resize(bootstrap_.size(), aux);
        for (size_t i = 0; i < bootstrap_.size(); i ++){
            
            for (size_t j = 0; j < number_of_plates; j ++){
                
                // pick plate
                picked_plates[i][j] = plates_list_[rand() % number_of_plates];
                
            }
        }
    }

    
    // loop over plates
    size_t plates_computed = 0;
    #pragma omp parallel for ordered schedule(dynamic)
    for (size_t i = skip_plates_; i < plates_list_.size(); i++){
        
        if (input.flag_verbose() >= 3){
            std::cout << "Loading plate " << plates_list_[i] << std::endl;
        }
        CorrelationPlate plate(input, plates_list_[i], kPlateNeighbours.GetNeighboursList(plates_list_[i]));
        
        #pragma omp critical (plates_computed)
        {
            plates_computed ++;
            if (flag_verbose_correlation_results_ >= 2 or (flag_verbose_correlation_results_ >= 1 and plates_computed == plates_computed/100*100)){
                #pragma omp critical (cout)
                {
                    std::cout << plates_computed << " out of " << plates_list_.size() << " plates computed" << std::endl;
                }
            }
            else{
                plate.set_flag_verbose_correlation_plate(0);
            }
        }
        
        // compute cross-correlation in selected plate
        plate.ComputeCrossCorrelation(object_list, spectra_list, input);
        
        // save the cross-correlation in selected plate
        if (flag_write_partial_results_ >= 1){
            if (input.flag_verbose() >= 1){
                std::cout << "saving the cross-correlation" << std::endl;
            }
            plate.SaveCrossCorrelation(input);
        }
        
        // add to total value
        int thread_num = omp_get_thread_num();
        correlation_threads_[thread_num] += plate;
        
        // add to bootstrap realizations
        if (flag_compute_bootstrap_){
            AddToBootstrapRealizations(picked_plates, plate);
        }
                                       
    }
    
    // normalize cross-correlation
    NormalizeCrossCorrelation();
    
    // save cross-correlation measurements
    SaveCrossCorrelation();
    
}

void CorrelationResults::AddToBootstrapRealizations(const std::vector<std::vector<size_t> >& picked_plates, const CorrelationPlate& plate){
    /**
     EXPLANATION:
     Adds the contribution of a single plate to the corresponding bootstrap realizations
     
     INPUTS:
     picked_plates - a vector containing the selected plates for each of the bootstrap realizations
     plate - a CorrelationPlate with the correlation measured in a single plate
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     CorrelationResults
     
     FUNCITONS USED:
     NONE
     */
    
    if (flag_verbose_correlation_results_ >= 2){
        std::cout << "Adding plate to bootstrap realizations" << std::endl;
    }
    
    #pragma omp parallel for ordered schedule(dynamic)
    for (size_t i = 0; i < bootstrap_.size(); i ++){
        
        for (size_t j = 0; j < picked_plates[i].size(); j ++){
            
            // add plate to bootstrap realization
            if (picked_plates[i][j] == plate.plate_number()){
                bootstrap_[i] += plate;
            }
        }
    }
    
    
}

void CorrelationResults::CreateBinFiles(){
    /**
     EXPLANATION:
     Creates the bin files in where the pair information will be stored. If files already exist, resets them
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationResults
     
     FUNCITONS USED:
     ToStr
     */
    
    for (int bin = 0; bin < num_bins_; bin ++){
        
        std::string filename = detailed_results_ + ToStr(bin) + ".dat";
        if (flag_verbose_correlation_results_ >= 1){
            std::cout << "creating file:" << std::endl << filename << std::endl;
        }
        
        // writing headers in file (open the file erasing the previous content)
        std::ofstream bin_file(filename.c_str(),std::ofstream::trunc);
        if (bin_file.is_open()){
            bin_file << "# obj_plate obj_num spec_plate spec_fiber spec_MJD pixel_number pixel_delta pixel_z pixel_w \n";
            bin_file.close();
        }
        else{
            std::cout << "Error : In CorrelationResults::CreateBinFiles : Unable to open file:" << std::endl << filename << std::endl;
        }
    }
}

void CorrelationResults::NormalizeCrossCorrelation(){
    /**
     EXPLANATION:
     Normalizes the cross correlation results
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     CorrelationResults
     
     FUNCITONS USED:
     NONE
     */
    
    if (flag_verbose_correlation_results_ >= 1){
        std::cout << "Normalizing cross-correlation" << std::endl;
    }
    
    for (size_t i = 0; i < correlation_threads_.size(); i++){
        normalized_correlation_ += correlation_threads_[i];
    }
    
    normalized_correlation_.Normalize();

    // noremalize bootstrap realizations
    if (flag_compute_bootstrap_){
        for (size_t i = 0; i < bootstrap_.size(); i ++){
            bootstrap_[i].Normalize();
        }
    }
    
    
}

void CorrelationResults::SaveCrossCorrelation(){
    /**
     EXPLANATION:
     Save the cross correlation results
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     CorrelationResults
     
     FUNCITONS USED:
     ToStr
     */
    
    std::string filename;
    
    // save normalized cross-correlation
    if (flag_verbose_correlation_results_ >= 1){
        std::cout << "Saving cross-correlation" << std::endl;
    }
    if (flag_projection_correction_){
        filename = output_base_name_ + ".data";
        {
            std::ofstream file(filename.c_str(),std::ofstream::trunc);
            if (file.is_open()){
                for (size_t i = 0; i < num_bins_; i++){
                    
                    file << i << " " << normalized_correlation_.xi(i) - normalized_correlation_.xi(i) << std::endl;
                }
                
                file.close();
            }
            else{
                std::cout << "Error : In CorrelationResults::SaveCrossCorrelation : Unable to open file:" << std::endl << filename << std::endl;
            }
        }
        filename = output_base_name_ + ".uncorrected.data";
        {
            std::ofstream file(filename.c_str(),std::ofstream::trunc);
            if (file.is_open()){
                for (size_t i = 0; i < num_bins_; i++){
                    
                    file << i << " " << normalized_correlation_.xi(i) << std::endl;
                }
                
                file.close();
            }
            else{
                std::cout << "Error : In CorrelationResults::SaveCrossCorrelation : Unable to open file:" << std::endl << filename << std::endl;
            }
        }
        filename = output_base_name_ + ".correction.data";
        {
            std::ofstream file(filename.c_str(),std::ofstream::trunc);
            if (file.is_open()){
                for (size_t i = 0; i < num_bins_; i++){
                    
                    file << i << " " << normalized_correlation_.xi_correction(i) << std::endl;
                }
                
                file.close();
            }
            else{
                std::cout << "Error : In CorrelationResults::SaveCrossCorrelation : Unable to open file:" << std::endl << filename << std::endl;
            }
        }
    }
    else{
        filename = output_base_name_ + ".data";
        {
            std::ofstream file(filename.c_str(),std::ofstream::trunc);
            if (file.is_open()){
                for (size_t i = 0; i < num_bins_; i++){
                    
                    file << i << " " << normalized_correlation_.xi(i) << std::endl;
                }
                
                file.close();
            }
            else{
                std::cout << "Error : In CorrelationResults::SaveCrossCorrelation : Unable to open file:" << std::endl << filename << std::endl;
            }
        }
    }
    filename = output_base_name_ + ".full.data";
    {
        std::ofstream file(filename.c_str(),std::ofstream::trunc); 
        if (file.is_open()){
            file << "# bin_index " << CorrelationPlate::InfoHeader() << std::endl;
        
            for (size_t i = 0; i < num_bins_; i++){
            
                file << i << " " << normalized_correlation_.Info(i) << std::endl;
            }
        
            file.close();
        }
        else{
            std::cout << "Error : In CorrelationResults::SaveCrossCorrelation : Unable to open file:" << std::endl << filename << std::endl;
        }
    }
    
    // save grid
    filename = output_base_name_ + ".grid";
    {
        std::ofstream file(filename.c_str(),std::ofstream::trunc);
        if (file.is_open()){
            
            for (size_t i = 0; i < num_bins_; i++){
                
                file << i << " " << normalized_correlation_.mean_pi(i) << " " << normalized_correlation_.mean_sigma(i) << " " << normalized_correlation_.mean_z_in_bin(i) << std::endl;
            }
            
            file.close();
        }
        else{
            std::cout << "Error : In CorrelationResults::SaveCrossCorrelation : Unable to open file:" << std::endl << filename << std::endl;
        }
    }
    
    // save mean redshift
    if (flag_verbose_correlation_results_ >= 1){
        std::cout << "Saving mean z" << std::endl;
    }
    filename = output_base_name_ + ".z";
    {
        std::ofstream file(filename.c_str(),std::ofstream::trunc);
        if (file.is_open()){
            
            file << normalized_correlation_.mean_z();
            
            file.close();
        }
        else{
            std::cout << "Error : In CorrelationResults::SaveCrossCorrelation : Unable to open file:" << std::endl << filename << std::endl;
        }
    }
    // save bootstrap realizations
    if (flag_compute_bootstrap_){
        if (flag_verbose_correlation_results_ >= 1){
            std::cout << "Saving bootstrap realizations" << std::endl;
        }
        for (size_t i = 0; i < bootstrap_.size(); i ++){
            filename = bootstrap_results_ + ".bootstrap" + ToStr(i) + ".data";
            {
                std::ofstream file(filename.c_str(),std::ofstream::trunc); 
                if (file.is_open()){
                    
                    for (size_t j = 0; j < num_bins_; j++){
                        
                        file << j << " " << bootstrap_[i].xi(j) << std::endl;
                    }
                    
                    file.close();
                }
                else{
                    std::cout << "Error : In CorrelationResults::SaveCrossCorrelation : Unable to open file:" << std::endl << filename << std::endl;
                }
            }
            
            filename = bootstrap_results_ + ".bootstrap" + ToStr(i) + ".grid";
            {
                std::ofstream file(filename.c_str(),std::ofstream::trunc);
                if (file.is_open()){
                    
                    for (size_t j = 0; j < num_bins_; j++){
                        
                        file << j << " " << bootstrap_[i].mean_pi(j) << " " << bootstrap_[i].mean_sigma(j) << " " << bootstrap_[i].mean_z_in_bin(j) << std::endl;
                    }
                    
                    file.close();
                }
                else{
                    std::cout << "Error : In CorrelationResults::SaveCrossCorrelation : Unable to open file:" << std::endl << filename << std::endl;
                }
            }
            
        }

    }

}


