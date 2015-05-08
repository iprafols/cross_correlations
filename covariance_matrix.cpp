/**
 covariance_matrix.cpp
 Purpose: This files contains the body for the functions defined in covariance_matrix.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 11/17/2014
 */

#include "covariance_matrix.h"

CovarianceMatrix::CovarianceMatrix(const Input& input, const PlateNeighbours& kPlateNeighbours){
    /**
     EXPLANATION:
     Cosntructs a CovarianceMatrix instance and initializes all its variables
     
     INPUTS:
     input - object of type Input
     kPlateNeighbours - a PlateNeighbours instance containing the list of plates
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovarianceMatrix
     Input
     
     FUNCITONS USED:
     NONE
     */
    
    // set flags from input
    flag_verbose_covariance_matrix_ = input.flag_verbose_covariance_matrix();
    
    // setting the number of bins from input
    num_bins_ = input.num_bins();
    
    // setting the results directory and the pairs file name from input
    output_base_name_ = input.output() + input.output_base_name();

    // initializing covariance matrix, all elements set to 0
    if (input.flag_covariance_matrix_from_file()){
        if (flag_verbose_covariance_matrix_ >= 2){
            std::cout << "Initializig covariance matrix" << std::endl;
        }
        for (size_t i = 0; i < num_bins_; i++){
            for (size_t j = i; j < num_bins_; j++){
                cov_mat_[std::pair<size_t,size_t>(i,j)] = 0.0;
            }
        }
    }
    // initialization of the plates map
    else{
        // initialization of the normalized cross-correlation variable
        normalized_cov_mat_ = CorrelationPlate(input, _NORM_, kPlateNeighbours.GetNeighboursList(_NORM_),true);
        
        plates_list_ = kPlateNeighbours.GetPlatesList();
        for (size_t i = 0; i < plates_list_.size(); i ++){
            covariance_plates_[plates_list_[i]] = CorrelationPlate(input, plates_list_[i], kPlateNeighbours.GetNeighboursList(plates_list_[i]), true);
        }
        skip_plates_ = input.skip_plates();
    }
    
    // initializing bootstrap covariance matrix, all elements set to 0
    if (input.flag_compute_bootstrap()){
        if (flag_verbose_covariance_matrix_ >= 2){
            std::cout << "Initializig bootstrap covariance matrix" << std::endl;
        }
        for (size_t i = 0; i < num_bins_; i++){
            for (size_t j = i; j < num_bins_; j++){
                bootstrap_cov_mat_[std::pair<size_t,size_t>(i,j)] = 0.0;
            }
        }        
    }    
    
}

double CovarianceMatrix::bootstrap_cov_mat(size_t i,size_t j) const{
    /**
     EXPLANATION:
     Access function for bootstrap_cov_mat_
     
     INPUTS:
     i,j - Indexes of the selected bootstrap_cov_mat_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovarianceMatrix
     
     FUNCITONS USED:
     NONE
     */    
    if (i <= j){
        CovMat::const_iterator it = bootstrap_cov_mat_.find(std::pair<size_t,size_t>(i,j));
        if (it == bootstrap_cov_mat_.end()){
            return _BAD_DATA_;
        }
        else{
            return (*it).second;
        }
    }
    else{
        CovMat::const_iterator it = bootstrap_cov_mat_.find(std::pair<size_t,size_t>(j,i));
        if (it == bootstrap_cov_mat_.end()){
            return _BAD_DATA_;
        }
        else{
            return (*it).second;
        }
    }
}

double CovarianceMatrix::cov_mat(size_t i,size_t j) const{
    /**
     EXPLANATION:
     Access function for covmat_
     
     INPUTS:
     i,j - Indexes of the selected cov_mat_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovarianceMatrix
     
     FUNCITONS USED:
     NONE
     */
    if (i <= j){
        CovMat::const_iterator it = cov_mat_.find(std::pair<size_t,size_t>(i,j));
        if (it == cov_mat_.end()){
            return _BAD_DATA_;
        }
        else{
            return (*it).second;
        }
    }
    else{
        CovMat::const_iterator it = cov_mat_.find(std::pair<size_t,size_t>(j,i));
        if (it == cov_mat_.end()){
            return _BAD_DATA_;
        }
        else{
            return (*it).second;
        }
    }
}

void CovarianceMatrix::ComputeBootstrapCovMat(const std::vector<CorrelationPlate>& bootstrap){
    /**
     EXPLANATION:
     Computes the covariance matrix using the bootstrap realizations
     
     INPUTS:
     bootstrap - a vector of CorrelationPlate instances containing the cross-correlation for the different bootstrap realizations
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     CovarianceMatrix
     
     FUNCITONS USED:
     NONE
     */
    
    if (flag_verbose_covariance_matrix_ >= 1){
        std::cout << "Computing the bootstrap covariance matrix" << std::endl;
    }
    
    // checking that we have at least two bootstrap realization
    if (bootstrap.size() < 2){
        std::cout << "There are " << bootstrap.size() << "bootstrap samples. This is not enough to do the calculation, skipping..." << std::endl;
        return;
    }
        
    // computing mean values
    std::vector<double> mean;
    mean.resize(num_bins_);
    for (size_t i = 0; i < num_bins_; i++){
        for (size_t j = 0; j < bootstrap.size(); j++){
            mean[i] += bootstrap[j].xi(i);
        }
        
        mean[i] /= double(bootstrap.size());
    }
    
    // computing covariance matrix
    for (CovMat::iterator it = bootstrap_cov_mat_.begin(); it != bootstrap_cov_mat_.end(); it ++){
        for (size_t j = 0; j < bootstrap.size(); j++){
            (*it).second += (bootstrap[j].xi((*it).first.first)-mean[(*it).first.first])*(bootstrap[j].xi((*it).first.second)-mean[(*it).first.second]);
        }
        (*it).second /= double(bootstrap.size())-1;
    }

    // saving covariance matrix
    SaveBootstrapCovMat();
}

void CovarianceMatrix::ComputeCovMat(const AstroObjectDataset& object_list, const LyaSpectraDataset& spectra_list, const Input& input, const PlateNeighbours& kPlateNeighbours){
    /**
     EXPLANATION:
     Computes the covariance matrix
     
     INPUTS:
     object_list - an AstroObjectDataset instance
     spectra_list - a LyaSpectraDataset instance
     input - object of type Input
     kPlateNeighbours - a PlateNeighbours instance containing the list of plates
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovarianceMatrix
     PlateNeighbours
     
     FUNCITONS USED:
     NONE
     */
    
    if (input.flag_covariance_matrix_from_file()){
        std::vector<int> plates = kPlateNeighbours.GetPlatesList();
        std::vector<double> total_weight;
        total_weight.reserve(num_bins_);
        
        if (flag_verbose_covariance_matrix_ >= 1){
            std::cout << "Computing the covariance matrix" << std::endl;
        }
        if (flag_verbose_covariance_matrix_ >= 2){
            std::cout << std::setprecision(8);
        }
        // load interpolation map
        LyaAutoInterpolationMap lya_auto_correlation_map(input);
        
        // loop over regions (1st index: i)
        for (size_t i = 0; i < num_bins_; i++){
            
            // reading pairs from bin i
            PairDataset pair_dataset_i(input, i, plates);
            std::vector<int> plates_bin_i = pair_dataset_i.GetPlatesList();
            
            // computing weight if necessary
            if (total_weight.size() <= i){
                if (flag_verbose_covariance_matrix_ >= 1){
                    std::cout << "computing weight for bin " << i << "\n";
                }
                total_weight.push_back(ComputeTotalWeight(pair_dataset_i, plates_bin_i));
                if (flag_verbose_covariance_matrix_ >= 1){
                    std::cout << "done\n";
                }
            }
            
            // loop over regions (2nd index: j)
            for (size_t j = i; j < num_bins_; j++){
                
                // checking that the desired covariance matrix element is existent
                CovMat::iterator it = cov_mat_.find(std::pair<int, int>(i, j));
                if (it == cov_mat_.end()){
                    std::cout << "Warning : In CovarianceMatrix::ComputeCovMat : Element " << i << ", " << j << " of the covariance matrix not found. Ignoring..."  << std::endl;
                    continue;
                }
                
                // reading pairs from region j
                PairDataset pair_dataset_j(input, j, plates);
                std::vector<int> plates_bin_j = pair_dataset_j.GetPlatesList();
                
                // computing weight if necessary
                if (total_weight.size() <= j){
                    if (flag_verbose_covariance_matrix_ >= 1){
                        std::cout << "computing weight for bin " << j << "\n";
                    }
                    total_weight.push_back(ComputeTotalWeight(pair_dataset_j, plates_bin_j));
                    if (flag_verbose_covariance_matrix_ >= 1){
                        std::cout << "done\n";
                    }
                }
                
                // computing covariance matrix
                #pragma omp parallel for schedule(dynamic)
                // loop over plates in bin i
                for (size_t plates_i = 0; plates_i < plates_bin_i.size(); plates_i ++){
                    double add,weight;
                    
                    std::vector<Pair> list_i = pair_dataset_i.list(plates_bin_i[plates_i]);
                    if (list_i.size() == 0){
                        continue;
                    }
                    
                    // loop over plates in bin j
                    for (size_t plates_j = 0; plates_j < plates_bin_j.size(); plates_j ++){
                        
                        // check that the plates are neighbours
                        if (kPlateNeighbours.AreNeighbours(plates_bin_i[plates_i], plates_bin_j[plates_j])){
                            
                            std::vector<Pair> list_j = pair_dataset_j.list(plates_bin_j[plates_j]);
                            if (list_j.size() == 0){
                                continue;
                            }
                            
                            // loop over pairs in bin i
                            for (size_t pairs_i = 0; pairs_i < list_i.size(); pairs_i ++){
                                
                                Pair object_i = list_i[pairs_i];
                                if (object_i.pixel_weight() == 0.0){
                                    continue;
                                }
                                
                                double sigma_aux = 2.0*object_i.pixel_dist(); // auxiliar variable to compute sigma values
                                
                                
                                // loop over pairs in bin j
                                for (size_t pairs_j = 0; pairs_j < pair_dataset_j.GetNumberPairs(plates_bin_j[plates_j]); pairs_j ++){
                                    
                                    Pair object_j = list_j[pairs_j];
                                    if (object_j.pixel_weight() == 0.0){
                                        continue;
                                    }
                                    
                                    double cos_theta = object_i.spectrum_angle().CosAngularDistance(object_j.spectrum_angle());
                                    
                                    // compute sigma
                                    double sigma;
                                    if (cos_theta == 1.0){
                                        sigma = 0.0;
                                    }
                                    else{
                                        sigma = sqrt(sigma_aux*(1.0-cos_theta)*object_j.pixel_dist()); // auxiliar variable to compute sigma values
                                    }
                                    
                                    // compute pi
                                    double pi = object_i.pixel_dist() - object_j.pixel_dist();
                                    
                                    // add to covariance matrix
                                    if (sigma < 0.01){
                                        if (object_i.pixel_number() == object_j.pixel_number() and object_i.pixel_weight() != 0.0){
                                            weight = object_i.pixel_weight()*object_j.pixel_weight();
                                            add = pow(1+object_i.pixel_z(),CorrelationPlate::half_gamma())/object_i.pixel_weight()/CorrelationPlate::one_plus_z0_to_the_half_gamma();
                                            
                                            #pragma omp critical(covariance)
                                            {
                                                (*it).second += add*weight;
                                            }
                                        }
                                        /*if (pi >= 0.0){
                                         add = lya_auto_correlation_map.LinearInterpolation(pi);
                                         }
                                         else{
                                         add = lya_auto_correlation_map.LinearInterpolation(-pi);
                                         }
                                         if (add != _BAD_DATA_){
                                         weight = object_i.pixel_weight()*object_j.pixel_weight();
                                         if (object_i.pixel_number() == object_j.pixel_number() and weight != 0.0){
                                         add += 1/object_i.pixel_weight();
                                         }
                                         (*it).second += add*weight;
                                         total_weight += weight;
                                         //std::cout << "prova: " << (*it).second << " " << add << " " << weight << std::endl;
                                         }*/
                                    }
                                }
                            }
                        }
                    }
                }
                
                // normalizing element
                if (total_weight[i] > 0.0 and total_weight[j] > 0.0){
                    (*it).second /= total_weight[i]*total_weight[j];
                }
                if (flag_verbose_covariance_matrix_ >= 2){
                    std::cout << i << " " << j << " " << (*it).second << std::endl;
                }

            }
            
        }
    }
    else{
        
        if (flag_verbose_covariance_matrix_ >= 1){
            #pragma omp critical (cout)
            {
                std::cout << "Computing the covariance matrix" << std::endl;
            }
        }
        // loop over plates
        size_t plates_computed = 0;
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = skip_plates_; i < plates_list_.size(); i++){
            
            PlatesMapSimple<CorrelationPlate>::map::iterator it = covariance_plates_.find(plates_list_[i]);
            
            #pragma omp critical (plates_computed)
            {
                plates_computed ++;
                if (flag_verbose_covariance_matrix_ >= 2 or (flag_verbose_covariance_matrix_ >= 1 and plates_computed == plates_computed/100*100)){
                    #pragma omp critical (cout)
                    {
                        std::cout << plates_computed << " out of " << plates_list_.size() << " plates computed" << std::endl;
                    }
                }
                else{
                    (*it).second.set_flag_verbose_correlation_plate(0);
                }
            }

            // compute covariance matrix in selected plate
            (*it).second.ComputeCovMat(object_list, spectra_list, input);

        }
        
        // normalize covariance matrix
        NormalizeCovMat();
    }
    
    // saving covariance matrix
    SaveCovMat();

    

}

double CovarianceMatrix::ComputeTotalWeight(const PairDataset& pair_dataset, const std::vector<int>& plates_list){
    /**
     EXPLANATION:
     Computes the total weight of a given pairDataset bin
     
     INPUTS:
     pair_dataset - pairDataset to compute the weight of
     plates_list - list of plates contained in pair_dataset
     
     OUTPUTS:
     total_weight - total weight in the dataset
     
     CLASSES USED:
     Pair
     PairDataset
     CovarianceMatrix
     
     FUNCITONS USED:
     NONE
     */
    
    double total_weight = 0.0;
    
    // loop over plates
    for (size_t plate = 0; plate < plates_list.size(); plate ++){
        
        std::vector<Pair> list = pair_dataset.list(plates_list[plate]);
        
        // loop over pairs
        for (size_t pair = 0; pair < list.size(); pair ++){
            total_weight += list[pair].pixel_weight();
        }
    }
}

void CovarianceMatrix::NormalizeCovMat(){
    /**
     EXPLANATION:
     Normalizes the covariance matrix results
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     CovarianceMatrix
     
     FUNCITONS USED:
     NONE
     */
    
    if (flag_verbose_covariance_matrix_ >= 1){
        std::cout << "Normalizing covariance matrix" << std::endl;
    }
    
    for (size_t i = 0; i < plates_list_.size(); i++){
        
        normalized_cov_mat_ += (*covariance_plates_.find(plates_list_[i])).second;
        
    }
    normalized_cov_mat_.Normalize();

    cov_mat_ = normalized_cov_mat_.cov_mat();
}

void CovarianceMatrix::SaveBootstrapCovMat(){
    /**
     EXPLANATION:
     Saves the bootstrap covariance matrix
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     CovarianceMatrix
     
     FUNCITONS USED:
     NONE
     */
    std::string filename;
    
    if (flag_verbose_covariance_matrix_ >= 1){
        std::cout << "Saving bootstrap covariance matrix" << std::endl;
    }
    
    filename = output_base_name_ + ".bootstrap.cov";
    {
        std::ofstream file(filename.c_str(),std::ofstream::trunc); 
        if (file.is_open()){
            if (flag_verbose_covariance_matrix_ >= 2){
                std::cout << "Saving full bootstrap covariance matrix" << std::endl;
            }
            for (CovMat::iterator it = bootstrap_cov_mat_.begin(); it != bootstrap_cov_mat_.end(); it ++){
                file << (*it).first.first << " " << (*it).first.second << " " << (*it).second << std::endl;
            }
            
            file.close();
        }
        else{
            std::cout << "Error : In CovarianceMatrix::SaveBootstrapCovMat : Unable to open file:" << std::endl << filename << std::endl;
        }
    }

    filename = output_base_name_ + ".bootstrap.diag.cov";
    {
        std::ofstream file(filename.c_str(),std::ofstream::trunc); 
        if (file.is_open()){
            if (flag_verbose_covariance_matrix_ >= 2){
                std::cout << "Saving diagonal bootstrap covariance matrix" << std::endl;
            }
            for (size_t i = 0; i < num_bins_; i++){
                CovMat::iterator it = bootstrap_cov_mat_.find(std::pair<size_t,size_t>(i,i));
                file << (*it).first.first << " " << (*it).first.second << " " << (*it).second << std::endl;
            }
            
            file.close();
        }
        else{
            std::cout << "Error : In CovarianceMatrix::SaveBootstrapCovMat : Unable to open file:" << std::endl << filename << std::endl;
        }
    }
    
    
}

void CovarianceMatrix::SaveCovMat(){
    /**
     EXPLANATION:
     Saves the covariance matrix
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     CovarianceMatrix
     
     FUNCITONS USED:
     NONE
     */
    std::string filename;
    
    if (flag_verbose_covariance_matrix_ >= 1){
        std::cout << "Saving covariance matrix" << std::endl;
    }
    
    filename = output_base_name_ + ".cov";
    {
        std::ofstream file(filename.c_str(),std::ofstream::trunc);
        if (file.is_open()){
            if (flag_verbose_covariance_matrix_ >= 2){
                std::cout << "Saving full covariance matrix" << std::endl;
            }
            for (CovMat::iterator it = cov_mat_.begin(); it != cov_mat_.end(); it ++){
                file << (*it).first.first << " " << (*it).first.second << " " << (*it).second << std::endl;
            }
            
            file.close();
        }
        else{
            std::cout << "Error : In CovarianceMatrix::SaveCovMat : Unable to open file:" << std::endl << filename << std::endl;
        }
    }
    
}

    