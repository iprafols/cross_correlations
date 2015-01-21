/**
 covariance_matrix.cpp
 Purpose: This files contains the body for the functions defined in covariance_matrix.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 11/17/2014
 */

#include "covariance_matrix.h"

CovarianceMatrix::CovarianceMatrix(const Input& input){
    /**
     EXPLANATION:
     Cosntructs a CovarianceMatrix instance and initializes all its variables
     
     INPUTS:
     input - object of type Input
     
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
    if (flag_verbose_covariance_matrix_ >= 2){
        std::cout << "Initializig covariance matrix" << std::endl;
    }
    for (size_t i = 0; i < num_bins_; i++){
        for (size_t j = i; j < num_bins_; j++){
            cov_mat_[std::pair<size_t,size_t>(i,j)] = 0.0;
        }
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

void CovarianceMatrix::ComputeCovMat(const Input& input, const PlateNeighbours& kPlateNeighbours){
    /**
     EXPLANATION:
     Computes the covariance matrix
     
     INPUTS:
     kPlateNeighbours - a PlateNeighbours instance containing the list of plates
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovarianceMatrix
     PlateNeighbours
     
     FUNCITONS USED:
     NONE
     */
    
    std::vector<int> plates = kPlateNeighbours.GetPlatesList();
    
    if (flag_verbose_covariance_matrix_ >= 1){
        std::cout << "Computing the covariance matrix" << std::endl;
    }
    // load interpolation map
    LyaAutoInterpolationMap lya_auto_correlation_map(input);
    double add;
    
    // loop over regions (1st index: i)
    for (size_t i = 0; i < num_bins_; i++){
        
        // reading pairs from bin i
        PairDataset pair_dataset_i(input, i, plates);
        std::vector<int> plates_bin_i = pair_dataset_i.GetPlatesList();
        
        // loop over regions (2n index: j)
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
            
            // computing covariance matrix
            // loop over plates in bin i
            for (size_t plates_i = 0; plates_i < plates_bin_i.size(); plates_i ++){
                // loop over plates in bin j
                for (size_t plates_j = 0; plates_j < plates_bin_j.size(); plates_j ++){
                    
                    // check that the plates are neighbours
                    if (kPlateNeighbours.AreNeighbours(plates_bin_i[plates_i], plates_bin_j[plates_j])){
                        
                        // loop over pairs in bin i
                        for (size_t pairs_i = 0; pairs_i < pair_dataset_i.GetNumberPairs(plates_bin_i[plates_i]); pairs_i ++){
                            
                            Pair object_i = pair_dataset_i.list(plates_bin_i[plates_i], pairs_i);
                            
                            double sigma_aux = 2.0*object_i.pixel_dist(); // auxiliar variable to compute sigma values
                            
                            
                            // loop over pairs in bin j
                            for (size_t pairs_j = 0; pairs_j < pair_dataset_j.GetNumberPairs(plates_bin_j[plates_j]); pairs_j ++){
                                
                                Pair object_j = pair_dataset_j.list(plates_bin_j[plates_j], pairs_j);
                                
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
                                    if (pi >= 0.0){
                                        add = lya_auto_correlation_map.LinearInterpolation(pi);
                                    }
                                    else{
                                       add =  lya_auto_correlation_map.LinearInterpolation(-pi);
                                    }
                                    if (add == _BAD_DATA_){
                                        (*it).second += add;
                                        if (object_i.pixel_number() == object_j.pixel_number()){
                                            (*it).second += 1/object_i.pixel_weight();
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
    }
    
    // saving covariance matrix
    //SaveCovMat();
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





    