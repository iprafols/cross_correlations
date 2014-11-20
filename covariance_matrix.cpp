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
    
    // setting the number of bins and plates from input
    num_bins_ = input.num_bins();
    
    // output settings
    output_base_name_ = input.output() + input.output_base_name();

    // initializing covariance matrix, all elements set to 0
    std::cout << "TEST: initializig covariance matrix" << std::endl;
    for (size_t i = 0; i < num_bins_; i++){
        for (size_t j = i; j < num_bins_; j++){
            cov_mat_[std::pair<size_t,size_t>(i,j)] = 0.0;
        }
    }
    
    // initializing bootstrap covariance matrix, all elements set to 0
    if (input.flag_compute_bootstrap()){
        std::cout << "TEST: initializig bootstrap covariance matrix" << std::endl;
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
        return (*bootstrap_cov_mat_.find(std::pair<size_t,size_t>(i,j))).second;
    }
    else{
        return (*bootstrap_cov_mat_.find(std::pair<size_t,size_t>(j,i))).second;
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
        return (*cov_mat_.find(std::pair<size_t,size_t>(i,j))).second;
    }
    else{
        return (*cov_mat_.find(std::pair<size_t,size_t>(j,i))).second;
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
    std::cout << "Computing the covariance matrix using the bootstrap realizations" << std::endl;
    
    // checking that we have at least two bootstrap realization
    if (bootstrap.size() < 2){
        std::cout << "There are " << bootstrap.size() << "bootstrap samples. This is not enough to do the calculation, skipping..." << std::endl;
        return;
    }
        
    // computing mean values
    //std::cout << "TEST: printing num_bins_: " << num_bins_ << std::endl;
    //std::cout << "TEST: computing mean values, enter character to proceed" << std::endl;
    std::vector<double> mean;
    //std::cout << "TEST: vector created to store mean values, enter character to proceed" << std::endl;
    //std::cin >> test;
    mean.resize(num_bins_);
    //std::cout << "TEST: vector resized" << std::endl;
    //std::cout << "TEST: bootstrap.size() = " << bootstrap.size() << std::endl;
    for (size_t i = 0; i < num_bins_; i++){
        for (size_t j = 0; j < bootstrap.size(); j++){
            mean[i] += bootstrap[j].xi(i);
        }
        
        mean[i] /= double(bootstrap.size());
    }
    
    // computing covariance matrix
    //std::cout << "TEST: computing bootstrap covariance matriz" << std::endl;
    for (CovMat::iterator it = bootstrap_cov_mat_.begin(); it != bootstrap_cov_mat_.end(); it ++){
        for (size_t j = 0; j < bootstrap.size(); j++){
            (*it).second += (bootstrap[j].xi((*it).first.first)-mean[(*it).first.first])*(bootstrap[j].xi((*it).first.second)-mean[(*it).first.second]);
        }
        (*it).second /= double(bootstrap.size())-1;
    }

    //std::cout << "TEST: saving bootstrap covariance matrix" << std::endl;
    SaveBootstrapCovMat();
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
    
    filename = output_base_name_ + ".bootstrap.cov";
    {
        std::ofstream file(filename.c_str(),std::ofstream::trunc); 
        if (file.is_open()){
            std::cout << "TEST: saving full bootstrap covariance matrix" << std::endl;
            for (CovMat::iterator it = bootstrap_cov_mat_.begin(); it != bootstrap_cov_mat_.end(); it ++){
                file << (*it).first.first << " " << (*it).first.second << " " << (*it).second << std::endl;
            }
            
            file.close();
        }
        else{
            std::cout << "Unable to open file:" << std::endl << filename << std::endl;
        }
    }

    filename = output_base_name_ + ".bootstrap.diag.cov";
    {
        std::ofstream file(filename.c_str(),std::ofstream::trunc); 
        if (file.is_open()){
            std::cout << "TEST: saving diagonal bootstrap covariance matrix" << std::endl;
            for (size_t i = 0; i < num_bins_; i++){
                CovMat::iterator it = bootstrap_cov_mat_.find(std::pair<size_t,size_t>(i,i));
                file << (*it).first.first << " " << (*it).first.second << " " << (*it).second << std::endl;
            }
            
            file.close();
        }
        else{
            std::cout << "Unable to open file:" << std::endl << filename << std::endl;
        }
    }
    
    
}





    