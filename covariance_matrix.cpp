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


    // initialization of the plates map
    if (flag_verbose_covariance_matrix_ >= 2){
        std::cout << "Initializig covariance matrix" << std::endl;
    }
    
    // initialization of the normalized cross-correlation variable
    normalized_cov_mat_ = CovariancePlate(input, _NORM_, kPlateNeighbours.GetNeighboursList(_NORM_));
    
    plates_list_ = kPlateNeighbours.GetPlatesList();
    int num_threads = atoi(std::getenv("OMP_NUM_THREADS"));
    covariance_threads_.reserve(num_threads);
    for (size_t i = 0; i < num_threads; i ++){
        covariance_threads_.push_back(CovariancePlate(input, _NORM_, kPlateNeighbours.GetNeighboursList(_NORM_)));
    }
    skip_plates_ = input.skip_plates();
    
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
     bootstrap - a vector of CovariancePlate instances containing the cross-correlation for the different bootstrap realizations
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovariancePlate
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

void CovarianceMatrix::ComputeCovMat(const AstroObjectDataset& object_list, const SpectraDataset& spectra_list, const Input& input, const PlateNeighbours& kPlateNeighbours){
    /**
     EXPLANATION:
     Computes the covariance matrix
     
     INPUTS:
     object_list - an AstroObjectDataset instance
     spectra_list - a SpectraDataset instance
     input - object of type Input
     kPlateNeighbours - a PlateNeighbours instance containing the list of plates
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovarianceMatrix
     PlateNeighbours
     SpectraDataset
     
     FUNCITONS USED:
     NONE
     */
        
    // compute the 1D lyman-alpha auto-correlation
    std::vector<LyaAutoInterpolationMap> lya_auto;
    for (size_t i = 0; i <= input.pixels_separation(); i++){
        lya_auto.push_back(LyaAutoInterpolationMap(input, i));
    }
    
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
        
        CovariancePlate plate (input, plates_list_[i], kPlateNeighbours.GetNeighboursList(plates_list_[i]));
        
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
                plate.set_flag_verbose_covariance_plate(0);
            }
        }

        // compute covariance matrix in selected plate
        plate.ComputeCovMat(object_list, spectra_list, input, lya_auto);
        
        // add to total value
        int thread_num = omp_get_thread_num();
        covariance_threads_[thread_num] += plate;

    }
    
    // normalize covariance matrix
    NormalizeCovMat();
    
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
     CovariancePlate
     CovarianceMatrix
     
     FUNCITONS USED:
     NONE
     */
    
    if (flag_verbose_covariance_matrix_ >= 1){
        std::cout << "Normalizing covariance matrix" << std::endl;
    }
    
    /*for (size_t i = 0; i < plates_list_.size(); i++){
        
        normalized_cov_mat_ += (*covariance_plates_.find(plates_list_[i])).second;
        
    }*/
    for (size_t i = 0; i < covariance_threads_.size(); i++){
        
        normalized_cov_mat_ += covariance_threads_[i];
        
    }
    
    // test tocheck the total weights of the bins
    if (flag_verbose_covariance_matrix_ >= 3){
        std::cout << "TEST TO CHECK THE TOTAL WEIGHTS OF THE BINS: compare the values with those of *.full.data" << std::endl;
        for (size_t i = 0; i < num_bins_; i++){
            std::cout << i << " " << normalized_cov_mat_.weight(i) << std::endl;
        }
        std::cout << "END OF TEST" << std::endl;
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
     CovariancePlate
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
     CovariancePlate
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

    