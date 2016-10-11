/**
 distortion_matrix.cpp
 Purpose: This files contains the body for the functions defined in distortion_matrix.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 11/17/2014
 */

#include "distortion_matrix.h"

DistortionMatrix::DistortionMatrix(const Input& input, const PlateNeighbours& kPlateNeighbours){
    /**
     EXPLANATION:
     Cosntructs a DistortionMatrix instance and initializes all its variables
     
     INPUTS:
     input - object of type Input
     kPlateNeighbours - a PlateNeighbours instance containing the list of plates
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     DistortionMatrix
     Input
     
     FUNCITONS USED:
     NONE
     */
    
    // set flags from input
    flag_verbose_distortion_matrix_ = input.flag_verbose_distortion_matrix();
    
    // setting the number of bins from input
    num_bins_ = input.num_bins();
    
    // setting the results directory and the pairs file name from input
    output_base_name_ = input.output() + input.output_base_name();

    if (flag_verbose_distortion_matrix_ >= 2){
        std::cout << "Initializig distortion matrix" << std::endl;
    }
    
    // initialization of the normalized cross-correlation variable
    normalized_dist_mat_ = DistortionPlate(input, _NORM_, kPlateNeighbours.GetNeighboursList(_NORM_));
    
    plates_list_ = kPlateNeighbours.GetPlatesList();
    int num_threads = atoi(std::getenv("OMP_NUM_THREADS"));
    distortion_threads_.reserve(num_threads);
    for (size_t i = 0; i < num_threads; i ++){
        distortion_threads_.push_back(DistortionPlate(input, _NORM_, kPlateNeighbours.GetNeighboursList(_NORM_)));
    }
    skip_plates_ = input.skip_plates();
    
}

double DistortionMatrix::dist_mat(size_t i,size_t j) const{
    /**
     EXPLANATION:
     Access function for dist_mat_
     
     INPUTS:
     i,j - Indexes of the selected dist_mat_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     DistortionMatrix
     
     FUNCITONS USED:
     NONE
     */
    CovMat::const_iterator it = dist_mat_.find(std::pair<size_t,size_t>(i,j));
    if (it == dist_mat_.end()){
        return _BAD_DATA_;
    }
    else{
        return (*it).second;
    }
}

void DistortionMatrix::ComputeDistMat(const AstroObjectDataset& object_list, const SpectraDataset& spectra_list, const Input& input, const PlateNeighbours& kPlateNeighbours){
    /**
     EXPLANATION:
     Computes the distortion matrix
     
     INPUTS:
     object_list - an AstroObjectDataset instance
     spectra_list - a SpectraDataset instance
     input - object of type Input
     kPlateNeighbours - a PlateNeighbours instance containing the list of plates
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     DistortionMatrix
     PlateNeighbours
     SpectraDataset
     
     FUNCITONS USED:
     NONE
     */

    if (flag_verbose_distortion_matrix_ >= 1){
        #pragma omp critical (cout)
        {
            std::cout << "Computing the distortion matrix" << std::endl;
        }
    }
    // loop over plates
    size_t plates_computed = 0;
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = skip_plates_; i < plates_list_.size(); i++){
        
        DistortionPlate plate (input, plates_list_[i], kPlateNeighbours.GetNeighboursList(plates_list_[i]));
        
        #pragma omp critical (plates_computed)
        {
            plates_computed ++;
            if (flag_verbose_distortion_matrix_ >= 2 or (flag_verbose_distortion_matrix_ >= 1 and plates_computed == plates_computed/100*100)){
                #pragma omp critical (cout)
                {
                    std::cout << plates_computed << " out of " << plates_list_.size() << " plates computed" << std::endl;
                }
            }
            else{
                plate.set_flag_verbose_distortion_plate(0);
            }
        }

        // compute distortion matrix in selected plate
        plate.ComputeDistMat(object_list, spectra_list, input);
        
        // add to total value
        int thread_num = omp_get_thread_num();
        distortion_threads_[thread_num] += plate;

    }
        
    // normalize distortion matrix
    NormalizeDistMat();
    
    // saving distortion matrix
    SaveDistMat();

}

void DistortionMatrix::NormalizeDistMat(){
    /**
     EXPLANATION:
     Normalizes the distortion matrix results
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovariancePlate
     DistortionMatrix
     
     FUNCITONS USED:
     NONE
     */
    
    if (flag_verbose_distortion_matrix_ >= 1){
        std::cout << "Normalizing distortion matrix" << std::endl;
    }
    
    for (size_t i = 0; i < distortion_threads_.size(); i++){
        
        normalized_dist_mat_ += distortion_threads_[i];
        
    }
    
    normalized_dist_mat_.Normalize();

    dist_mat_ = normalized_dist_mat_.dist_mat();
}

void DistortionMatrix::SaveDistMat(){
    /**
     EXPLANATION:
     Saves the distortion matrix
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovariancePlate
     DistortionMatrix
     
     FUNCITONS USED:
     NONE
     */
    std::string filename;
    
    if (flag_verbose_distortion_matrix_ >= 1){
        std::cout << "Saving distortion matrix" << std::endl;
    }
    
    filename = output_base_name_ + ".dmat";
    {
        std::ofstream file(filename.c_str(),std::ofstream::trunc);
        if (file.is_open()){
            for (CovMat::iterator it = dist_mat_.begin(); it != dist_mat_.end(); it ++){
                file << (*it).first.first << " " << (*it).first.second << " " << (*it).second << std::endl;
            }
            
            file.close();
        }
        else{
            std::cout << "Error : In DistortionMatrix::SaveDistMat : Unable to open file:" << std::endl << filename << std::endl;
        }
    }
    
}

    