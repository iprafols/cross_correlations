/**
 metal_grid.cpp
 Purpose: This files contains the body for the functions defined in metal_grid.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 28/09/2017
 */

#include "correlation_plate.h"

MetalGrid::MetalGrid(int bad_data){
    /**
     EXPLANATION:
     Constructs a MetalGrid instance and initializes all its variables
     
     INPUTS:
     bad_data - an int valued _BAD_DATA_INT_
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     MetalGrid
     
     FUNCITONS USED:
     ToStr
     */
    if (bad_data != _BAD_DATA_INT_){
        #pragma omp critical (cout)
        {
        std::cout << "Error while initializing a MetalGrid 'bad data' instance" << std::endl;
        }
        std::exit(EXIT_FAILURE);
    }
    
    num_bins_ = _BAD_DATA_INT_;
    
}

MetalGrid::MetalGrid(const size_t& num_bins, const double& metal_wl, const std::string& metal_name, const size_t& flag_verbose_metal_grid){
    /**
     EXPLANATION:
     Constructs a MetalGrid instance and initializes all its variables
     
     INPUTS:
     num_bins - a size_t with the number of bins in the grid
     metal_wl - a double with the metal's wavelength
     metal_name a string with the metal's name
     flag_verbose_metal_grid - a size_t indicating the degree of verbosity of the instance
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     MetalGrid
     
     FUNCITONS USED:
     ToStr
     */
    
    // set flags from input
    flag_verbose_metal_grid_ = flag_verbose_metal_grid;
    
    num_bins_ = num_bins;
    metal_wl_ = metal_wl;
    metal_name_ = metal_name;
    
    // initialize grid
    if (flag_verbose_metal_grid_ >= 3){
        #pragma omp critical (cout)
        {
            std::cout << "MetalGrid: initialize grid" << std::endl;
        }
    }
    mean_pi_.resize(num_bins_,0.0);
    mean_sigma_.resize(num_bins_,0.0);
    mean_z_.resize(num_bins_, 0.0);
    weight_.resize(num_bins_,0.0);
}

double MetalGrid::mean_pi(size_t index) const {
    /**
     EXPLANATION:
     Access function for mean_pi_
     
     INPUTS:
     index - index of the selected mean_pi_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     MetalGrid
     
     FUNCITONS USED:
     NONE
     */

    if (index < mean_pi_.size()){
        return mean_pi_[index];
    }
    else{
        return _BAD_DATA_;
    }
}

double MetalGrid::mean_sigma(size_t index) const {
    /**
     EXPLANATION:
     Access function for mean_sigma_
     
     INPUTS:
     index - index of the selected mean_sigma_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     MetalGrid
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < mean_sigma_.size()){
        return mean_sigma_[index];
    }
    else{
        return _BAD_DATA_;
    }
}

double MetalGrid::mean_z(size_t index) const {
    /**
     EXPLANATION:
     Access function for mean_z_
     
     INPUTS:
     index - index of the selected mean_z_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     MetalGrid
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < mean_z_.size()){
        return mean_z_[index];
    }
    else{
        return _BAD_DATA_;
    }
}

double MetalGrid::weight(size_t index) const {
    /**
     EXPLANATION:
     Access function for weight_
     
     INPUTS:
     index - index of the selected weight_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     MetalGrid
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < weight_.size()){
        return weight_[index];
    }
    else{
        return _BAD_DATA_;
    }
}

void MetalGrid::set_mean_pi(size_t index, double value){
    /**
     EXPLANATION:
     Set function for mean_pi_
     
     INPUTS:
     index - index of the selected mean_pi_ element
     value - element's new value
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     MetalGrid
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < mean_pi_.size()){
        mean_pi_[index] = value;
    }
    else{
        #pragma omp critical (cout)
        {
        std::cout << "Warining: in MetalGrid::set_mean_pi(index, value): The given index is out of bouds, ignoring..." << std::endl;
        }
    }
}

void MetalGrid::set_mean_sigma(size_t index, double value){
    /**
     EXPLANATION:
     Set function for mean_sigma_
     
     INPUTS:
     index - index of the selected mean_sigma_ element
     value - element's new value
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     MetalGrid
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < mean_sigma_.size()){
        mean_sigma_[index] = value;
    }
    else{
        #pragma omp critical (cout)
        {
        std::cout << "Warining: in MetalGrid::set_mean_sigma(index, value): The given index is out of bouds, ignoring..." << std::endl;
        }
    }
}

void MetalGrid::set_mean_z(size_t index, double value){
    /**
     EXPLANATION:
     Set function for mean_z_
     
     INPUTS:
     index - index of the selected mean_sigma_ element
     value - element's new value
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     MetalGrid
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < mean_z_.size()){
        mean_z_[index] = value;
    }
    else{
        #pragma omp critical (cout)
        {
        std::cout << "Warining: in MetalGrid::set_mean_z_in_bin(index, value): The given index is out of bouds, ignoring..." << std::endl;
        }
    }
}

void MetalGrid::set_weight(size_t index, double value){
    /**
     EXPLANATION:
     Set function for weight_
     
     INPUTS:
     index - index of the selected weight_ element
     value - element's new value
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     MetalGrid
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < weight_.size()){
        weight_[index] = value;
    }
    else{
        #pragma omp critical (cout)
        {
        std::cout << "Warining: in MetalGrid::set_weight(index, value): The given index is out of bouds, ignoring..." << std::endl;
        }
    }
}

void MetalGrid::AddPair(const size_t& k_index, const double& weight, const double& pi, const double& sigma, const double& z){
    /**
     EXPLANATION:
     Adds pair contribution to the specified bin
     
     INPUTS:
     k_index - an integer specifying the bin to add the contribution to
     weight - an double specifying the weight of the pair
     pi - a double specifying the parallel separation of the pair
     sigma - a double specifying the perpendicular separation of the pair
     z - a double specifying the average redshift of the pair
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     MetalGrid
     LyaPixel
     
     FUNCITONS USED:
     NONE
     */
    
    if (k_index > mean_pi_.size()){
        #pragma omp critical (cout)
        {
        std::cout << "In function MetalGrid::AddPair, the given index is out of bounds. Ignoring..." << std::endl;
        }
    }
    
    mean_pi_[k_index] += pi*weight;
    mean_sigma_[k_index] += sigma*weight;
    mean_z_[k_index] += z*weight;
    weight_[k_index] += weight;
}

void MetalGrid::Normalize(){
    /**
     EXPLANATION:
     Normalizes the values in the grid
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     MetalGrid
     
     FUNCITONS USED:
     NONE
     */
    
    for (size_t i = 0; i < num_bins_; i ++){
        
        if (weight_[i] == 0.0){
            // if the weight is zero, complain and exit the function
            std::cerr << "Zero Division Error in MetalGrid::Normalize. Ignoring..." << std::endl;
            continue;
        }
        // if the weight is 1.0, normalization is not needed
        else if (weight_[i] != 1.0){
            
            #pragma omp critical (normalize)
            {
            mean_pi_[i] /= weight_[i];
            mean_sigma_[i] /= weight_[i];
            mean_z_[i] /= weight_[i];
            weight_[i] = 1.0;
            }
        }
        
    }
}

void MetalGrid::operator+= (const MetalGrid& other){
    /**
     EXPLANATION:
     Overloads the += operator. Updates the mean_sigma_, mean_pi_, mean_z_, and weight_ values by adding those in other.
     
     INPUTS:
     other - a MetalGrid instance to be added
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     MetalGrid
     
     FUNCITONS USED:
     NONE
     */
    
    // check that both instances have the same number of bins
    if (num_bins_ != other.num_bins()){
        #pragma omp critical (cout)
        {
            std::cout << "Warning : In MetalGrid::operator+= : Trying to add MetalGrid with different number of bins. Ignoring..." << std::endl;
        }
        return;
    }
    
    for (size_t i = 0; i < num_bins_; i ++){
        mean_pi_[i] += other.mean_pi(i);
        mean_sigma_[i] += other.mean_sigma(i);
        mean_z_[i] += other.mean_z(i);
        weight_[i] += other.weight(i);
    }
}

MetalGrid MetalGrid::operator- (const MetalGrid& other){
    /**
     EXPLANATION:
     Overloads the - operator. Returns a new instance with the mean_sigma_, mean_pi_, mean_z_, and weight_ values subtracted by those in other.
     
     INPUTS:
     other - a MetalGrid instance to be subtracted
     
     OUTPUTS:
     temp - a MetalGrid instance
     
     CLASSES USED:
     MetalGrid
     
     FUNCITONS USED:
     NONE
     */

    MetalGrid temp;
    temp = MetalGrid(num_bins_, metal_wl_, metal_name_, flag_verbose_metal_grid_);
    
    // check that both instances have the same number of bins
    if (num_bins_ != other.num_bins()){
        #pragma omp critical (cout)
        {
            std::cout << "Warning : In MetalGrid::operator- : Trying to add MetalGrid with different number of bins. Returning zero filled MetalGrid..." << std::endl;
        }
        return temp;
    }
    
    for (size_t i = 0; i < num_bins_; i ++){
        temp.set_mean_pi(i, mean_pi_[i] - other.mean_pi(i));
        temp.set_mean_sigma(i, mean_sigma_[i] - other.mean_sigma(i));
        temp.set_mean_z(i, mean_z_[i] - other.mean_z(i));
        temp.set_weight(i, weight_[i] - other.weight(i));
    }
    
    return temp;
}
        
MetalGrid MetalGrid::operator* (const MetalGrid& other){
    /**
     EXPLANATION:
     Overloads the * operator. Returns a new instance with the mean_sigma_, mean_pi_, mean_z, and weight_ values multiplied by those in other.
     
     INPUTS:
     other - a MetalGrid instance to be subtracted
     
     OUTPUTS:
     temp - a MetalGrid instance
     
     CLASSES USED:
     MetalGrid
     
     FUNCITONS USED:
     NONE
     */
    MetalGrid temp;
    temp = MetalGrid(num_bins_, metal_wl_, metal_name_, flag_verbose_metal_grid_);
    
    // check that both instances have the same number of bins
    if (num_bins_ != other.num_bins()){
        #pragma omp critical (cout)
        {
            std::cout << "Warning : In MetalGRid::operator* : Trying to add MetalGrid with different number of bins. Returning zero filled MetalGrid..." << std::endl;
        }
        return temp;
    }
    
    for (size_t i = 0; i < num_bins_; i ++){
        temp.set_mean_pi(i, mean_pi_[i]*other.mean_pi(i));
        temp.set_mean_sigma(i, mean_sigma_[i]*other.mean_sigma(i));
        temp.set_mean_z(i, mean_z_[i]*other.mean_z(i));
        temp.set_weight(i, weight_[i]*other.weight(i));
    }
    
    return temp;
}

