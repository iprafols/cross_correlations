/**
 correlation_plate.cpp
 Purpose: This files contains the body for the functions defined in correlation_plate.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 10/02/2014
 */

#include "correlation_plate.h"

CorrelationPlate::CorrelationPlate(int bad_data){
    /**
     EXPLANATION:
     Cosntructs a CorrelationPlate instance and initializes all its variables
     
     INPUTS:
     bad_data - a double valued _BAD_DATA_INT_
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     ToStr
     */
    if (bad_data != _BAD_DATA_INT_){
        #pragma omp critical (cout)
        {
        std::cout << "Error while initializing a CorrelationPlate 'bad data' instance" << std::endl;
        }
        std::exit(EXIT_FAILURE);
    }
    
    plate_number_ = _BAD_DATA_INT_;
    
}

CorrelationPlate::CorrelationPlate(const Input& input, const int plate_number, const std::vector<int>& plate_neighbours){
    /**
     EXPLANATION:
     Cosntructs a CorrelationPlate instance and initializes all its variables
     
     INPUTS:
     input - a Input instance
     plate_number - an integer with the plate number
     plate_neighbours - a vector containing the plate numbers of the neighbouring plates
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     ToStr
     */
    
    // set flags from input
    flag_verbose_correlation_plate_ = input.flag_verbose_correlation_plate();
    flag_write_partial_results_ = input.flag_write_partial_results();
    
    plate_number_ = plate_number;
    plate_neighbours_ = plate_neighbours;
    num_bins_ = input.num_bins();
    if (plate_number_ == _NORM_){
        results_ = "";
        pairs_file_name_ = "";
    }
    else{
        results_ = input.detailed_results();
        pairs_file_name_ = ToStr(plate_number_);
    }
    
    // initialize cross-correlation computations
    if (flag_verbose_correlation_plate_ >= 3){
        #pragma omp critical (cout)
        {
            std::cout << "CorrelationPlate: initialize cross-correlation" << std::endl;
        }
    }
    xi_.resize(num_bins_,0.0);
    mean_pi_.resize(num_bins_,0.0);
    mean_sigma_.resize(num_bins_,0.0);
    mean_z_in_bin_.resize(num_bins_, 0.0);
    weight_.resize(num_bins_,0.0);
    num_averaged_pairs_.resize(num_bins_,0);
    mean_z_ = 0.0;
    weight_z_ = 0.0;
    
    // pair storage settings
    if (flag_write_partial_results_ >= 2){
        if (flag_verbose_correlation_plate_ >= 3){
            #pragma omp critical (cout)
            {
                std::cout << "CorrelationPlate: pair storage settings" << std::endl;
            }
        }
        max_pairs_ = 100;
        position_.resize(num_bins_, 0);
        if (flag_verbose_correlation_plate_ >= 3){
            #pragma omp critical (cout)
            {
            std::cout << "CorrelationPlate: reserving space" << std::endl;
            }
        }
        
        std::vector<Pair> aux;
        aux.resize(max_pairs_, Pair(_BAD_DATA_));
        pairs_information_.resize(num_bins_, aux);
        
    }
    
    // correction to the projection of deltas
    flag_projection_correction_ = input.flag_projection_correction();
    if (flag_projection_correction_){
        // initialize cross-correlation
        xi_correction_.resize(num_bins_,0.0);
    }
}

CorrelationPlate::CorrelationPlate(const int plate_number, const int num_bins, const std::string& results, const std::string& pairs_file_name, const std::vector<int>& plate_neighbours, const size_t& flag_verbose_correlation_plate, const size_t& flag_write_partial_results, const bool& flag_projection_correction){
    /**
     EXPLANATION:
     Cosntructs a CorrelationPlate instance and initializes all its variables
     
     INPUTS:
     input - a Input instance
     results - name of the folder where detailed information will be stored
     plate_number - an integer with the plate number
     plate_neighbours - a vector containing the plate numbers of the neighbouring plates
     flag_verbose_correlation_plate - correlation_plate verbose flag
     flag_write_partial_results - flag to write partial results
     flag_projection_correction - flag to compute the correction to the cross-correlation due to the projection of the delta field
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     ToStr
     */
    
    flag_verbose_correlation_plate_ = flag_verbose_correlation_plate;
    flag_write_partial_results_ = flag_write_partial_results;
    
    plate_number_ = plate_number;
    plate_neighbours_ = plate_neighbours;
    num_bins_ = num_bins;
    if (plate_number_ == _NORM_){
        results_ = "";
        pairs_file_name_ = "";
    }
    else{
        results_ = results;
        pairs_file_name_ = "detailed_info_plate_" + ToStr(plate_number_);
    }
    
    // initialize cross-correlation computations
    xi_.resize(num_bins_,0.0);
    mean_pi_.resize(num_bins_,0.0);
    mean_sigma_.resize(num_bins_,0.0);
    mean_z_in_bin_.resize(num_bins_, 0.0);
    weight_.resize(num_bins_,0.0);
    num_averaged_pairs_.resize(num_bins_,0);
    mean_z_ = 0.0;
    weight_z_ = 0.0;
    
    // correction to the projection of deltas
    flag_projection_correction_ = flag_projection_correction;
    if (flag_projection_correction_){
        // initialize cross-correlation
        xi_correction_.resize(num_bins_,0.0);
    }
    
}

double CorrelationPlate::mean_pi(size_t index) const {
    /**
     EXPLANATION:
     Access function for mean_pi_
     
     INPUTS:
     index - index of the selected mean_pi_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
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

double CorrelationPlate::mean_sigma(size_t index) const {
    /**
     EXPLANATION:
     Access function for mean_sigma_
     
     INPUTS:
     index - index of the selected mean_sigma_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
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

double CorrelationPlate::mean_z_in_bin(size_t index) const {
    /**
     EXPLANATION:
     Access function for mean_z_in_bin_
     
     INPUTS:
     index - index of the selected mean_sigma_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < mean_z_in_bin_.size()){
        return mean_z_in_bin_[index];
    }
    else{
        return _BAD_DATA_;
    }
}

int CorrelationPlate::num_averaged_pairs(size_t index) const {
    /**
     EXPLANATION:
     Access function for num_average_pairs_
     
     INPUTS:
     index - index of the selected num_average_pairs_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < num_averaged_pairs_.size()){
        return num_averaged_pairs_[index];
    }
    else{
        return _BAD_DATA_INT_;
    }
}

double CorrelationPlate::weight(size_t index) const {
    /**
     EXPLANATION:
     Access function for weight_
     
     INPUTS:
     index - index of the selected weight_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
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

double CorrelationPlate::xi(size_t index) const {
    /**
     EXPLANATION:
     Access function for xi_
     
     INPUTS:
     index - index of the selected xi_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < xi_.size()){
        return xi_[index];
    }
    else{
        return _BAD_DATA_;
    }
}

double CorrelationPlate::xi_correction(size_t index) const {
    /**
     EXPLANATION:
     Access function for xi_correction_
     
     INPUTS:
     index - index of the selected xi_correction_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < xi_correction_.size()){
        return xi_correction_[index];
    }
    else{
        return _BAD_DATA_;
    }
}


void CorrelationPlate::set_mean_pi(size_t index, double value){
    /**
     EXPLANATION:
     Set function for mean_pi_
     
     INPUTS:
     index - index of the selected mean_pi_ element
     value - element's new value
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < mean_pi_.size()){
        mean_pi_[index] = value;
    }
    else{
        #pragma omp critical (cout)
        {
        std::cout << "Warining: in CorrelationPlate::set_mean_pi(index, value): The given index is out of bouds, ignoring..." << std::endl;
        }
    }
}

void CorrelationPlate::set_mean_sigma(size_t index, double value){
    /**
     EXPLANATION:
     Set function for mean_sigma_
     
     INPUTS:
     index - index of the selected mean_sigma_ element
     value - element's new value
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < mean_sigma_.size()){
        mean_sigma_[index] = value;
    }
    else{
        #pragma omp critical (cout)
        {
        std::cout << "Warining: in CorrelationPlate::set_mean_sigma(index, value): The given index is out of bouds, ignoring..." << std::endl;
        }
    }
}

void CorrelationPlate::set_mean_z_in_bin(size_t index, double value){
    /**
     EXPLANATION:
     Set function for mean_z_in_bin_
     
     INPUTS:
     index - index of the selected mean_sigma_ element
     value - element's new value
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < mean_sigma_.size()){
        mean_z_in_bin_[index] = value;
    }
    else{
        #pragma omp critical (cout)
        {
        std::cout << "Warining: in CorrelationPlate::set_mean_z_in_bin(index, value): The given index is out of bouds, ignoring..." << std::endl;
        }
    }
}

void CorrelationPlate::set_num_averaged_pairs(size_t index, int value){
    /**
     EXPLANATION:
     Set function for num_average_pairs_
     
     INPUTS:
     index - index of the selected num_average_pairs_ element
     value - element's new value
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < num_averaged_pairs_.size()){
        num_averaged_pairs_[index] = value;
    }
    else{
        #pragma omp critical (cout)
        {
        std::cout << "Warining: in CorrelationPlate::set_num_average_pairs(index, value): The given index is out of bouds, ignoring..." << std::endl;
        }
    }
}

void CorrelationPlate::set_weight(size_t index, double value){
    /**
     EXPLANATION:
     Set function for weight_
     
     INPUTS:
     index - index of the selected weight_ element
     value - element's new value
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < weight_.size()){
        weight_[index] = value;
    }
    else{
        #pragma omp critical (cout)
        {
        std::cout << "Warining: in CorrelationPlate::set_weight(index, value): The given index is out of bouds, ignoring..." << std::endl;
        }
    }
}

void CorrelationPlate::set_xi(size_t index, double value){
    /**
     EXPLANATION:
     Set function for xi_
     
     INPUTS:
     index - index of the selected xi_ element
     value - element's new value
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < xi_.size()){
        xi_[index] = value;
    }
    else{
        #pragma omp critical (cout)
        {
        std::cout << "Warining: in CorrelationPlate::set_xi(index, value): The given index is out of bouds, ignoring..." << std::endl;
        }
    }
}

void CorrelationPlate::set_xi_correction(size_t index, double value){
    /**
     EXPLANATION:
     Set function for xi_correction_
     
     INPUTS:
     index - index of the selected xi_correction_ element
     value - element's new value
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    
    if (index < xi_.size()){
        xi_correction_[index] = value;
    }
    else{
        #pragma omp critical (cout)
        {
        std::cout << "Warining: in CorrelationPlate::set_xi_correction(index, value): The given index is out of bouds, ignoring..." << std::endl;
        }
    }
}

void CorrelationPlate::AddPair(const size_t& k_index, const LyaPixel& pixel, const double& pi, const double& sigma, const LyaMeanProjectedDeltasInterpolationMap& mean_proj_deltas){
    /**
     EXPLANATION:
     Adds pair contribution to xi in the specified bin
     
     INPUTS:
     k_index - an integer specifying the bin to add the contribution to
     pixel - an LyaPixel instance to add the contribution from
     pi - a double specifying the parallel separation of the pair
     sigma - a double specifying the perpendicular separation of the pair
     mean_proj_deltas - a LyaMeanProjectedDeltasInterpolationMap with the correction to the cross-correlation due to the projection
                        of the delta field. Ignored if flag_projection_correction_ is set to 'false'
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     LyaPixel
     
     FUNCITONS USED:
     NONE
     */
    
    if (k_index > xi_.size()){
        #pragma omp critical (cout)
        {
        std::cout << "In function CorrelationPlate::AddPair, the given index is out of bounds. Ignoring..." << std::endl;
        }
    }
    if (k_index > xi_correction_.size()){
        #pragma omp critical (cout)
        {
        std::cout << "In function CorrelationPlate::AddPair, something strange is going on... Ignoring..." << std::endl;
        }
    }
    
    xi_[k_index] += pixel.delta()*pixel.weight();
    mean_pi_[k_index] += pi*pixel.weight();
    mean_sigma_[k_index] += sigma*pixel.weight();
    mean_z_in_bin_[k_index] += pixel.z()*pixel.weight();
    weight_[k_index] += pixel.weight();
    num_averaged_pairs_[k_index] ++;
    if (flag_projection_correction_){
        double correction = mean_proj_deltas.LinearInterpolation(pixel.z());
        if (correction != _BAD_DATA_){
            xi_correction_[k_index] += correction*pixel.weight();
        }
    }
}

void CorrelationPlate::ComputeCrossCorrelation(const AstroObjectDataset& object_list, const SpectraDataset& spectra_list, const Input& input, const LyaMeanProjectedDeltasInterpolationMap& mean_proj_deltas){
    /**
     EXPLANATION:
     Computes the cross-correlation
     
     INPUTS:
     object_list - an AstroObjectDataset instance
     spectra_list - a SpectraDataset instance
     input - a Input instance to load bin settings
     mean_proj_deltas - a LyaMeanProjectedDeltasInterpolationMap with the correction to the cross-correlation due to the projection
                        of the delta field. Ignored if flag_projection_correction_ is set to 'false'
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObjectDataset
     CorrelationPlate
     Input
     LyaPixel
     LyaSpectrum
     SpectraDataset
     
     FUNCITONS USED:
     NONE
     */
    if (plate_number_ == _NORM_){
        #pragma omp critical (cout)
        {
            std::cout << "Warning : In CorrelationPlate::ComputeCrossCorrelation : Plate number is set to _NORM_. The cross-correlation should not be computed in this CorrelationPlate instance. Ignoring..." << std::endl;
        }
        return;
    }
    
    // load bin settings
    double max_pi = input.max_pi();
    double max_sigma = input.max_sigma();
    double step_pi = input.step_pi();
    double step_sigma = input.step_sigma();
    double num_pi_bins = input.num_pi_bins();
    double num_sigma_bins = input.num_sigma_bins();    
    
    size_t number_of_objects = object_list.num_objects_in_plate(plate_number_);
    
    if (flag_verbose_correlation_plate_ >= 1){
        #pragma omp critical (cout)
        {
            std::cout << "Computing cross-correlation in plate " << plate_number_ << ". In this plate there are " << number_of_objects << " AstroObjects" << std::endl;
        }
    }

    // loop over AstroObjects
    for (size_t i = 0; i < number_of_objects; i ++){
        
        AstroObject object = object_list.list(plate_number_, i);
        
        // checking that the object was loaded successfully
        if (object.dist() == _BAD_DATA_){
            if (flag_verbose_correlation_plate_ >= 2){
                #pragma omp critical (cout)
                {
                    std::cout << "_BAD_DATA_ AstroObject found. Ignoring..." << std::endl;
                }
            }
            continue;
        }
        
        double sigma_aux_o = 2.0*object.dist(); // auxiliar variable to compute sigma values
        
        // loop over neighbouring plates
        for (size_t j = 0; j < plate_neighbours_.size(); j ++){
            
            size_t number_of_spectra = spectra_list.num_objects_in_plate(plate_neighbours_[j]);
            
            // loop over LyaSpectra
            for (size_t k = 0; k < number_of_spectra; k ++){
                
                LyaSpectrum lya_spectrum = spectra_list.list(plate_neighbours_[j], k);
                // checking that the object was loaded successfully
                if (lya_spectrum.dist() == _BAD_DATA_){
                    if (flag_verbose_correlation_plate_ >= 2){
                        #pragma omp critical (cout)
                        {
                            std::cout << "_BAD_DATA_ LyaSpectrum found. Ignoring..." << std::endl;
                        }
                    }
                    continue;
                }
                std::vector<LyaPixel> spectrum = lya_spectrum.spectrum();
                
                // compute angular separation
                double cos_theta = object.angle().CosAngularDistance(lya_spectrum.angle());
                
                double sigma_aux;
                if (cos_theta == 1.0){
                    sigma_aux = 0.0;
                }
                else{
                    sigma_aux = sigma_aux_o*(1.0-cos_theta); // auxiliar variable to compute sigma values
                }
                
                // check if all pixels in the spectrum are too far apart
                double pair_min_sigma = sqrt(sigma_aux*spectrum[0].dist()); // minimum distance obtained for lowest redshift pixel
                if (pair_min_sigma > max_sigma){ // if the minimum value for sigma (r_perp) is too large, the whole spectra is discarded
                    if (flag_verbose_correlation_plate_ >= 2){
                        #pragma omp critical (cout)
                        {
                            std::cout << "Spectrum rejected: min_sigma is too large" << std::endl;
                        }
                        if (flag_verbose_correlation_plate_ >= 3){
                            #pragma omp critical (cout)
                            {
                                std::cout << "min_sigma = " << pair_min_sigma << std::endl;
                                std::cout << "ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                                std::cout << object.angle() << " " << lya_spectrum.angle() << " " << cos_theta << " " << acos(cos_theta) << " " << pair_min_sigma << std::endl << std::endl;
                            }
                        }
                    }
                    continue;
                } 
                if (pair_min_sigma <= 0.1){ // if the spectrum is being paired with its quasar, the whole spectra is discarded
                    if (flag_verbose_correlation_plate_ >= 2){
                        #pragma omp critical (cout)
                        {
                            std::cout << "Spectrum rejected: It is being paired with its origin quasar" << std::endl;
                        }
                    }
                    continue;
                }
                
                double pair_max_pi = spectrum.back().dist() - object.dist(); // maximum distance obtained for highest redshift pixel
                if (pair_max_pi < -max_pi) { // if maximum value for pi (r_par) is too small, the whole spectra is discarded
                    if (flag_verbose_correlation_plate_ >= 2){
                        #pragma omp critical (cout)
                        {
                            std::cout << "Spectrum rejected: max_pi is too small" << std::endl;
                        }
                        if (flag_verbose_correlation_plate_ >= 3){
                            #pragma omp critical (cout)
                            {
                                std::cout << "max_pi = " << pair_max_pi << std::endl;
                                std::cout << "ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                                std::cout << object.angle() << " " << lya_spectrum.angle() << " " << cos_theta << " " << acos(cos_theta) << " " << pair_max_pi << std::endl << std::endl;
                            }
                        }
                    }
                    continue;
                }
                
                double pair_min_pi = spectrum[0].dist() - object.dist(); // minimum distance obtained for lowest redshift pixel
                if (pair_min_pi > max_pi){ // if minimum value for pi (r_par) is too high, the whole spectrum is discarded
                    if (flag_verbose_correlation_plate_ >= 2){
                        #pragma omp critical (cout)
                        {
                            std::cout << "Spectrum rejected: min_pi is too large" << std::endl;
                        }
                        if (flag_verbose_correlation_plate_ >= 3){
                            #pragma omp critical (cout)
                            {
                                std::cout << "min_pi = " << pair_min_pi << std::endl;
                                std::cout << "ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                                std::cout << object.angle() << " " << lya_spectrum.angle() << " " << cos_theta << " " << acos(cos_theta) << " " << pair_min_pi << std::endl << std::endl;
                            }
                        }
                    }
                    continue;
                }
                
                // loop over LyaPixels
                for (int p = 0; p < spectrum.size(); p ++){
                    
                    // cehck that the weight is not zero
                    if (spectrum[p].weight() == 0.0){
                        continue;
                    }
                    
                    // compute pi and sigma
                    double sigma = sqrt(sigma_aux*spectrum[p].dist());
                    if (sigma > max_sigma){ // if sigma (r_perp) is too large, discard pixel
                        if (flag_verbose_correlation_plate_ >= 3){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Pixel rejected: sigma is too large" << std::endl;
                                std::cout << "sigma = " << sigma << std::endl;
                            }
                        }
                        continue;
                    }

                    double pi = spectrum[p].dist()-object.dist();
                    if ((pi > max_pi) or (pi <- max_pi)){ // if pi (r_par) is too large or too small, discard pixel
                        if (flag_verbose_correlation_plate_ >= 3){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Pixel rejected: abs(pi) is too large" << std::endl;
                                std::cout << "pi = " << pi << std::endl;
                            }
                        }
                        continue;
                    }
                    
                    // locate pi pixel (i_index)
                    int i_index = int(pi/step_pi)+num_pi_bins/2;
                    if (pi<0.0){
                        i_index -= 1;
                    }
                    if (i_index < 0){
                        if (flag_verbose_correlation_plate_ >= 2){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Pixel rejected: Bad indexing: i_index = " << i_index << std::endl;
                            }
                        }
                        continue;
                    }

                    
                    // locate sigma pixel (j_index)
                    int j_index = int(sigma/step_sigma);
                    if (j_index < 0){
                        if (flag_verbose_correlation_plate_ >= 2){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Pixel rejected: Bad indexing: j_index = " << j_index << std::endl;
                            }
                        }
                        continue;
                    }
                    
                    // locate xi pixel (k)
                    int k_index = i_index*num_sigma_bins+j_index;
                    if (k_index < 0){
                        if (flag_verbose_correlation_plate_ >= 2){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Pixel rejected: Bad indexing: k_index = " << k_index << std::endl;
                            }
                        }
                        continue;
                        
                    }
                    // add contribution to xi in the specified bin
                    if (flag_verbose_correlation_plate_ >= 3){
                        #pragma omp critical (cout)
                        {
                            std::cout << "Pixel accepted: assign contribution to bin " << k_index << std::endl;
                        }
                    }
                    AddPair(k_index, spectrum[p], pi, sigma, mean_proj_deltas);
                    
                    // add contribution to mean redshift
                    if (flag_verbose_correlation_plate_ >= 3){
                        #pragma omp critical (cout)
                        {
                            std::cout << "Pixel accepted: assign contribution to mean redshift " << k_index << std::endl;
                        }
                    }
                    mean_z_ += spectrum[p].z()*spectrum[p].weight();
                    weight_z_ += spectrum[p].weight();
                    
                    // write down pair information in bin file
                    if (flag_write_partial_results_ >= 2){
                        if (flag_verbose_correlation_plate_ >= 3){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Pixel accepted: keeping pair" << std::endl;
                            }
                        }
                        KeepPair(k_index, lya_spectrum, p, plate_number_, i);
                    }
                    
                }
            }
        }
    }
    if (flag_verbose_correlation_plate_ >= 3){
        #pragma omp critical (cout)
        {
        std::cout << "CorrelationPlate: correlation computed" << std::endl;
        }
    }
    
    if (flag_write_partial_results_ >= 2){
        for (size_t i = 0; i < num_bins_; i++){
            if (position_[i] > 0){
                if (flag_verbose_correlation_plate_ >= 3){
                    #pragma omp critical (cout)
                    {
                    std::cout << "CorrelationPlate: saving pairs" << std::endl;
                    }
                }
                SavePairs(i);
            }
        }
    }
    
    if (flag_verbose_correlation_plate_ >= 1){
        #pragma omp critical (cout)
        {
            std::cout << "Computed cross-correlation in plate " << plate_number_ << std::endl;
        }
    }
}

std::string CorrelationPlate::Info(size_t bin){
    /**
     EXPLANATION:
     Returns a string with the information for the selected bin
     
     INPUTS:
     bin - an unsigned integral specifying the bin
     
     OUTPUTS:
     out - a string with the information
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     ToStr
     */
    string out = "";
    
    if (flag_projection_correction_){
        out += ToStr(xi_[bin] - xi_correction_[bin]) + " " + ToStr(xi_correction_[bin]) + " ";
    }
    out += ToStr(xi_[bin]) + " " + ToStr(mean_pi_[bin]) + " " + ToStr(mean_sigma_[bin]) + " " + ToStr(mean_z_in_bin_[bin]) + " " + ToStr(weight_[bin]) + " " + ToStr(num_averaged_pairs_[bin]) + " " + ToStr(plate_number_);
    
    return out;
}

std::string CorrelationPlate::InfoHeader(const bool& flag_projection_correction){
    /**
     EXPLANATION:
     Returns a string with the columns information
     
     INPUTS:
     flag_projection_correction - a boolean specifying of the correction to the cross-correlation due to projecting deltas is being computed
     
     OUTPUTS:
     a string with the information
     
     CLASSES USED:
     NONE
     
     FUNCITONS USED:
     NONE
     */
    
    if (flag_projection_correction){
        return "xi xi_correction xi_uncorrected mean_pi mean_sigma mean_z weight num_averaged_pairs plate";
    }
    else{
        return "xi mean_pi mean_sigma mean_z weight num_averaged_pairs plate";
    }
}

void CorrelationPlate::KeepPair(const int& k_index, const LyaSpectrum& lya_spectrum, const size_t& pixel_number, const size_t obj_plate, const size_t obj_num){
    /**
     EXPLANATION:
     Keeps the pair information to save at an appropiate time
     
     INPUTS:
     k_index - an integer specifying the pair's bin
     lya_spectrum - a LyaSpectrum instance
     pixel_number - an unsigned integral with the position of the pixel contributing to the pair
     obj_plate - plate the object is found in
     obj_num - number of object in that plate list
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     CorrelationPlate
     LyaPixel
     LyaSpectrum
     SpherePoint
     
     FUNCITONS USED:
     ToStr
     */
    
    if (flag_verbose_correlation_plate_ >= 3){
        #pragma omp critical (cout)
        {
            std::cout << "CorrelationPlate: format pair to store" << std::endl;
        }
    }
    
    LyaPixel pixel = lya_spectrum.spectrum(pixel_number);
    
    Pair pair(obj_plate, obj_num, lya_spectrum.plate(), lya_spectrum.fiber(), lya_spectrum.mjd(), pixel_number, pixel.delta(), pixel.z(), pixel.weight());
    
    if (flag_verbose_correlation_plate_ >= 3){
        #pragma omp critical (cout)
        {
            std::cout << "CorrelationPlate: store pair" << std::endl;
        }
    }
    
    pairs_information_[k_index][position_[k_index]] = pair;
    position_[k_index] ++;
    
    if (position_[k_index] == max_pairs_){
        if (flag_verbose_correlation_plate_ >= 3){
            #pragma omp critical (cout)
            {
                std::cout << "CorrelationPlate: saving pairs" << std::endl;
            }
        }
        SavePairs(k_index);
    }
}

void CorrelationPlate::Normalize(){
    /**
     EXPLANATION:
     Normalizes the cross correlation results
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    
    if (plate_number_ == _NORM_){
        
        for (size_t i = 0; i < num_bins_; i ++){
                
            if (weight_[i] == 0.0){
                // if the weight is zero, complain and exit the function
                std::cerr << "Zero Division Error in CorrelationPlate::Normalize. Ignoring..." << std::endl;
                continue;
            }
            // if the weight is 1.0, normalization is not needed
            else if (weight_[i] != 1.0){
                
                #pragma opm critical (normalize)
                {
                xi_[i] /= weight_[i];
                mean_pi_[i] /= weight_[i];
                mean_sigma_[i] /= weight_[i];
                mean_z_in_bin_[i] /= weight_[i];
                if (flag_projection_correction_){
                    xi_correction_[i] /= weight_[i];
                }
                }
            }
            
        }
        mean_z_ /= weight_z_;
        
    }
    else{
        // if plate number is not _NORM_, the instance is not supposed to normalize
        #pragma omp critical (cout)
        {
            std::cout << "Warning : In CorrelationPlate::Normalize : Plate number is not set to _NORM_. This CorrelationPlate instance should not be normalized. Ignoring..." << std::endl;
        }
    }
    
}

void CorrelationPlate::SaveCrossCorrelation(const Input& input){
    /**
     EXPLANATION:
     Saves the cross correlation measured in a specific plate
     
     INPUTS:
     input - a Input instance 
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    
    std::string filename;
    
    if (plate_number_ == _NORM_){
        // if plate number is _NORM_, the instance is not supposed to be saved using this function
        #pragma omp critical (cout)
        {
            std::cout << "Warning : In CorrelationPlate::Normalize : Plate number is not set to _NORM_. This CorrelationPlate instance should not be normalized. Ignoring..." << std::endl;
        }
    }
    else{
        // save normalized cross-correlation
        if (flag_verbose_correlation_plate_ >= 1){
            #pragma omp critical (cout)
            {
                std::cout << "Saving cross-correlation for plate " << pairs_file_name_ << std::endl;
            }
        }
        if (flag_projection_correction_){
            filename = input.results() + "plate_" + pairs_file_name_ + ".data";
            {
                std::ofstream file(filename.c_str(),std::ofstream::trunc);
                if (file.is_open()){
                    for (size_t i = 0; i < num_bins_; i++){
                        
                        if (weight_[i] != 0.0){
                            file << i << " " << (xi_[i]-xi_correction_[i])/weight_[i] << std::endl;
                        }
                        else{
                            file << i << " NaN" << std::endl;
                        }
                    }
                    
                    file.close();
                }
                else{
                    #pragma omp critical (cout)
                    {
                        std::cout << "Error : In CorrelationPlate::SaveCrossCorrelation : Unable to open file:" << std::endl << filename << std::endl;
                    }
                }
            }
            filename = input.results() + "plate_" + pairs_file_name_ + ".uncorrected.data";
            {
                std::ofstream file(filename.c_str(),std::ofstream::trunc);
                if (file.is_open()){
                    for (size_t i = 0; i < num_bins_; i++){
                        
                        if (weight_[i] != 0.0){
                            file << i << " " << xi_[i]/weight_[i] << std::endl;
                        }
                        else{
                            file << i << " NaN" << std::endl;
                        }
                    }
                    
                    file.close();
                }
                else{
                    #pragma omp critical (cout)
                    {
                        std::cout << "Error : In CorrelationPlate::SaveCrossCorrelation : Unable to open file:" << std::endl << filename << std::endl;
                    }
                }
            }
            filename = input.results() + "plate_" + pairs_file_name_ + ".correction.data";
            {
                std::ofstream file(filename.c_str(),std::ofstream::trunc);
                if (file.is_open()){
                    for (size_t i = 0; i < num_bins_; i++){
                        
                        if (weight_[i] != 0.0){
                            file << i << " " << xi_correction_[i]/weight_[i] << std::endl;
                        }
                        else{
                            file << i << " NaN" << std::endl;
                        }
                    }
                    
                    file.close();
                }
                else{
                    #pragma omp critical (cout)
                    {
                        std::cout << "Error : In CorrelationPlate::SaveCrossCorrelation : Unable to open file:" << std::endl << filename << std::endl;
                    }
                }
            }
        }
        else{
            filename = input.results() + "plate_" + pairs_file_name_ + ".data";
            {
                std::ofstream file(filename.c_str(),std::ofstream::trunc);
                if (file.is_open()){
                    for (size_t i = 0; i < num_bins_; i++){
                        
                        if (weight_[i] != 0.0){
                            file << i << " " << xi_[i]/weight_[i] << std::endl;
                        }
                        else{
                            file << i << " NaN" << std::endl;
                        }
                    }
                    
                    file.close();
                }
                else{
                    #pragma omp critical (cout)
                    {
                        std::cout << "Error : In CorrelationPlate::SaveCrossCorrelation : Unable to open file:" << std::endl << filename << std::endl;
                    }
                }
            }
        }
        // save normalized grid
        filename = input.results() + "plate_" + pairs_file_name_ + ".grid";
        {
            std::ofstream file(filename.c_str(),std::ofstream::trunc);
            if (file.is_open()){
                
                for (size_t i = 0; i < num_bins_; i++){
                    
                    if (weight_[i] != 0.0){
                        file << i << " " <<  mean_pi_[i]/weight_[i] << " " << mean_sigma_[i]/weight_[i] << " " << mean_z_in_bin_[i]/weight_[i] << std::endl;
                    }
                    else{
                        file << i << " NaN NaN Nan" << std::endl;
                    }
                }
                
                file.close();
            }
            else{
                #pragma omp critical (cout)
                {
                    std::cout << "Error : In CorrelationResults::SaveCrossCorrelation : Unable to open file:" << std::endl << filename << std::endl;
                }
            }
        }

    }
    
}

void CorrelationPlate::SavePairs(const int& k_index){
    /**
     EXPLANATION:
     Writes down pair information in bin file
     
     INPUTS:
     k_index - an integer specifying the pair's bin
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     CorrelationPlate
     LyaPixel
     LyaSpectrum
     SpherePoint
     
     FUNCITONS USED:
     ToStr
     */
    
    std::string filename;
    
    filename = results_ + ToStr(k_index) + ".dat";
    std::ofstream bin_file(filename.c_str(),std::ofstream::app);
    
    if (bin_file.is_open()){
        for (size_t i = 0; i < position_[k_index]; i++){
            bin_file << pairs_information_[k_index][i].obj_plate() << " ";
            bin_file << pairs_information_[k_index][i].obj_num() << " ";
            bin_file << pairs_information_[k_index][i].spec_plate() << " ";
            bin_file << pairs_information_[k_index][i].spec_fiber() << " ";
            bin_file << pairs_information_[k_index][i].spec_mjd() << " ";
            bin_file << pairs_information_[k_index][i].pixel_number() << " ";
            bin_file << pairs_information_[k_index][i].pixel_delta() << " ";
            bin_file << pairs_information_[k_index][i].pixel_z() << " ";
            bin_file << pairs_information_[k_index][i].pixel_weight() << " \n";
            
        }
        bin_file.close();
    }
    else{
        #pragma omp critical (cout)
        {
            std::cout << "Error : In CorrelationPlate::SavePairs : Unable to open file:" << std::endl << filename << std::endl;
        }
    }
    
    CorrelationPlate::position_[k_index] = 0;
}

void CorrelationPlate::operator+= (const CorrelationPlate& other){
    /**
     EXPLANATION:
     Overloads the += operator. Updates the xi_, mean_sigma, mean_pi and weight_ or cov_mat_ values by adding those in other.
     
     INPUTS:
     other - a CorrelationPlate instance to be added
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    
    // check that both instances have the same number of bins
    if (num_bins_ != other.num_bins()){
        #pragma omp critical (cout)
        {
            std::cout << "Warning : In CorrelationPlate::operator+= : Trying to add CorrelationPlates with different number of bins. Ignoring..." << std::endl;
        }
        return;
    }
    // check that both instances have the same value for flag_projection_correction
    if (flag_projection_correction_ != other.flag_projection_correction()){
        #pragma omp critical (cout)
        {
            std::cout << "Warning : In CorrelationPlate::operator+= : Trying to add CorrelationPlates with different flag_projection_correction. Ignoring..." << std::endl;
        }
        return;
    }
    
    for (size_t i = 0; i < xi_.size(); i ++){
        xi_[i] += other.xi(i);
        mean_pi_[i] += other.mean_pi(i);
        mean_sigma_[i] += other.mean_sigma(i);
        mean_z_in_bin_[i] += other.mean_z_in_bin(i);
        weight_[i] += other.weight(i);
        num_averaged_pairs_[i] += other.num_averaged_pairs(i);
        if (flag_projection_correction_){
            xi_correction_[i] += other.xi_correction(i);
        }
    }
    mean_z_ += other.mean_z();
    weight_z_ += other.weight_z();
    
}

CorrelationPlate CorrelationPlate::operator- (const CorrelationPlate& other){
    /**
     EXPLANATION:
     Overloads the - operator. Returns a new instance with the xi_, mean_sigma, mean_pi and weight_ values or cov_mat subtracted by those in other.
     
     INPUTS:
     other - a CorrelationPlate instance to be subtracted
     
     OUTPUTS:
     temp - a CorrelationPlate instance
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */

    CorrelationPlate temp;
    temp = CorrelationPlate(plate_number_, num_bins_, results_, pairs_file_name_, plate_neighbours_, flag_verbose_correlation_plate_, flag_write_partial_results_, flag_projection_correction_);
    
    // check that both instances have the same number of bins
    if (num_bins_ != other.num_bins()){
        #pragma omp critical (cout)
        {
            std::cout << "Warning : In CorrelationPlate::operator- : Trying to add CorrelationPlates with different number of bins. Returning zero filled CorrelationPlates..." << std::endl;
        }
        return temp;
    }
    // check that both instances have the same value for flag_projection_correction
    if (flag_projection_correction_ != other.flag_projection_correction()){
        #pragma omp critical (cout)
        {
            std::cout << "Warning : In CorrelationPlate::operator- : Trying to add CorrelationPlates with different flag_projection_correction. Ignoring..." << std::endl;
        }
        return temp;
    }
    
    for (size_t i = 0; i < xi_.size(); i ++){
        temp.set_xi(i, xi_[i] - other.xi(i));
        temp.set_mean_pi(i, mean_pi_[i] - other.mean_pi(i));
        temp.set_mean_sigma(i, mean_sigma_[i] - other.mean_sigma(i));
        temp.set_mean_z_in_bin(i, mean_z_in_bin_[i] - other.mean_z_in_bin(i));
        temp.set_weight(i, weight_[i] - other.weight(i));
        temp.set_num_averaged_pairs(i, num_averaged_pairs_[i] - other.num_averaged_pairs(i));
        if (flag_projection_correction_){
            temp.set_xi_correction(i, xi_correction_[i] - other.xi_correction(i));
        }
    }
    temp.set_mean_z(mean_z_ - other.mean_z());
    temp.set_weight_z(weight_z_ - other.weight_z());
    
    return temp;
    
}
        
CorrelationPlate CorrelationPlate::operator* (const CorrelationPlate& other){
    /**
     EXPLANATION:
     Overloads the * operator. Returns a new instance with the xi_, mean_sigma, mean_pi and weight_ values multiplied by those in other.
     
     INPUTS:
     other - a CorrelationPlate instance to be subtracted
     
     OUTPUTS:
     temp - a CorrelationPlate instance
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    CorrelationPlate temp;
    temp = CorrelationPlate(plate_number_, num_bins_, results_, pairs_file_name_, plate_neighbours_, flag_verbose_correlation_plate_, flag_write_partial_results_, flag_projection_correction_);
    
    // check that both instances have the same number of bins
    if (num_bins_ != other.num_bins()){
        #pragma omp critical (cout)
        {
            std::cout << "Warning : In CorrelationPlate::operator* : Trying to add CorrelationPlates with different number of bins. Returning zero filled CorrelationPlates..." << std::endl;
        }
        return temp;
    }
    // check that both instances have the same value for flag_projection_correction
    if (flag_projection_correction_ != other.flag_projection_correction()){
        #pragma omp critical (cout)
        {
            std::cout << "Warning : In CorrelationPlate::operator* : Trying to add CorrelationPlates with different flag_projection_correction. Ignoring..." << std::endl;
        }
        return temp;
    }
    
    for (size_t i = 0; i < xi_.size(); i ++){
        temp.set_xi(i, xi_[i]*other.xi(i));
        temp.set_mean_pi(i, mean_pi_[i]*other.mean_pi(i));
        temp.set_mean_sigma(i, mean_sigma_[i]*other.mean_sigma(i));
        temp.set_mean_z_in_bin(i, mean_z_in_bin_[i]*other.mean_z_in_bin(i));
        temp.set_weight(i, weight_[i]*other.weight(i));
        temp.set_num_averaged_pairs(i, num_averaged_pairs_[i]*other.num_averaged_pairs(i));
        if (flag_projection_correction_){
            temp.set_xi_correction(i, xi_correction_[i]*other.xi_correction(i));
        }
    }
    temp.set_mean_z(mean_z_*other.mean_z());
    temp.set_weight_z(weight_z_*other.weight_z());

            
    return temp;
    
}

