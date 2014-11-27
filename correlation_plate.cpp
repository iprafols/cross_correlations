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
        std::cout << "Error while initializing a CorrelationPlate 'bad data' instance" << std::endl;
        std::exit;
    }
    
    plate_number_ = _BAD_DATA_INT_;
    
}

CorrelationPlate::CorrelationPlate(const int plate_number, const size_t num_bins, const std::string& results, const std::string& pairs_file_name, const std::vector<int>& plate_neighbours, const size_t flag_write_partial_results){
    /**
     EXPLANATION:
     Cosntructs a CorrelationPlate instance and initializes all its variables
     
     INPUTS:
     plate_number - an integer with the plate number
     num_bins - an unsigned integer with the number of bins
     results - a string with the path to the results directory (missing the bin number)
     pairs_file_name - a string with the base name of the file storing the pairs information
     plate_neighbours - a vector containing the plate numbers of the neighbouring plates
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     ToStr
     */
    flag_write_partial_results_ = flag_write_partial_results;
    
    plate_number_ = plate_number;
    plate_neighbours_ = plate_neighbours;
    num_bins_ = num_bins;
    if (plate_number_ == _NORM_){
        results_ = "";
    }
    else{
        results_ = results;
    }
    pairs_file_name_ = pairs_file_name;
    
    xi_.resize(num_bins_,0.0);
    mean_pi_.resize(num_bins_,0.0);
    mean_sigma_.resize(num_bins_,0.0);
    weight_.resize(num_bins_,0.0);
    num_averaged_pairs_.resize(num_bins_,0);

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
        std::cout << "the given index is out of bouds, ignoring..." << std::endl;
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
        std::cout << "the given index is out of bouds, ignoring..." << std::endl;
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
        std::cout << "the given index is out of bouds, ignoring..." << std::endl;
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
        std::cout << "the given index is out of bouds, ignoring..." << std::endl;
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
        std::cout << "the given index is out of bouds, ignoring..." << std::endl;
    }
}

void CorrelationPlate::AddPair(const int& k_index, const LyaPixel& pixel, const double& pi, const double& sigma){
    /**
     EXPLANATION:
     Adds pair contribution to xi in the specified bin
     
     INPUTS:
     k_index - an integer specifying the bin to add the contribution to
     pixel - an LyaPixel instance to add the contribution from
     pi - a double specifying the parallel separation of the pair
     sigma - a double specifying the perpendicular separation of the pair
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     LyaPixel
     
     FUNCITONS USED:
     NONE
     */
    
    xi_[k_index] += (pixel.forest()-1.0)*pixel.weight();
    mean_pi_[k_index] += pi*pixel.weight();
    mean_sigma_[k_index] += sigma*pixel.weight();
    weight_[k_index] += pixel.weight();
    num_averaged_pairs_[k_index] ++;
    
}

void CorrelationPlate::ComputeCrossCorrelation(const AstroObjectDataset& object_list, const LyaSpectraDataset& spectra_list, const Input& input){
    /**
     EXPLANATION:
     Computes the cross-correlation
     
     INPUTS:
     object_list - an AstroObjectDataset instance
     spectra_list - a LyaSpectraDataset instance
     input - a Input instance to load bin settings
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObjectDataset
     CorrelationPlate
     Input
     LyaPixel
     LyaSpectraDataset
     LyaSpectrum
     
     FUNCITONS USED:
     NONE
     */
    if (plate_number_ == _NORM_){
        std::cout << "Warning: Plate number is set to _NORM_. The cross-correlation should not be computed in this CorrelationPlate instance. Ignoring..." << std::endl;
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
    
    std::cout << "computing cross-correlation in plate " << plate_number_ << "; in this plate there are " << number_of_objects << " AstroObjects" << std::endl;

    // loop over AstroObjects
    for (size_t i = 0; i < number_of_objects; i ++){
        
        AstroObject object = object_list.list(plate_number_, i);
        
        // checking that the object was loaded successfully
        if (object.dist() == _BAD_DATA_){
            std::cout << "TEST: _BAD_DATA_ AstroObject found" << std::endl;
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
                    std::cout << "TEST: _BAD_DATA_ LyaSpectrum found" << std::endl;
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
                    //std::cout << "TEST: spectra rejected because min_sigma is too large. value = " << pair_min_sigma << std::endl;
                    //std::cout << "TEST: ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                    //std::cout << "TEST: " << object.angle() << " " << lya_spectrum.angle() << " " << cos_theta << " " << acos(cos_theta) << " " << pair_min_sigma << std::endl;
                    continue;
                } 
                if (pair_min_sigma <= 0.1){ // if the spectrym is being paired with its quasar, the whole spectra is discarded
                    continue;
                }
                
                double pair_max_pi = object.dist()-spectrum[0].dist(); // minimum distance obtained for lowest redshift pixel
                if (pair_max_pi < -max_pi) { // if maximum value for pi (r_par) is too small, the whole spectra is discarded
                    //std::cout << "TEST: spectra rejected because min_pi is too small" << std::endl;
                    continue;
                }
                
                double pair_min_pi = object.dist()-spectrum.back().dist(); // maximum distance obtained for highest redshift pixel
                if (pair_min_pi > max_pi){ // if minimum value for pi (r_par) is too high, the whole spectrum is discarded
                    //std::cout << "TEST: spectra rejected because min_pi is too large" << std::endl;
                    continue;
                }
                
                // loop over LyaPixels
                for (int p = 0; p < spectrum.size(); p ++){
                    
                    // compute pi and sigma
                    double sigma = sqrt(sigma_aux*spectrum[p].dist());
                    if (sigma > max_sigma){ // if sigma (r_perp) is too large, discard pixel
                        continue;
                    }

                    double pi = object.dist()-spectrum[p].dist();
                    if ((pi > max_pi) or (pi <- max_pi)){ // if pi (r_par) is too large or too small, discard pixel
                        continue;
                    }
                    
                    //std::cout << "TEST: accepted pixel" << std::endl;
                    // locate pi pixel (i_index)
                    int i_index = int(pi/step_pi)+num_pi_bins/2;
                    if (pi<0.0){
                        i_index -= 1;
                    }
                    if (i_index < 0){
                        //std::cout << "TEST: pi = " << pi << "; step_pi = " << step_pi << "; pi/step_pi = " << pi/step_pi << "; int(pi/step_pi) = " << int(pi/step_pi) << std::endl;
                    }

                    
                    // locate sigma pixel (j_index)
                    int j_index = int(sigma/step_sigma);
                    if (j_index < 0){
                        //std::cout << "TEST: sigma = " << sigma << "; step_sigma = " << step_sigma << "; sigma/step_sigma = " << sigma/step_sigma << "; int(sigma/step_sigma) = " << int(sigma/step_sigma) << "; sqrt(sigma_aux*spectrum[p].dist()) = " << sqrt(sigma_aux*spectrum[p].dist()) << "; sigma_aux*spectrum[p].dist() = " << sigma_aux*spectrum[p].dist() << "; sigma_aux = " << sigma_aux << "; spectrum[p].dist() = " << spectrum[p].dist() << "; sigma_aux_o = " << sigma_aux_o << "; cos_theta = " << cos_theta <<                        std::endl;
                    }
                    
                    // locate xi pixel (k)
                    int k_index = i_index*num_sigma_bins+j_index;
                    if (k_index < 0 or i_index < 0 or j_index < 0){
                        //std::cout << "TEST: i_index = " << i_index << " j_index = " << j_index << " k_index = " << k_index << std::endl;
                        continue;
                        
                    }
                    // add contribution to xi in the specified bin
                    AddPair(k_index, spectrum[p], pi, sigma);
                    
                    // write down pair information in bin file
                    if (flag_write_partial_results_ >= 2){
                        SavePair(k_index, object, lya_spectrum, p, pi, sigma);
                    }
                    
                }
            }
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
    
    out += ToStr(xi_[bin]) + " " + ToStr(mean_pi_[bin]) + " " + ToStr(mean_sigma_[bin]) + " " + ToStr(weight_[bin]) + " " + ToStr(num_averaged_pairs_[bin]) + " " + ToStr(plate_number_);
    
    return out;
}

std::string CorrelationPlate::InfoHeader(){
    /**
     EXPLANATION:
     Returns a string with the columns information
     
     INPUTS:
     NONE
     
     OUTPUTS:
     a string with the information
     
     CLASSES USED:
     NONE
     
     FUNCITONS USED:
     NONE
     */
    
    return "xi mean_pi mean_sigma weight num_averaged_pairs plate";
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
                std::cout << "Zero Division Error in CorrelationPlate::Normalize, aborting..." << std::endl;
                return;
            }
            // if the weight is 1.0, normalization is not needed
            else if (weight_[i] != 1.0){
                
                xi_[i] /= weight_[i];
                mean_pi_[i] /= weight_[i];
                mean_sigma_[i] /= weight_[i];
                weight_[i] = 1.0;
            }
            
        }
        
    }
    else{
        // if plate number is not _NORM_, the instance is not supposed to normalize
        std::cout << "Warning: Plate number is not set to _NORM_. This CorrelationPlate instance should not be normalized. Ignoring..." << std::endl;
    }
    
}

void CorrelationPlate::SavePair(const int& k_index, const AstroObject& object, const LyaSpectrum& lya_spectrum, const size_t& p, const double& pi, const double& sigma){
    /**
     EXPLANATION:
     Writes down pair information in bin file
     
     INPUTS:
     k_index - an integer specifying the pair's bin
     object - an AstroObject instance
     lya_spectrum - a LyaSpectrum instance
     p - an unsigned integral with the position of the pixel contributing to the pair
     pi - a double specifying the parallel separation of the pair
     sigma - a double specifying the perpendicular separation of the pair
     
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
    filename = results_ + ToStr(k_index) + pairs_file_name_;
    
    std::ofstream bin_file;
    bin_file.open(filename.c_str(),std::ofstream::app); // opens the file to append content
    if (bin_file.is_open()){
        LyaPixel pixel = lya_spectrum.spectrum(p);
        // checking that the object was loaded successfully
        if (pixel.z() == _BAD_DATA_){
            std::cout << "TEST: _BAD_DATA_ LyaPixel found" << std::endl;
            return;
        }
        
        bin_file << object.angle() << " " << object.z() << " " << lya_spectrum.angle() << " " << pixel.z() << " " << pixel.dist() << " " << p << " " << pixel.forest()-1.0 << " " << pixel.weight() << " " << pi << " " << sigma << std::endl;
        bin_file.close();
    }
    else{
        std::cout << "Unable to open file:" << std::endl << filename << std::endl;
    }

}

void CorrelationPlate::operator+= (const CorrelationPlate& other){
    /**
     EXPLANATION:
     Overloads the += operator. Updates the xi_, mean_sigma, mean_pi and weight_ values by adding those in other.
     
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
    if (num_bins_ == other.num_bins()){
        for (size_t i = 0; i < num_bins_; i ++){
            xi_[i] += other.xi(i);
            mean_pi_[i] += other.mean_pi(i);
            mean_sigma_[i] += other.mean_sigma(i);
            weight_[i] += other.weight(i);
            num_averaged_pairs_[i] += other.num_averaged_pairs(i);
        }
        
    }
    // if they don't, complain
    else{
        std::cout << "Warning: Trying to add CorrelationPlates with different number of bins. Ignoring..." << std::endl;
    }
    
    
}

CorrelationPlate CorrelationPlate::operator- (const CorrelationPlate& other){
    /**
     EXPLANATION:
     Overloads the - operator. Returns a new instance with the xi_, mean_sigma, mean_pi and weight_ values subtracted by those in other.
     
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
    temp = CorrelationPlate(plate_number_, num_bins_, results_, pairs_file_name_, plate_neighbours_, flag_write_partial_results_);
    
    // check that both instances have the same number of bins
    if (num_bins_ == other.num_bins()){
        
        for (size_t i = 0; i < num_bins_; i ++){
            temp.set_xi(i, xi_[i] - other.xi(i));
            temp.set_mean_pi(i, mean_pi_[i] - other.mean_pi(i));
            temp.set_mean_sigma(i, mean_sigma_[i] - other.mean_sigma(i));
            temp.set_weight(i, weight_[i] - other.weight(i));
            temp.set_num_averaged_pairs(i, num_averaged_pairs_[i] - other.num_averaged_pairs(i));
        }
        
    }
    
    // if they don't, complain
    else{
        std::cout << "Warning: Trying to subtract CorrelationPlates with different number of bins. Returning zero filled CorrelationPlates..." << std::endl;
    }
    
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
    temp = CorrelationPlate(plate_number_, num_bins_, results_, pairs_file_name_, plate_neighbours_, flag_write_partial_results_);
    
    // check that both instances have the same number of bins
    if (num_bins_ == other.num_bins()){
        
        for (size_t i = 0; i < num_bins_; i ++){
            temp.set_xi(i, xi_[i]*other.xi(i));
            temp.set_mean_pi(i, mean_pi_[i]*other.mean_pi(i));
            temp.set_mean_sigma(i, mean_sigma_[i]*other.mean_sigma(i));
            temp.set_weight(i, weight_[i]*other.weight(i));
            temp.set_num_averaged_pairs(i, num_averaged_pairs_[i]*other.num_averaged_pairs(i));
        }
        
    }
    // if they don't, complain
    else{
        std::cout << "Warning: Trying to multiply CorrelationPlates with different number of bins. Returning zero filled CorrelationPlates..." << std::endl;
    }
            
    return temp;
    
}
   
