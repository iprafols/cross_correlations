/**
 distortion_plate.cpp
 Purpose: This files contains the body for the functions defined in distortion_plate.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 04/04/2016
 */

#include "distortion_plate.h"

DistortionPlate::DistortionPlate(int bad_data){
    /**
     EXPLANATION:
     Constructs a DistortionPlate instance and initializes all its variables
     
     INPUTS:
     bad_data - a double valued _BAD_DATA_INT_
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     DistortionPlate
     
     FUNCITONS USED:
     ToStr
     */
    if (bad_data != _BAD_DATA_INT_){
        std::cout << "Error while initializing a DistortionPlate 'bad data' instance" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    
    plate_number_ = _BAD_DATA_INT_;
    
}

DistortionPlate::DistortionPlate(const Input& input, const int plate_number, const std::vector<int>& plate_neighbours){
    /**
     EXPLANATION:
     Cosntructs a DistortionPlate instance and initializes all its variables
     
     INPUTS:
     input - a Input instance
     plate_number - an integer with the plate number
     plate_neighbours - a vector containing the plate numbers of the neighbouring plates
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     DistortionPlate
     
     FUNCITONS USED:
     ToStr
     */
    
    // set flags from input
    flag_verbose_distortion_plate_ = input.flag_verbose_distortion_plate();
    
    plate_number_ = plate_number;
    plate_neighbours_ = plate_neighbours;
    num_bins_ = input.num_bins();
    
    // initialize distortion matrix computation
    for (size_t i = 0; i < num_bins_; i++){
        for (size_t j = 0; j < num_bins_; j++){
            dist_mat_[std::pair<size_t,size_t>(i,j)] = 0.0;
        }
    }
    weight_.resize(num_bins_,0.0);
    
}

DistortionPlate::DistortionPlate(const int plate_number, const int num_bins, const std::vector<int>& plate_neighbours, size_t flag_verbose_distortion_plate){
    /**
     EXPLANATION:
     Cosntructs a DistortionPlate instance and initializes all its variables
     
     INPUTS:
     input - a Input instance
     plate_number - an integer with the plate number
     plate_neighbours - a vector containing the plate numbers of the neighbouring plates
     flag_verbose_distortion_plate - correlation_plate verbose flag
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     DistortionPlate
     
     FUNCITONS USED:
     ToStr
     */
    
    flag_verbose_distortion_plate_ = flag_verbose_distortion_plate;
    
    plate_number_ = plate_number;
    plate_neighbours_ = plate_neighbours;
    num_bins_ = num_bins;
    
    
    // initialize distortion matrix computation
    for (size_t i = 0; i < num_bins_; i++){
        for (size_t j = 0; j < num_bins_; j++){
            dist_mat_[std::pair<size_t,size_t>(i,j)] = 0.0;
        }
    }
    weight_.resize(num_bins_,0.0);
}

double DistortionPlate::weight(size_t index) const {
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

void DistortionPlate::set_dist_mat(size_t i, size_t j, double value){
    /**
     EXPLANATION:
     Set function for dist_mat_
     
     INPUTS:
     i,j - indexs of the selected dist_mat_ element
     value - element's new value
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     DistortionPlate
     
     FUNCITONS USED:
     NONE
     */
    
    CovMat::iterator it = dist_mat_.find(std::pair<size_t, size_t>(i,j));
    if (it != dist_mat_.end()){
        (*it).second = value;
    }
    else{
        std::cout << "Warining: in DistortionPlate::set_dist_mat(i, j, value): The given index is out of bouds, ignoring..." << std::endl;
    }
}

void DistortionPlate::set_weight(size_t index, double value){
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
        std::cout << "Warining: in CorrelationPlate::set_weight(index, value): The given index is out of bouds, ignoring..." << std::endl;
    }
}

void DistortionPlate::AddPair(const LyaPixel& pixel, const LyaPixel& pixel2, const size_t& i, const size_t& j, const double& forest_total_weight, const double& forest_mean_loglam, const double& forest_aux){
    /**
     EXPLANATION:
     Adds pair contribution to the distortion matrix in the specified bin
     
     INPUTS:
     pixel,pixel2 - LyaPixel instances to add the contribution from
     i,j - distortion matrix element to add the contributions to
     forest_total_weight - sum of the weights of all the pixels in the forest
     forest_mean_loglam - mean logarithm of wavelength of the forest
     forest_aux - sum of (loglam-forest_mean_loglam)*weigth of all the pixels in the forest
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     DistortionPlate
     LyaAutoInterpolationMap
     LyaPixel
     
     FUNCITONS USED:
     NONE
     */
    
    CovMat::iterator it;
    it = dist_mat_.find(std::pair<int, int>(i, j));
    
    if (it == dist_mat_.end()){
        std::cout << "Warning : In CovarianceMatrix::ComputeDistMat : Element " << i << ", " << j << " of the distortion matrix not found. Ignoring..."  << std::endl;
        return;
    }
    
    if (forest_total_weight == 0.0 or forest_aux == 0.0){
        std::cout << "Warning : In CovarianceMatrix::ComputeDistMat : Zero division error encounterd. Ignoring..."  << std::endl;
        return;
    }
    
    // add to distortion matrix
    double weight = pixel.weight();
    double add;
    if (i == j){
        add = 1.0-pixel2.weight()/forest_total_weight-(pixel2.loglam()-forest_mean_loglam)*(pixel1.loglam()-forest_mean_loglam)*pixel2.weight()/forest_aux;
    }
    else{
        add = -pixel2.weight()/forest_total_weight-(pixel2.loglam()-forest_mean_loglam)*(pixel1.loglam()-forest_mean_loglam)*pixel2.weight()/forest_aux;
    }
    
    (*it).second += add*weight;
    
}

void DistortionPlate::ComputeDistMat(const AstroObjectDataset& object_list, const SpectraDataset& spectra_list, const Input& input){
    /**
     EXPLANATION:
     Computes the distortion matrix
     
     INPUTS:
     object_list - an AstroObjectDataset instance
     spectra_list - a SpectraDataset instance
     input - a Input instance to load bin settings
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObjectDataset
     DistortionPlate
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
            std::cout << "Warning : In DistortionPlate::ComputeDistMatrix : Plate number is set to _NORM_. The distortion should not be computed in this DistortionPlate instance. Ignoring..." << std::endl;
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
    size_t pixels_separation = input.pixels_separation();
    
    size_t number_of_spectra = spectra_list.num_objects_in_plate(plate_number_);
    
    if (flag_verbose_distortion_plate_ >= 1){
        #pragma omp critical (cout)
        {
            std::cout << "Computing distortion_matrix in plate " << plate_number_ << ". In this plate there are " << number_of_spectra << " LyaSpectra" << std::endl;
        }
    }
    if (number_of_spectra == 0){
        return;
    }
    
    // loop over LyaSpectra
    for (size_t lya_spectrum_num = 0; lya_spectrum_num < number_of_spectra; lya_spectrum_num ++){
        
        LyaSpectrum lya_spectrum = spectra_list.list(plate_number_, lya_spectrum_num);
        // checking that the object was loaded successfully
        if (lya_spectrum.dist() == _BAD_DATA_){
            if (flag_verbose_distortion_plate_ >= 2){
                #pragma omp critical (cout)
                {
                    std::cout << "_BAD_DATA_ LyaSpectrum found. Ignoring..." << std::endl;
                }
            }
            continue;
        }
        std::vector<LyaPixel> spectrum = lya_spectrum.spectrum();
        
        // compute forest variables (loop over forest pixels, twice)
        double forest_aux = 0.0; //forest_aux: sum of (loglam-forest_mean_loglam)*weigth of all the pixels in the forest
        double forest_total_weight = 0.0;
        double forest_mean_loglam = 0.0;
        
        for (int pixel = 0; pixel < spectrum.size(); pixel ++){
            forest_total_weight += spectrum[pixel].weight();
            forest_mean_loglam += spectrum[pixel].loglam()*spectrum[pixel].weight();
        }
        forest_mean_loglam /= forest_total_weight;
        
        for (int pixel = 0; pixel < spectrum.size(); pixel ++){
            forest_aux += (spectrum[pixel].loglam()-forest_mean_loglam)*(spectrum[pixel].loglam()-forest_mean_loglam)*spectrum[pixel].weight();
        }
        
        // loop over neighbouring plates
        for (size_t plate_neighbours_num = 0; plate_neighbours_num < plate_neighbours_.size(); plate_neighbours_num ++){
            
            size_t number_of_objects = object_list.num_objects_in_plate(plate_neighbours_[plate_neighbours_num]);
            if (number_of_objects == 0){
                continue;
            }
            
            // loop over AstroObjects
            for (size_t object_num = 0; object_num < number_of_objects; object_num ++){
                
                AstroObject object = object_list.list(plate_neighbours_[plate_neighbours_num], object_num);
                
                // checking that the object was loaded successfully
                if (object.dist() == _BAD_DATA_){
                    if (flag_verbose_distortion_plate_ >= 2){
                        #pragma omp critical (cout)
                        {
                            std::cout << "_BAD_DATA_ AstroObject found. Ignoring..." << std::endl;
                        }
                    }
                    continue;
                }
                
                double sigma_aux_o = 2.0*object.dist(); // auxiliar variable to compute sigma values
                
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
                    if (flag_verbose_distortion_plate_ >= 2){
                        #pragma omp critical (cout)
                        {
                            std::cout << "Spectrum rejected: min_sigma is too large" << std::endl;
                        }
                        if (flag_verbose_distortion_plate_ >= 3){
                            #pragma omp critical (cout)
                            {
                                std::cout << "min_sigma = " << pair_min_sigma << std::endl;
                                std::cout << "ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                                std::cout << object.angle() << " " << lya_spectrum.angle() << " " << cos_theta << " " << acos(cos_theta) << " " << pair_min_sigma << std::endl;
                            }
                        }
                    }
                    continue;
                }
                // if the spectrum is being paired with its quasar, the whole spectra is discarded
                if (pair_min_sigma <= 0.1){
                    if (flag_verbose_distortion_plate_ >= 2){
                        #pragma omp critical (cout)
                        {
                            std::cout << "Spectrum rejected: It is being paired with its origin quasar" << std::endl;
                        }
                    }
                    continue;
                }
                
                double pair_max_pi = spectrum.back().dist() - object.dist(); // minimum distance obtained for lowest redshift pixel
                if (pair_max_pi < -max_pi) { // if maximum value for pi (r_par) is too small, the whole spectra is discarded
                    if (flag_verbose_distortion_plate_ >= 2){
                        #pragma omp critical (cout)
                        {
                            std::cout << "Spectrum rejected: max_pi is too small" << std::endl;
                        }
                        if (flag_verbose_distortion_plate_ >= 3){
                            #pragma omp critical (cout)
                            {
                                std::cout << "max_pi = " << pair_max_pi << std::endl;
                                std::cout << "ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                                std::cout << object.angle() << " " << lya_spectrum.angle() << " " << cos_theta << " " << acos(cos_theta) << " " << pair_max_pi << std::endl;
                            }
                        }
                    }
                    continue;
                }
                
                double pair_min_pi = spectrum[0].dist() - object.dist(); // maximum distance obtained for highest redshift pixel
                if (pair_min_pi > max_pi){ // if minimum value for pi (r_par) is too high, the whole spectrum is discarded
                    if (flag_verbose_distortion_plate_ >= 2){
                        #pragma omp critical (cout)
                        {
                            std::cout << "Spectrum rejected: min_pi is too large" << std::endl;
                        }
                        if (flag_verbose_distortion_plate_ >= 3){
                            #pragma omp critical (cout)
                            {
                                std::cout << "min_pi = " << pair_min_pi << std::endl;
                                std::cout << "ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                                std::cout << object.angle() << " " << lya_spectrum.angle() << " " << cos_theta << " " << acos(cos_theta) << " " << pair_min_pi << std::endl;
                            }
                        }
                    }
                    continue;
                }
                
                // loop over LyaPixels 1
                for (int pixel1 = 0; pixel1 < spectrum.size(); pixel1 ++){
                    
                    // cehck that the weight is not zero
                    if (spectrum[pixel1].weight() == 0.0){
                        continue;
                    }
                    
                    // compute pi and sigma
                    double sigma1 = sqrt(sigma_aux*spectrum[pixel1].dist());
                    if (sigma1 > max_sigma){ // if sigma (r_perp) is too large, discard pixel
                        if (flag_verbose_distortion_plate_ >= 3){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Pixel rejected: sigma is too large" << std::endl;
                                std::cout << "sigma = " << sigma1 << std::endl;
                            }
                        }
                        continue;
                    }
                    
                    double pi1 = spectrum[pixel1].dist()-object.dist();
                    if ((pi1 > max_pi) or (pi1 <- max_pi)){ // if pi (r_par) is too large or too small, discard pixel
                        if (flag_verbose_distortion_plate_ >= 3){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Pixel rejected: abs(pi) is too large" << std::endl;
                                std::cout << "pi = " << pi1 << std::endl;
                            }
                        }
                        continue;
                    }
                    
                    // locate pi pixel (i_index)
                    int i_index1 = int(pi1/step_pi)+num_pi_bins/2;
                    if (pi1<0.0){
                        i_index1 -= 1;
                    }
                    if (i_index1 < 0){
                        if (flag_verbose_distortion_plate_ >= 2){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Pixel rejected: Bad indexing: i_index = " << i_index1 << std::endl;
                            }
                        }
                        continue;
                    }
                    
                    // locate sigma pixel (j_index)
                    int j_index1 = int(sigma1/step_sigma);
                    if (j_index1 < 0){
                        if (flag_verbose_distortion_plate_ >= 2){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Pixel rejected: Bad indexing: j_index = " << j_index1 << std::endl;
                            }
                        }
                        continue;
                    }
                    
                    // locate xi pixel1 (k)
                    int k_index1 = i_index1*num_sigma_bins+j_index1;
                    if (k_index1 < 0){
                        if (flag_verbose_distortion_plate_ >= 2){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Pixel rejected: Bad indexing: k_index = " << k_index1 << std::endl;
                            }
                        }
                        continue;
                        
                    }
                    
                    // add pixel's weight to the total weight of the bin
                    weight_[k_index1] += spectrum[pixel1].weight();
                    
                    
                    // loop over LyaPixels 2
                    for (int pixel2 = pixel1; pixel2 < spectrum.size(); pixel2 ++){
                    
                        // compute pi and sigma
                        double sigma2 = sqrt(sigma_aux*spectrum[pixel2].dist());
                        if (sigma2 > max_sigma){ // if sigma (r_perp) is too large, discard pixel
                            if (flag_verbose_distortion_plate_ >= 3){
                                #pragma omp critical (cout)
                                {
                                    std::cout << "Pixel rejected: sigma is too large" << std::endl;
                                    std::cout << "sigma = " << sigma2 << std::endl;
                                }
                            }
                            continue;
                        }
                        
                        double pi2 = spectrum[pixel2].dist()-object.dist();
                        if ((pi2 > max_pi) or (pi2 <- max_pi)){ // if pi (r_par) is too large or too small, discard pixel
                            if (flag_verbose_distortion_plate_ >= 3){
                                #pragma omp critical (cout)
                                {
                                    std::cout << "Pixel rejected: abs(pi) is too large" << std::endl;
                                    std::cout << "pi = " << pi2 << std::endl;
                                }
                            }
                            continue;
                        }
                        
                        // locate pi pixel (i_index)
                        int i_index2 = int(pi2/step_pi)+num_pi_bins/2;
                        if (pi2<0.0){
                            i_index2 -= 1;
                        }
                        if (i_index2 < 0){
                            if (flag_verbose_distortion_plate_ >= 2){
                                #pragma omp critical (cout)
                                {
                                    std::cout << "Pixel rejected: Bad indexing: i_index = " << i_index2 << std::endl;
                                }
                            }
                            continue;
                        }
                        
                        // locate sigma pixel (j_index)
                        int j_index2 = int(sigma2/step_sigma);
                        if (j_index2 < 0){
                            if (flag_verbose_distortion_plate_ >= 2){
                                #pragma omp critical (cout)
                                {
                                    std::cout << "Pixel rejected: Bad indexing: j_index = " << j_index2 << std::endl;
                                }
                            }
                            continue;
                        }
                        
                        // locate xi pixel2 (k)
                        int k_index2 = i_index2*num_sigma_bins+j_index2;
                        if (k_index2 < 0){
                            if (flag_verbose_distortion_plate_ >= 2){
                                #pragma omp critical (cout)
                                {
                                    std::cout << "Pixel rejected: Bad indexing: k_index = " << k_index2 << std::endl;
                                }
                            }
                            continue;
                            
                        }
                    
                        AddPair(spectrum[pixel1], spectrum[pixel2], k_index1, k_index2, forest_total_weight, forest_mean_loglam, forest_aux);
                        
                        if (pixel1 != pixel2){
                            AddPair(spectrum[pixel2], spectrum[pixel1], k_index2, k_index1, forest_total_weight, forest_mean_loglam, forest_aux);
                        }
                        
                    }
            
                }
            }
        }
    }
    
    
}

void DistortionPlate::Normalize(){
    /**
     EXPLANATION:
     Normalizes the distortion matrix
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     DistortionPlate
     
     FUNCITONS USED:
     NONE
     */
    
    if (plate_number_ == _NORM_){
        
        CovMat::iterator it;
        for (size_t i = 0; i < num_bins_; i ++){
            for (size_t j = 0; j < num_bins_; j ++){
                it = dist_mat_.find(std::pair<size_t,size_t>(i,j));
                if (it != dist_mat_.end() and weight_[i] != 0.0){
                    (*it).second /= weight_[i];
                }
            }
            //weight_[i] = 1.0;
        }
        
    }
    else{
        // if plate number is not _NORM_, the instance is not supposed to normalize
        std::cout << "Warning : In DistortionPlate::Normalize : Plate number is not set to _NORM_. This DistortionPlate instance should not be normalized. Ignoring..." << std::endl;
    }
    
}

void DistortionPlate::operator+= (const DistortionPlate& other){
    /**
     EXPLANATION:
     Overloads the += operator. Updates the dist_mat_ values by adding those in other.
     
     INPUTS:
     other - a DistortionPlate instance to be added
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     DistortionPlate
     
     FUNCITONS USED:
     NONE
     */
    
    
    // check that both instances have the same number of bins
    if (num_bins_ != other.num_bins()){
        std::cout << "Warning : In DistortionPlate::operator+= : Trying to add DistortionPlates with different number of bins. Ignoring..." << std::endl;
        return;
    }
    
    CovMat::iterator it;
    CovMat::const_iterator other_it;
        
    for (size_t i = 0; i < num_bins_; i++){
        for (size_t j = 0; j < num_bins_; j++){
            it = dist_mat_.find(std::pair<size_t,size_t>(i,j));
            other_it = other.dist_mat_.find(std::pair<size_t,size_t>(i,j));
            if (it != dist_mat_.end() and other_it != other.dist_mat_.end()){
                (*it).second += (*other_it).second;
            }
        }
        weight_[i] += other.weight(i);
    }
    
}

DistortionPlate DistortionPlate::operator- (const DistortionPlate& other){
    /**
     EXPLANATION:
     Overloads the - operator. Returns a new instance with the dist_mat subtracted by those in other.
     
     INPUTS:
     other - a DistortionPlate instance to be subtracted
     
     OUTPUTS:
     temp - a DistortionPlate instance
     
     CLASSES USED:
     DistortionPlate
     
     FUNCITONS USED:
     NONE
     */

    DistortionPlate temp;
    temp = DistortionPlate(plate_number_, num_bins_, plate_neighbours_, flag_verbose_distortion_plate_);
    
    // check that both instances have the same number of bins
    if (num_bins_ != other.num_bins()){
        std::cout << "Warning : In DistortionPlate::operator+= : Trying to add DistortionPlates with different number of bins. Returning zero filled DistortionPlates..." << std::endl;
        return temp;
    }
    
    CovMat::iterator it;
    CovMat::const_iterator other_it;
    
    for (size_t i = 0; i < num_bins_; i++){
        for (size_t j = 0; j < num_bins_; j++){
            it = dist_mat_.find(std::pair<size_t,size_t>(i,j));
            other_it = other.dist_mat_.find(std::pair<size_t,size_t>(i,j));
            if (it != dist_mat_.end() and other_it != other.dist_mat_.end()){
                temp.set_dist_mat(i, j, (*it).second - (*other_it).second);
            }
        }
        temp.set_weight(i, weight_[i] - other.weight(i));
    }
    
    return temp;
    
}
        
DistortionPlate DistortionPlate::operator* (const DistortionPlate& other){
    /**
     EXPLANATION:
     Overloads the * operator. Returns a new instance with the xi_, mean_sigma, mean_pi and weight_ values multiplied by those in other.
     
     INPUTS:
     other - a DistortionPlate instance to be subtracted
     
     OUTPUTS:
     temp - a DistortionPlate instance
     
     CLASSES USED:
     DistortionPlate
     
     FUNCITONS USED:
     NONE
     */
    DistortionPlate temp;
    temp = DistortionPlate(plate_number_, num_bins_, plate_neighbours_, flag_verbose_distortion_plate_);
    
    // check that both instances have the same number of bins
    if (num_bins_ != other.num_bins()){
        std::cout << "Warning : In DistortionPlate::operator+= : Trying to add DistortionPlates with different number of bins. Returning zero filled DistortionPlates..." << std::endl;
        return temp;
    }
    
    CovMat::iterator it;
    CovMat::const_iterator other_it;
    
    for (size_t i = 0; i < num_bins_; i++){
        for (size_t j = 0; j < num_bins_; j++){
            it = dist_mat_.find(std::pair<size_t,size_t>(i,j));
            other_it = other.dist_mat_.find(std::pair<size_t,size_t>(i,j));
            if (it != dist_mat_.end() and other_it != other.dist_mat_.end()){
                temp.set_dist_mat(i, j, (*it).second*(*other_it).second);
            }
        }
        temp.set_weight(i, weight_[i]*other.weight(i));
    }
    
    return temp;
    
}

