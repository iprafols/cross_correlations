/**
 covariance_plate.cpp
 Purpose: This files contains the body for the functions defined in covariance_plate.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 10/02/2014
 */

#include "covariance_plate.h"

CovariancePlate::CovariancePlate(int bad_data){
    /**
     EXPLANATION:
     Constructs a CovariancePlate instance and initializes all its variables
     
     INPUTS:
     bad_data - a double valued _BAD_DATA_INT_
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovariancePlate
     
     FUNCITONS USED:
     ToStr
     */
    if (bad_data != _BAD_DATA_INT_){
        std::cout << "Error while initializing a CovariancePlate 'bad data' instance" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    
    plate_number_ = _BAD_DATA_INT_;
    
}

CovariancePlate::CovariancePlate(const Input& input, const int plate_number, const std::vector<int>& plate_neighbours){
    /**
     EXPLANATION:
     Cosntructs a CovariancePlate instance and initializes all its variables
     
     INPUTS:
     input - a Input instance
     plate_number - an integer with the plate number
     plate_neighbours - a vector containing the plate numbers of the neighbouring plates
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovariancePlate
     
     FUNCITONS USED:
     ToStr
     */
    
    // set flags from input
    flag_verbose_covariance_plate_ = input.flag_verbose_covariance_plate();
    
    plate_number_ = plate_number;
    plate_neighbours_ = plate_neighbours;
    num_bins_ = input.num_bins();
    pairs_file_name_ = ToStr(plate_number_);
    
    // initialize covariance matrix computation
    for (size_t i = 0; i < num_bins_; i++){
        for (size_t j = i; j < num_bins_; j++){
            cov_mat_[std::pair<size_t,size_t>(i,j)] = 0.0;
            //weight_[std::pair<size_t,size_t>(i,j)] = 0.0;
        }
    }
    weight_.resize(num_bins_, 0.0);
    
}

CovariancePlate::CovariancePlate(const int plate_number, const int num_bins, const std::vector<int>& plate_neighbours, size_t flag_verbose_covariance_plate){
    /**
     EXPLANATION:
     Cosntructs a CovariancePlate instance and initializes all its variables
     
     INPUTS:
     input - a Input instance
     plate_number - an integer with the plate number
     plate_neighbours - a vector containing the plate numbers of the neighbouring plates
     flag_verbose_covariance_plate_ - correlation_plate verbose flag
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovariancePlate
     
     FUNCITONS USED:
     ToStr
     */
    
    flag_verbose_covariance_plate_ = flag_verbose_covariance_plate;
    
    plate_number_ = plate_number;
    plate_neighbours_ = plate_neighbours;
    num_bins_ = num_bins;
    
    // initialize covariance matrix computation
    for (size_t i = 0; i < num_bins_; i++){
        for (size_t j = i; j < num_bins_; j++){
            cov_mat_[std::pair<size_t,size_t>(i,j)] = 0.0;
            //weight_[std::pair<size_t,size_t>(i,j)] = 0.0;
        }
    }
    weight_.resize(num_bins_, 0.0);
}

double CovariancePlate::weight(size_t i) const{
    /**
     EXPLANATION:
     Access function for weight_
     
     INPUTS:
     index - index of the selected weight_ element
     
     OUTPUTS:
     the weight for the selected bin
     
     CLASSES USED:
     CorrelationPlate
     
     FUNCITONS USED:
     NONE
     */
    if (i < weight_.size()){
        return weight_[i];
    }
    else{
        std::cout << "Warining: in CovariancePlate::set_weight(i, value): The given index is out of bouds, returning BAD_DATA..." << std::endl;
        return _BAD_DATA_;
    }
}

void CovariancePlate::set_cov_mat(size_t i, size_t j, double value){
    /**
     EXPLANATION:
     Set function for cov_mat_
     
     INPUTS:
     i,j - indexs of the selected cov_mat_ element
     value - element's new value
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovariancePlate
     
     FUNCITONS USED:
     NONE
     */
    
    CovMat::iterator it = cov_mat_.find(std::pair<size_t, size_t>(i,j));
    if (it != cov_mat_.end()){
        (*it).second = value;
    }
    else{
        std::cout << "Warining: in CovariancePlate::set_cov_mat(i, j, value): The given index is out of bouds, ignoring..." << std::endl;
    }
}

//void CovariancePlate::set_weight(size_t i, size_t j, double value){
void CovariancePlate::set_weight(size_t i, double value){
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
    
    if (i < weight_.size()){
        weight_[i] = value;
    }
    else{
        std::cout << "Warining: in CovariancePlate::set_weight(i, value): The given index is out of bouds, ignoring..." << std::endl;
    }
    /*CovMat::iterator it = weight_.find(std::pair<size_t, size_t>(i,j));
    if (it != weight_.end()){
        (*it).second = value;
    }
    else{
        std::cout << "Warining: in CovariancePlate::set_weight(i, j, value): The given index is out of bouds, ignoring..." << std::endl;
    }*/
}

void CovariancePlate::AddPair(const LyaPixel& pixel1, const LyaPixel& pixel2, const size_t& i, const size_t& j, const LyaAutoInterpolationMap& lya_auto){
    /**
     EXPLANATION:
     Adds pair contribution to the covariance matrix in the specified bin
     
     INPUTS:
     pixel1,pixel2 - LyaPixel instances to add the contribution from
     i,j - covariance matrix element to add the contributions to
     lya_auto - LyaAutoInterpolationMap for the pixel's separation
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovariancePlate
     LyaAutoInterpolationMap
     LyaPixel
     
     FUNCITONS USED:
     NONE
     */
    
    CovMat::iterator it, it_weight;
    if (i <= j){
        it = cov_mat_.find(std::pair<size_t, size_t>(i, j));
        //it_weight = weight_.find(std::pair<size_t, size_t>(i, j));
    }
    else{
        it = cov_mat_.find(std::pair<size_t, size_t>(j, i));
        //it_weight = weight_.find(std::pair<size_t, size_t>(j, i));
    }
    if (it == cov_mat_.end()){ //or it_weight == weight_.end()){
        std::cout << "Warning : In CovarianceMatrix::ComputeCovMat : Element " << i << ", " << j << " of the covariance matrix not found. Ignoring..."  << std::endl;
        return;
    }
    
    // add to covariance matrix
    double weight = pixel1.weight()*pixel2.weight();
    double add = lya_auto.LinearInterpolation((pixel1.z()+pixel2.z())/2.0);

    if (add != _BAD_DATA_){
        (*it).second += add*weight;
    }
    
}

void CovariancePlate::ComputeCovMat(const AstroObjectDataset& object_list, const SpectraDataset& spectra_list, const Input& input, const std::vector<LyaAutoInterpolationMap>& lya_auto){
    /**
     EXPLANATION:
     Computes the covariance matrix
     
     INPUTS:
     object_list - an AstroObjectDataset instance
     spectra_list - a SpectraDataset instance
     input - a Input instance to load bin settings
     lya_auto - an IntepolationMap containing the 1D Lya forest auto-correlation
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObjectDataset
     CovariancePlate
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
            std::cout << "Warning : In CovariancePlate::ComputeCovMat : Plate number is set to _NORM_. The covariance should not be computed in this CovariancePlate instance. Ignoring..." << std::endl;
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
    
    if (flag_verbose_covariance_plate_ >= 1){
        #pragma omp critical (cout)
        {
            std::cout << "Computing covariance_matrix in plate " << plate_number_ << ". In this plate there are " << number_of_spectra << " LyaSpectra" << std::endl;
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
            if (flag_verbose_covariance_plate_ >= 2){
                #pragma omp critical (cout)
                {
                    std::cout << "_BAD_DATA_ LyaSpectrum found. Ignoring..." << std::endl;
                }
            }
            continue;
        }
        std::vector<LyaPixel> spectrum = lya_spectrum.spectrum();
        
        // loop over neighbouring plates 1
        for (size_t plate_neighbours_num1 = 0; plate_neighbours_num1 < plate_neighbours_.size(); plate_neighbours_num1 ++){
            
            size_t number_of_objects1 = object_list.num_objects_in_plate(plate_neighbours_[plate_neighbours_num1]);
            if (number_of_objects1 == 0){
                continue;
            }
            
            // loop over AstroObjects 1
            for (size_t object_num1 = 0; object_num1 < number_of_objects1; object_num1 ++){
                
                AstroObject object1 = object_list.list(plate_neighbours_[plate_neighbours_num1], object_num1);
                
                // checking that the object was loaded successfully
                if (object1.dist() == _BAD_DATA_){
                    if (flag_verbose_covariance_plate_ >= 2){
                        #pragma omp critical (cout)
                        {
                            std::cout << "_BAD_DATA_ AstroObject found. Ignoring..." << std::endl;
                        }
                    }
                    continue;
                }
                
                double sigma_aux_o_1 = 2.0*object1.dist(); // auxiliar variable to compute sigma values
                
                // compute angular separation
                double cos_theta1 = object1.angle().CosAngularDistance(lya_spectrum.angle());
                
                double sigma_aux1;
                if (cos_theta1 == 1.0){
                    sigma_aux1 = 0.0;
                }
                else{
                    sigma_aux1 = sigma_aux_o_1*(1.0-cos_theta1); // auxiliar variable to compute sigma values
                }
                
                // check if all pixels in the spectrum are too far apart
                double pair_min_sigma1 = sqrt(sigma_aux1*spectrum[0].dist()); // minimum distance obtained for lowest redshift pixel
                if (pair_min_sigma1 > max_sigma){ // if the minimum value for sigma (r_perp) is too large, the whole spectra is discarded
                    if (flag_verbose_covariance_plate_ >= 2){
                        #pragma omp critical (cout)
                        {
                            std::cout << "Spectrum rejected: min_sigma is too large" << std::endl;
                        }
                        if (flag_verbose_covariance_plate_ >= 3){
                            #pragma omp critical (cout)
                            {
                                std::cout << "min_sigma = " << pair_min_sigma1 << std::endl;
                                std::cout << "ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                                std::cout << object1.angle() << " " << lya_spectrum.angle() << " " << cos_theta1 << " " << acos(cos_theta1) << " " << pair_min_sigma1 << std::endl;
                            }
                        }
                    }
                    continue;
                }
                // if the spectrum is being paired with its quasar, the whole spectra is discarded
                if (pair_min_sigma1 <= 0.1){
                    if (flag_verbose_covariance_plate_ >= 2){
                        #pragma omp critical (cout)
                        {
                            std::cout << "Spectrum rejected: It is being paired with its origin quasar" << std::endl;
                        }
                    }
                    continue;
                }
                
                double pair_max_pi1 = spectrum.back().dist() - object1.dist(); // minimum distance obtained for lowest redshift pixel
                if (pair_max_pi1 < -max_pi) { // if maximum value for pi (r_par) is too small, the whole spectra is discarded
                    if (flag_verbose_covariance_plate_ >= 2){
                        #pragma omp critical (cout)
                        {
                            std::cout << "Spectrum rejected: max_pi is too small" << std::endl;
                        }
                        if (flag_verbose_covariance_plate_ >= 3){
                            #pragma omp critical (cout)
                            {
                                std::cout << "max_pi = " << pair_max_pi1 << std::endl;
                                std::cout << "ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                                std::cout << object1.angle() << " " << lya_spectrum.angle() << " " << cos_theta1 << " " << acos(cos_theta1) << " " << pair_max_pi1 << std::endl;
                            }
                        }
                    }
                    continue;
                }
                
                double pair_min_pi1 = spectrum[0].dist() - object1.dist(); // maximum distance obtained for highest redshift pixel
                if (pair_min_pi1 > max_pi){ // if minimum value for pi (r_par) is too high, the whole spectrum is discarded
                    if (flag_verbose_covariance_plate_ >= 2){
                        #pragma omp critical (cout)
                        {
                            std::cout << "Spectrum rejected: min_pi is too large" << std::endl;
                        }
                        if (flag_verbose_covariance_plate_ >= 3){
                            #pragma omp critical (cout)
                            {
                                std::cout << "min_pi = " << pair_min_pi1 << std::endl;
                                std::cout << "ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                                std::cout << object1.angle() << " " << lya_spectrum.angle() << " " << cos_theta1 << " " << acos(cos_theta1) << " " << pair_min_pi1 << std::endl;
                            }
                        }
                    }
                    continue;
                }
                
                // loop over LyaPixels1
                for (int pixel1 = 0; pixel1 < spectrum.size(); pixel1 ++){
                    
                    // cehck that the weight is not zero
                    if (spectrum[pixel1].weight() == 0.0){
                        continue;
                    }
                    
                    // compute pi and sigma
                    double sigma1 = sqrt(sigma_aux1*spectrum[pixel1].dist());
                    if (sigma1 > max_sigma){ // if sigma (r_perp) is too large, discard pixel
                        if (flag_verbose_covariance_plate_ >= 3){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Pixel rejected: sigma is too large" << std::endl;
                                std::cout << "sigma = " << sigma1 << std::endl;
                            }
                        }
                        continue;
                    }
                    
                    double pi1 = spectrum[pixel1].dist()-object1.dist();
                    if ((pi1 > max_pi) or (pi1 <- max_pi)){ // if pi (r_par) is too large or too small, discard pixel
                        if (flag_verbose_covariance_plate_ >= 3){
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
                        if (flag_verbose_covariance_plate_ >= 2){
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
                        if (flag_verbose_covariance_plate_ >= 2){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Pixel rejected: Bad indexing: j_index = " << j_index1 << std::endl;
                            }
                        }
                        continue;
                    }
                    
                    // locate xi pixel (k)
                    int k_index1 = i_index1*num_sigma_bins+j_index1;
                    if (k_index1 < 0){
                        if (flag_verbose_covariance_plate_ >= 2){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Pixel rejected: Bad indexing: k_index = " << k_index1 << std::endl;
                            }
                        }
                        continue;
                        
                    }
                    
                    // pair found, add weight to denominator
                    weight_[k_index1] += spectrum[pixel1].weight();
                    
                    // loop over neighbouring plates 2
                    for (size_t plate_neighbours_num2 = plate_neighbours_num1; plate_neighbours_num2 < plate_neighbours_.size(); plate_neighbours_num2 ++){
                        
                        size_t number_of_objects2 = object_list.num_objects_in_plate(plate_neighbours_[plate_neighbours_num2]);
                        if (number_of_objects2 == 0){
                            continue;
                        }
                        
                        // loop over AstroObjects 2
                        for (size_t object_num2 = 0; object_num2 < number_of_objects2; object_num2 ++){
                            
                            AstroObject object2 = object_list.list(plate_neighbours_[plate_neighbours_num2], object_num2);
                            
                            // checking that the object was loaded successfully
                            if (object2.dist() == _BAD_DATA_){
                                if (flag_verbose_covariance_plate_ >= 2){
                                    #pragma omp critical (cout)
                                    {
                                        std::cout << "_BAD_DATA_ AstroObject found. Ignoring..." << std::endl;
                                    }
                                }
                                continue;
                            }
                            
                            double sigma_aux_o_2 = 2.0*object2.dist(); // auxiliar variable to compute sigma values
                            
                            
                            // compute angular separation
                            double cos_theta2 = object2.angle().CosAngularDistance(lya_spectrum.angle());
                            
                            double sigma_aux2;
                            if (cos_theta2 == 1.0){
                                sigma_aux2 = 0.0;
                            }
                            else{
                                sigma_aux2 = sigma_aux_o_2*(1.0-cos_theta2); // auxiliar variable to compute sigma values
                            }
                            
                            // check if all pixels in the spectrum are too far apart
                            double pair_min_sigma2 = sqrt(sigma_aux2*spectrum[0].dist()); // minimum distance obtained for lowest redshift pixel
                            if (pair_min_sigma2 > max_sigma){ // if the minimum value for sigma (r_perp) is too large, the whole spectra is discarded
                                if (flag_verbose_covariance_plate_ >= 2){
                                    #pragma omp critical (cout)
                                    {
                                        std::cout << "Spectrum rejected: min_sigma is too large" << std::endl;
                                    }
                                    if (flag_verbose_covariance_plate_ >= 3){
                                        #pragma omp critical (cout)
                                        {
                                            std::cout << "min_sigma = " << pair_min_sigma2 << std::endl;
                                            std::cout << "ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                                            std::cout << object2.angle() << " " << lya_spectrum.angle() << " " << cos_theta2 << " " << acos(cos_theta2) << " " << pair_min_sigma2 << std::endl;
                                        }
                                    }
                                }
                                continue;
                            }
                            // if the spectrum is being paired with its quasar, the whole spectra is discarded
                            if (pair_min_sigma2 <= 0.1){
                                if (flag_verbose_covariance_plate_ >= 2){
                                    #pragma omp critical (cout)
                                    {
                                        std::cout << "Spectrum rejected: It is being paired with its origin quasar" << std::endl;
                                    }
                                }
                                continue;
                            }
                            
                            double pair_max_pi2 = spectrum.back().dist() - object2.dist(); // minimum distance obtained for lowest redshift pixel
                            if (pair_max_pi2 < -max_pi) { // if maximum value for pi (r_par) is too small, the whole spectra is discarded
                                if (flag_verbose_covariance_plate_ >= 2){
                                    #pragma omp critical (cout)
                                    {
                                        std::cout << "Spectrum rejected: max_pi is too small" << std::endl;
                                    }
                                    if (flag_verbose_covariance_plate_ >= 3){
                                        #pragma omp critical (cout)
                                        {
                                            std::cout << "max_pi = " << pair_max_pi2 << std::endl;
                                            std::cout << "ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                                            std::cout << object2.angle() << " " << lya_spectrum.angle() << " " << cos_theta2 << " " << acos(cos_theta2) << " " << pair_max_pi2 << std::endl;
                                        }
                                    }
                                }
                                continue;
                            }
                            
                            double pair_min_pi2 = spectrum[0].dist() - object2.dist(); // maximum distance obtained for highest redshift pixel
                            if (pair_min_pi2 > max_pi){ // if minimum value for pi (r_par) is too high, the whole spectrum is discarded
                                if (flag_verbose_covariance_plate_ >= 2){
                                    #pragma omp critical (cout)
                                    {
                                        std::cout << "Spectrum rejected: min_pi is too large" << std::endl;
                                    }
                                    if (flag_verbose_covariance_plate_ >= 3){
                                        #pragma omp critical (cout)
                                        {
                                            std::cout << "min_pi = " << pair_min_pi2 << std::endl;
                                            std::cout << "ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                                            std::cout << object2.angle() << " " << lya_spectrum.angle() << " " << cos_theta2 << " " << acos(cos_theta2) << " " << pair_min_pi2 << std::endl;
                                        }
                                    }
                                }
                                continue;
                            }
                            
                            // loop over LyaPixels2
                            for (int pixel2 = pixel1; pixel2 < spectrum.size() and pixel2 <= pixel1 + pixels_separation; pixel2 ++){
                                
                                // cehck that the weight is not zero
                                if (spectrum[pixel2].weight() == 0.0){
                                    continue;
                                }
                                
                                // compute pi and sigma
                                double sigma2 = sqrt(sigma_aux2*spectrum[pixel2].dist());
                                if (sigma2 > max_sigma){ // if sigma (r_perp) is too large, discard pixel
                                    if (flag_verbose_covariance_plate_ >= 3){
                                        #pragma omp critical (cout)
                                        {
                                            std::cout << "Pixel rejected: sigma is too large" << std::endl;
                                            std::cout << "sigma = " << sigma2 << std::endl;
                                        }
                                    }
                                    continue;
                                }
                                
                                double pi2 = spectrum[pixel2].dist()-object2.dist();
                                if ((pi2 > max_pi) or (pi2 <- max_pi)){ // if pi (r_par) is too large or too small, discard pixel
                                    if (flag_verbose_covariance_plate_ >= 3){
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
                                    if (flag_verbose_covariance_plate_ >= 2){
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
                                    if (flag_verbose_covariance_plate_ >= 2){
                                        #pragma omp critical (cout)
                                        {
                                            std::cout << "Pixel rejected: Bad indexing: j_index = " << j_index2 << std::endl;
                                        }
                                    }
                                    continue;
                                }
                                
                                // locate xi pixel (k)
                                int k_index2 = i_index2*num_sigma_bins+j_index2;
                                if (k_index2 < 0){
                                    if (flag_verbose_covariance_plate_ >= 2){
                                        #pragma omp critical (cout)
                                        {
                                            std::cout << "Pixel rejected: Bad indexing: k_index = " << k_index2 << std::endl;
                                        }
                                    }
                                    continue;
                                    
                                }
                                
                                // avoid double counting
                                if (pixel2 == pixel1 and plate_neighbours_num2 == plate_neighbours_num1 and object_num2 < object_num1){
                                    continue;
                                }
                                
                                // compute covariance matrix contribution
                                AddPair(spectrum[pixel1],spectrum[pixel2],k_index1,k_index2, lya_auto[pixel2 - pixel1]);
                            }
                        }
                    }
                    

                }
            }
        }
    }
    
    
}

void CovariancePlate::Normalize(){
    /**
     EXPLANATION:
     Normalizes the covariance matrix
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovariancePlate
     
     FUNCITONS USED:
     NONE
     */
    
    if (plate_number_ == _NORM_){
        
        CovMat::iterator it, it_weight;
        for (size_t i = 0; i < num_bins_; i ++){
            for (size_t j = i; j < num_bins_; j ++){
                it = cov_mat_.find(std::pair<size_t,size_t>(i,j));
                if (weight_[i] != 0.0 and weight_[j] != 0.0){
                    (*it).second /= weight_[i];
                    (*it).second /= weight_[j];
                }
                else{
                    (*it).second == 0.0;
                }
                //it_weight = weight_.find(std::pair<size_t,size_t>(i,j));
                //if (it != cov_mat_.end() and it_weight != weight_.end()){
                //    (*it).second /= (*it_weight).second;
                //}
            }
        }
        
    }
    else{
        // if plate number is not _NORM_, the instance is not supposed to normalize
        std::cout << "Warning : In CovariancePlate::Normalize : Plate number is not set to _NORM_. This CovariancePlate instance should not be normalized. Ignoring..." << std::endl;
    }
    
}

void CovariancePlate::SaveCovMat(const Input& input){
    /**
     EXPLANATION:
     Saves the covariance matrix measured in a specific plate
     
     INPUTS:
     input - a Input instance
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovariancePlate
     
     FUNCITONS USED:
     NONE
     */
    if (plate_number_ == _NORM_){
        // if plate number is _NORM_, the instance is not supposed to be saved using this function
        std::cout << "Warning : In CovariancePlate::SaveCovMat : Plate number is not set to _NORM_. This CovariancePlate instance should not be saved using this function. Ignoring..." << std::endl;
    }
    else{
        std::string filename;
        CovMat::iterator it, it_weight;
        
        if (flag_verbose_covariance_plate_ >= 1){
            std::cout << "Saving covariance matrix from plate " << pairs_file_name_ << std::endl;
        }
        
        filename = input.results() + "plate_" + pairs_file_name_ + ".cov";
        {
            std::ofstream file(filename.c_str(),std::ofstream::trunc);
            if (file.is_open()){
                if (flag_verbose_covariance_plate_ >= 2){
                    std::cout << "Saving full covariance matrix" << std::endl;
                }
                for (size_t i = 0; i < num_bins_; i ++){
                    for (size_t j = i; j < num_bins_; j ++){
                        it = cov_mat_.find(std::pair<size_t, size_t>(i, j));
                        
                        if (it != cov_mat_.end() and weight_[i] != 0.0 and weight_[j] != 0.0){
                            file << i << " " << j << " " << (*it).second/weight_[i]/weight_[j] << std::endl;
                        }
                        else{
                            file << i << " " << j << " NaN" << std::endl;
                        }
                    }
                }
                
                file.close();
            }
            else{
                std::cout << "Error : In CovariancePlate::SaveCovMat : Unable to open file:" << std::endl << filename << std::endl;
            }
        }
        
    }
    
}

void CovariancePlate::operator+= (const CovariancePlate& other){
    /**
     EXPLANATION:
     Overloads the += operator. Updates the xi_, mean_sigma, mean_pi and weight_ or cov_mat_ values by adding those in other.
     
     INPUTS:
     other - a CovariancePlate instance to be added
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CovariancePlate
     
     FUNCITONS USED:
     NONE
     */
    
    
    // check that both instances have the same number of bins
    if (num_bins_ != other.num_bins()){
        std::cout << "Warning : In CovariancePlate::operator+= : Trying to add CovariancePlates with different number of bins. Ignoring..." << std::endl;
        return;
    }
    
    CovMat::iterator it, it_weight;
    CovMat::const_iterator other_it, other_it_weight;
        
    for (size_t i = 0; i < num_bins_; i++){
        for (size_t j = i; j < num_bins_; j++){
            it = cov_mat_.find(std::pair<size_t, size_t>(i, j));
            //it_weight = weight_.find(std::pair<size_t, size_t>(i, j));
            other_it = other.cov_mat_.find(std::pair<size_t, size_t>(i, j));
            //other_it_weight = other.weight_.find(std::pair<size_t, size_t>(i, j));
            
            if (it != cov_mat_.end() and other_it != other.cov_mat_.end()){
                (*it).second += (*other_it).second;
            }
            /*if (it_weight != weight_.end() and other_it_weight != other.weight_.end()){
                (*it_weight).second += (*other_it_weight).second;
            }*/
        }
    }
    for (size_t index = 0; index < weight_.size(); index ++){
        weight_[index] += other.weight(index);
    }
    
}

CovariancePlate CovariancePlate::operator- (const CovariancePlate& other){
    /**
     EXPLANATION:
     Overloads the - operator. Returns a new instance with the xi_, mean_sigma, mean_pi and weight_ values or cov_mat subtracted by those in other.
     
     INPUTS:
     other - a CovariancePlate instance to be subtracted
     
     OUTPUTS:
     temp - a CovariancePlate instance
     
     CLASSES USED:
     CovariancePlate
     
     FUNCITONS USED:
     NONE
     */

    CovariancePlate temp;
    temp = CovariancePlate(plate_number_, num_bins_, plate_neighbours_, flag_verbose_covariance_plate_);
    
    // check that both instances have the same number of bins
    if (num_bins_ != other.num_bins()){
        std::cout << "Warning : In CovariancePlate::operator+= : Trying to add CovariancePlates with different number of bins. Returning zero filled CovariancePlates..." << std::endl;
        return temp;
    }
    
    CovMat::iterator it, it_weight;
    CovMat::const_iterator other_it, other_it_weight;
    
    for (size_t i = 0; i < num_bins_; i++){
        for (size_t j = i; j < num_bins_; j++){
            it = cov_mat_.find(std::pair<size_t, size_t>(i, j));
            //it_weight = weight_.find(std::pair<size_t, size_t>(i, j));
            other_it = other.cov_mat_.find(std::pair<size_t, size_t>(i, j));
            //other_it_weight = other.weight_.find(std::pair<size_t, size_t>(i, j));
            
            if (it != cov_mat_.end() and other_it != other.cov_mat_.end()){
                temp.set_cov_mat(i, j, (*it).second - (*other_it).second);
            }
            //if (it_weight != weight_.end() and other_it_weight != other.weight_.end()){
            //    temp.set_weight(i, j, (*it_weight).second - (*other_it_weight).second);
            //}
        }
    }
    
    for (size_t index = 0; index < weight_.size(); index ++){
        temp.set_weight(index, other.weight(index) - weight_[index]);
    }
    
    return temp;
    
}
        
CovariancePlate CovariancePlate::operator* (const CovariancePlate& other){
    /**
     EXPLANATION:
     Overloads the * operator. Returns a new instance with the xi_, mean_sigma, mean_pi and weight_ values multiplied by those in other.
     
     INPUTS:
     other - a CovariancePlate instance to be subtracted
     
     OUTPUTS:
     temp - a CovariancePlate instance
     
     CLASSES USED:
     CovariancePlate
     
     FUNCITONS USED:
     NONE
     */
    CovariancePlate temp;
    temp = CovariancePlate(plate_number_, num_bins_, plate_neighbours_, flag_verbose_covariance_plate_);
    
    // check that both instances have the same number of bins
    if (num_bins_ != other.num_bins()){
        std::cout << "Warning : In CovariancePlate::operator+= : Trying to add CovariancePlates with different number of bins. Returning zero filled CovariancePlates..." << std::endl;
        return temp;
    }
    
    CovMat::iterator it, it_weight;
    CovMat::const_iterator other_it, other_it_weight;
    
    for (size_t i = 0; i < num_bins_; i++){
        for (size_t j = i; j < num_bins_; j++){
            it = cov_mat_.find(std::pair<size_t, size_t>(i, j));
            //it_weight = weight_.find(std::pair<size_t, size_t>(i, j));
            other_it = other.cov_mat_.find(std::pair<size_t, size_t>(i, j));
            //other_it_weight = other.weight_.find(std::pair<size_t, size_t>(i, j));
            
            if (it != cov_mat_.end() and other_it != other.cov_mat_.end()){
                temp.set_cov_mat(i, j, (*it).second*(*other_it).second);
            }
            //if (it_weight != weight_.end() and other_it_weight != other.weight_.end()){
            //    temp.set_weight(i, j, (*it_weight).second*(*other_it_weight).second);
            //}
        }
    }
    for (size_t index = 0; index < weight_.size(); index ++){
        temp.set_weight(index, other.weight(index)*weight_[index]);
    }
    
    return temp;
    
}

double CovariancePlate::half_gamma_ = 3.8/2.0;
double CovariancePlate::one_plus_z0_to_the_half_gamma_ = pow(1+2.25,CovariancePlate::half_gamma_);

