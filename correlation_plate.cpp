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
    flag_compute_covariance_ =  input.flag_compute_covariance();
    
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
        /*
         OLD STUFF
        pairs_file_name_ = "/" + input.pairs_file_name() + ToStr(plate_number_);
         */
    }
    
    xi_.resize(num_bins_,0.0);
    mean_pi_.resize(num_bins_,0.0);
    mean_sigma_.resize(num_bins_,0.0);
    weight_.resize(num_bins_,0.0);
    num_averaged_pairs_.resize(num_bins_,0);
    
}

CorrelationPlate::CorrelationPlate(const int plate_number, const int num_bins, const std::string& results, const std::string& pairs_file_name, const std::vector<int>& plate_neighbours, size_t flag_verbose_correlation_plate, size_t flag_write_partial_results){
    /**
     EXPLANATION:
     Cosntructs a CorrelationPlate instance and initializes all its variables
     
     INPUTS:
     input - a Input instance
     results - name of the folder where detailed information will be stored
     plate_number - an integer with the plate number
     plate_neighbours - a vector containing the plate numbers of the neighbouring plates
     flag_verbose_correlation_plate - correlation_plate verbose flag
     flag_write_partial_results_ - flag to write partial results
     
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
        std::cout << "Warining: in CorrelationPlate::set_mean_pi(index,value): The given index is out of bouds, ignoring..." << std::endl;
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
        std::cout << "Warining: in CorrelationPlate::set_mean_sigma(index,value): The given index is out of bouds, ignoring..." << std::endl;
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
        std::cout << "Warining: in CorrelationPlate::set_num_average_pairs(index,value): The given index is out of bouds, ignoring..." << std::endl;
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
        std::cout << "Warining: in CorrelationPlate::set_weight(index,value): The given index is out of bouds, ignoring..." << std::endl;
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
        std::cout << "Warining: in CorrelationPlate::set_xi(index,value): The given index is out of bouds, ignoring..." << std::endl;
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
        std::cout << "Warning : In CorrelationPlate::ComputeCrossCorrelation : Plate number is set to _NORM_. The cross-correlation should not be computed in this CorrelationPlate instance. Ignoring..." << std::endl;
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
        std::cout << "Computing cross-correlation in plate " << plate_number_ << ". In this plate there are " << number_of_objects << " AstroObjects" << std::endl;
    }

    // loop over AstroObjects
    for (size_t i = 0; i < number_of_objects; i ++){
        
        AstroObject object = object_list.list(plate_number_, i);
        
        // checking that the object was loaded successfully
        if (object.dist() == _BAD_DATA_){
            if (flag_verbose_correlation_plate_ >= 2){
                std::cout << "_BAD_DATA_ AstroObject found. Ignoring..." << std::endl;
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
                        std::cout << "_BAD_DATA_ LyaSpectrum found. Ignoring..." << std::endl;
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
                        std::cout << "Spectrum rejected: min_sigma is too large" << std::endl;
                        if (flag_verbose_correlation_plate_ >= 3){
                            std::cout << "min_sigma = " << pair_min_sigma << std::endl;
                            std::cout << "ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                            std::cout << object.angle() << " " << lya_spectrum.angle() << " " << cos_theta << " " << acos(cos_theta) << " " << pair_min_sigma << std::endl;
                        }
                    }
                    continue;
                } 
                if (pair_min_sigma <= 0.1){ // if the spectrum is being paired with its quasar, the whole spectra is discarded
                    if (flag_verbose_correlation_plate_ >= 2){
                        std::cout << "Spectrum rejected: It is being paired with its origin quasar" << std::endl;
                    }
                    continue;
                }
                
                double pair_max_pi = object.dist()-spectrum[0].dist(); // minimum distance obtained for lowest redshift pixel
                if (pair_max_pi < -max_pi) { // if maximum value for pi (r_par) is too small, the whole spectra is discarded
                    if (flag_verbose_correlation_plate_ >= 2){
                        std::cout << "Spectrum rejected: max_pi is too small" << std::endl;
                        if (flag_verbose_correlation_plate_ >= 3){
                            std::cout << "max_pi = " << pair_max_pi << std::endl;
                            std::cout << "ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                            std::cout << object.angle() << " " << lya_spectrum.angle() << " " << cos_theta << " " << acos(cos_theta) << " " << pair_min_sigma << std::endl;
                        }
                    }
                    continue;
                }
                
                double pair_min_pi = object.dist()-spectrum.back().dist(); // maximum distance obtained for highest redshift pixel
                if (pair_min_pi > max_pi){ // if minimum value for pi (r_par) is too high, the whole spectrum is discarded
                    if (flag_verbose_correlation_plate_ >= 2){
                        std::cout << "Spectrum rejected: min_pi is too large" << std::endl;
                        if (flag_verbose_correlation_plate_ >= 3){
                            std::cout << "min_pi = " << pair_min_pi << std::endl;
                            std::cout << "ra_object dec_object ra_spectra dec_spectra cos_theta theta min_sigma" << std::endl;
                            std::cout << object.angle() << " " << lya_spectrum.angle() << " " << cos_theta << " " << acos(cos_theta) << " " << pair_min_sigma << std::endl;
                        }
                    }
                    continue;
                }
                
                // loop over LyaPixels
                for (int p = 0; p < spectrum.size(); p ++){
                    
                    // compute pi and sigma
                    double sigma = sqrt(sigma_aux*spectrum[p].dist());
                    if (sigma > max_sigma){ // if sigma (r_perp) is too large, discard pixel
                        if (flag_verbose_correlation_plate_ >= 3){
                            std::cout << "Pixel rejected: sigma is too large" << std::endl;
                            std::cout << "sigma = " << sigma << std::endl;
                        }
                        continue;
                    }

                    double pi = object.dist()-spectrum[p].dist();
                    if ((pi > max_pi) or (pi <- max_pi)){ // if pi (r_par) is too large or too small, discard pixel
                        if (flag_verbose_correlation_plate_ >= 3){
                            std::cout << "Pixel rejected: abs(pi) is too large" << std::endl;
                            std::cout << "pi = " << pi << std::endl;
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
                            std::cout << "Pixel rejected: Bad indexing: i_index = " << i_index << std::endl;
                        }
                        continue;
                    }

                    
                    // locate sigma pixel (j_index)
                    int j_index = int(sigma/step_sigma);
                    if (j_index < 0){
                        if (flag_verbose_correlation_plate_ >= 2){
                            std::cout << "Pixel rejected: Bad indexing: j_index = " << j_index << std::endl;
                        }
                        continue;
                    }
                    
                    // locate xi pixel (k)
                    int k_index = i_index*num_sigma_bins+j_index;
                    if (k_index < 0){
                        if (flag_verbose_correlation_plate_ >= 2){
                            std::cout << "Pixel rejected: Bad indexing: k_index = " << k_index << std::endl;
                        }
                        continue;
                        
                    }
                    // add contribution to xi in the specified bin
                    AddPair(k_index, spectrum[p], pi, sigma);
                    
                    // write down pair information in bin file
                    if (flag_write_partial_results_ >= 1 or flag_compute_covariance_){
                        KeepPair(k_index, lya_spectrum, p);
                        /* 
                         OLD STUFF
                        SavePair(k_index, object, lya_spectrum, p, pi, sigma);
                         */
                    }
                    
                }
            }
        }
    }
    
    for (size_t i = 0; i < num_bins_; i++){
        if (CorrelationPlate::position_[i] > 0){
            SavePairs(i);
        }
        else{
            std::cout << "in bin " << i << " plate " << plate_number_ << " has either no pairs or else exactly " << max_pairs_ << std::endl;
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

void CorrelationPlate::InitializeStatic(const Input& input){
    /**
     EXPLANATION:
     Initializes the static variables
     
     INPUTS:
     input - a Input instance
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CorrelationPlate
     Pairs
     
     FUNCITONS USED:
     NONE
     */
    
    std::vector<Pair> v;
    v.resize(max_pairs_);
    
    for (size_t i = 0; i < input.num_bins(); i++){
        CorrelationPlate::pairs_information_.push_back(v);
        CorrelationPlate::position_.push_back(0);
    }
    
}

void CorrelationPlate::KeepPair(const int& k_index, const LyaSpectrum& lya_spectrum, const size_t& pixel_number){
    /**
     EXPLANATION:
     Keeps the pair information to save at an appropiate time
     
     INPUTS:
     k_index - an integer specifying the pair's bin
     lya_spectrum - a LyaSpectrum instance
     pixel_number - an unsigned integral with the position of the pixel contributing to the pair
     
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
    
    LyaPixel pixel = lya_spectrum.spectrum(pixel_number);
    SpherePoint angle = lya_spectrum.angle();
    
    Pair pair(angle.ra(), angle.dec(), pixel_number, pixel.dist(), pixel.weight());
    
    CorrelationPlate::pairs_information_[k_index][CorrelationPlate::position_[k_index]] = pair;
    CorrelationPlate::position_[k_index] ++;
    
    if (CorrelationPlate::position_[k_index] == CorrelationPlate::max_pairs_){
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
                std::cout << "Zero Division Error in CorrelationPlate::Normalize. Aborting..." << std::endl;
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
        std::cout << "Warning : In CorrelationPlate::Normalize : Plate number is not set to _NORM_. This CorrelationPlate instance should not be normalized. Ignoring..." << std::endl;
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
    filename = results_ + ToStr(k_index) + ".fits";
    
    // construct fits object
    std::auto_ptr<CCfits::FITS> pFits(0);
    
    try{
        pFits.reset(new CCfits::FITS(filename,CCfits::Write));
    }
    catch(CCfits::FITS::CantOpen){
        throw "Error : In CorrelationPlate::SavePairs : Unable to open file: \n " + ToStr(filename);
    }
    // read table from file
    CCfits::ExtHDU& table = (*pFits).extension(pairs_file_name_);
    
    long NAXIS2 = table.axis(1);
    size_t size = CorrelationPlate::position_[k_index];
    
    // prepare variables to write in the table
    
    std::valarray<double> spectrum_ra(size);
    std::valarray<double> spectrum_dec(size);
    std::valarray<double> pixel_dist(size);
    std::valarray<int> pixel_number(size);
    std::valarray<double> pixel_weight(size);
    
    for (size_t i = 0; i < size; i++){
        spectrum_ra[i] = CorrelationPlate::pairs_information_[k_index][i].spectrum_ra();
        spectrum_dec[i] = CorrelationPlate::pairs_information_[k_index][i].spectrum_dec();
        pixel_dist[i] = CorrelationPlate::pairs_information_[k_index][i].pixel_dist();
        pixel_number[i] = CorrelationPlate::pairs_information_[k_index][i].pixel_number();
        pixel_weight[i] = CorrelationPlate::pairs_information_[k_index][i].pixel_weight();
    }
    
    // write data in the table
    table.column("spectrum RA").write(spectrum_ra,NAXIS2+1);
    table.column("spectrum DEC").write(spectrum_dec,NAXIS2+1);
    table.column("pixel dist").write(pixel_dist,NAXIS2+1);
    table.column("pixel number").write(pixel_number,NAXIS2+1);
    table.column("pixel weight").write(pixel_weight,NAXIS2+1);
    
    
    CorrelationPlate::position_[k_index] = 0;
}

/*void CorrelationPlate::SavePair(const int& k_index, const AstroObject& object, const LyaSpectrum& lya_spectrum, const size_t& p, const double& pi, const double& sigma){
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
    

    /*
     OLD STUFF
     filename = results_ + ToStr(k_index) + pairs_file_name_;
    
    std::ofstream bin_file;
    bin_file.open(filename.c_str(),std::ofstream::app); // opens the file to append content
    if (bin_file.is_open()){
        LyaPixel pixel = lya_spectrum.spectrum(p);
        // checking that the object was loaded successfully
        if (pixel.z() == _BAD_DATA_){
            if (flag_verbose_correlation_plate_ >= 2){
                std::cout << "_BAD_DATA_ LyaPixel found. Ignoring..." << std::endl;
            }
            return;
        }
        
        bin_file << lya_spectrum.angle() << " " << pixel.dist() << " " << p << " " << pixel.weight() << std::endl;
        bin_file.close();
    }
    else{
        std::cout << "Error : In CorrelationPlate::SavePair : Unable to open file:" << std::endl << filename << std::endl;
    }

}*/

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
        std::cout << "Warning : In CorrelationPlate::operator+= : Trying to add CorrelationPlates with different number of bins. Ignoring..." << std::endl;
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
    temp = CorrelationPlate(plate_number_, num_bins_, results_, pairs_file_name_, plate_neighbours_, flag_verbose_correlation_plate_, flag_write_partial_results_);
    
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
        std::cout << "Warning : In CorrelationPlate::operator- : Trying to subtract CorrelationPlates with different number of bins. Returning zero filled CorrelationPlates..." << std::endl;
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
    temp = CorrelationPlate(plate_number_, num_bins_, results_, pairs_file_name_, plate_neighbours_, flag_verbose_correlation_plate_, flag_write_partial_results_);
    
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
        std::cout << "Warning : In CorrelationPlate::operator* : Trying to multiply CorrelationPlates with different number of bins. Returning zero filled CorrelationPlates..." << std::endl;
    }
            
    return temp;
    
}

size_t CorrelationPlate::max_pairs_ = 10000;
std::vector<size_t> CorrelationPlate::position_;
std::vector<std::vector<Pair> > CorrelationPlate::pairs_information_;

