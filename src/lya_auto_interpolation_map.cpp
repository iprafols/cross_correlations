/**
 lya_auto_interpolation_table.cpp
 Purpose: This files contains the body for the functions defined in lya_auto_interpolation_table.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 09/30/2014
 */

#include "lya_auto_interpolation_map.h"

LyaAutoInterpolationMap::LyaAutoInterpolationMap(const Input& input, const int& n){
    /**
     EXPLANATION:
     Cosntructs a LyaAutoInterpolationMap instance and initializes its variables
     
     INPUTS:
     input - object of type Input
     n - distance between pixels (in number of pixels)
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Input
     InterpolationMap
     
     FUNCITONS USED:
     NONE
     */
    
    // TODO: remove old stuff
    /*double lya_pixel_width = input.lya_pixel_width();
    double sigma_psf = input.sigma_psf();
    
    double aux_z = 0.0;
    double value;
    std::vector<std::pair<double, double> > k_pk;
     */
    std::string line;
    bool header = true;
    
    std::ifstream file (input.lya_auto_correlation().c_str());
    if (file.is_open()){
        while (getline(file,line)){
            std::vector<std::string> cols = Split(line," ");
            
            if (header){
                header = false;
                continue;
            }
            if (cols.size() < n + 2){
                throw std::invalid_argument( "invalid distance between pixels, check that compute_lya_1d.run was executed with the correct value for pixels_separation" );
            }
            interpolation_map_[atof(cols[0].c_str())] = atof(cols[n + 1].c_str());
            
            
            // TODO: remove old stuff
            /*if (double(atof(cols[0].c_str())) != aux_z){
                // integrate pk and add autocorrelation to interpolation map
                if (not k_pk.empty()){
                    value = IntegratePk(k_pk, n, lya_pixel_width, sigma_psf);
                    if (value != _BAD_DATA_){
                        interpolation_map_[aux_z] = value;
                    }
                    k_pk.clear();
                }
                // update redshift value
                aux_z = double(atof(cols[0].c_str()));
            
            }
            
            // add pair of (k, p(k)) to vector
            k_pk.push_back(std::pair<double, double>(double(atof(cols[1].c_str())), double(atof(cols[2].c_str()))));
            */
            
        }
        
        // TODO: remove old stuff
        /*if (input.flag_project_deltas()){
            ProjectLyaAuto();
        }*/
    }
    else{
        std::cout << "Error: in LyaAutoInterpolationMap::LyaAutoInterpolationMap : Could not read file" << std::endl;
        std::cout << input.lya_auto_correlation() << std::endl;
    }
}

double LyaAutoInterpolationMap::IntegratePk(const std::vector<std::pair<double, double> >& vec, const int& n, const double& lya_pixel_width, const double& sigma_psf){
    /**
     EXPLANATION:
     Computes the integral of the power spectrum convoluted by pixel's width function
     
     INPUTS:
     vec - a vector of (k, p(k)) pairs
     n - distance between pixels (in number of pixels)
     lya_pixel_width - lya pixel's widht (in km/s)
     sigma_psf - PSF (in km/s)
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Input
     InterpolationMap
     
     FUNCITONS USED:
     NONE
     */
    
    // check that the vector has at least three components
    if (vec.size() < 3){
        std::cout << "Error: in LyaAutoInterpolationMap::IntegratePk : The given vector has too few points" << std::endl;
        return _BAD_DATA_;
    }
    
    double order = double(n);
    double ak_halfs;
    double gaussian_amplitude = 1.0/sqrt(2.0*acos(-1.0))/sigma_psf;
    
    // integrate using Simpson's rule
    size_t n_half = vec.size()/2;
    double h = vec[1].first-vec[0].first;
    
    double sum = 0.0;
    for (size_t i = 1; i <= n_half; i ++){
        ak_halfs = lya_pixel_width*vec[2*i-2].first/2.0;
        sum += vec[2*i-2].second*sin(ak_halfs)/ak_halfs*cos(order*ak_halfs)*gaussian_amplitude*exp(-vec[2*i-2].first*vec[2*i-2].first/2.0/sigma_psf/sigma_psf);
        
        ak_halfs = lya_pixel_width*vec[2*i-1].first/2.0;
        sum += 4.0*vec[2*i-1].second*sin(ak_halfs)/ak_halfs*cos(order*ak_halfs)*gaussian_amplitude*exp(-vec[2*i-1].first*vec[2*i-1].first/2.0/sigma_psf/sigma_psf);
        
        ak_halfs = lya_pixel_width*vec[2*i].first/2.0;
        sum += vec[2*i].second*sin(ak_halfs)/ak_halfs*cos(order*ak_halfs)*gaussian_amplitude*exp(-vec[2*i].first*vec[2*i].first/2.0/sigma_psf/sigma_psf);
    }
    sum *= h/3.0/acos(-1.0);
    
    // add the contribution of 0 < k < k_{0}
    ak_halfs = lya_pixel_width*vec[0].first/2.0;
    sum += vec[0].first*vec[1].second*sin(ak_halfs)/ak_halfs*cos(order*ak_halfs)/acos(-1.0)*gaussian_amplitude*exp(-vec[0].first/2.0*vec[0].first/2.0/2.0/sigma_psf/sigma_psf);
    
    return sum;
}

void LyaAutoInterpolationMap::ProjectLyaAuto(){
    /**
     EXPLANATION:
     Projects the 1D autocorrelation
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaAutoInterpolationMap
     
     FUNCITONS USED:
     NONE
     */
    
    std::map<double,double>::const_iterator it;
    
    // compute forest variables (loop over forest pixels, twice)
    double forest_aux1 = 0.0; //forest_aux: sum of (z-forest_mean_z)**2 of all the pixels in the forest
    double forest_aux2 = 0.0; //forest_aux2: sum of xi*(z-forest_mean_z) of all the pixels in the forest
    double forest_total_weight = 0.0;
    double forest_mean_z = 0.0;
    double forest_mean_xi = 0.0;
    
    for (it = interpolation_map_.begin(); it != interpolation_map_.end(); it ++){
        forest_total_weight += 1.0;
        forest_mean_z += (*it).first;
        forest_mean_xi += (*it).second;
    }
    forest_mean_z /= forest_total_weight;
    forest_mean_xi /= forest_total_weight;
    
    for (it = interpolation_map_.begin(); it != interpolation_map_.end(); it ++){
        forest_aux1 += ((*it).first-forest_mean_z)*((*it).first-forest_mean_z);
        forest_aux2 += (*it).second*((*it).first-forest_mean_z);
    }
    double forest_aux = forest_aux2/forest_aux1;
    
    // project the delta field
    double projected_delta;
    for (it = interpolation_map_.begin(); it != interpolation_map_.end(); it ++){
        projected_delta = (*it).second-forest_mean_xi-forest_aux*((*it).first-forest_mean_z);
        interpolation_map_[(*it).first] = projected_delta;
    }

    
    
}

