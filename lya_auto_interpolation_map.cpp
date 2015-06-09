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
    
    double lya_pixel_width = input.lya_pixel_width();
    
    double aux_z = 0.0;
    double value;
    std::vector<std::pair<double, double> > k_pk;
    
    std::ifstream file (input.lya_auto_correlation().c_str());
    if (file.is_open()){
        while (getline(file,line)){
            std::vector<std::string> cols = Split(line," ");
            
            if (double(atof(cols[0])) != aux_z){
                // integrate pk and add autocorrelation to interpolation map
                if (not k_pk.empty()){
                    value = IntegratePk(k_pk, n, lya_pixel_width);
                    if (value != _BAD_DATA_){
                        interpolation_map_[aux_z] = value;
                    }
                    k_pk.clear();
                }
                // update redshift value
                aux_z = double(atof(cols[0]));
            
            }
            
            // add pair of (k, p(k)) to vector
            k_pk.push_back(std::pair(double(atof(cols[1])), double(atof(cols[2]))));
            
        }
    }
    else{
        std::cout << "Error: in LyaAutoInterpolationMap::LyaAutoInterpolationMap : Could not read file" << std::endl;
        std::cout << input.lya_auto_correlation() << std::endl;
    }
}

double LyaAutoInterpolationMap::IntegratePk(const std::vector<std::pair<double, double> >& vec, const int& n, const double& lya_pixel_width){
    /**
     EXPLANATION:
     Computes the integral of the power spectrum convoluted by pixel's width function
     
     INPUTS:
     vec - a vector of (k, p(k)) pairs
     n - distance between pixels (in number of pixels)
     lya_pixel_width - lya pixel's widht (in km/s)
     
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
    
    double order = 2.0*double(n)+1.0;
    double ak_halfs;
    
    // integrate using Simpson's rule
    size_t n_half = vec.size()/2;
    double h = vec[1].first-vec[0].first;
    
    double complex sum = 0.0;
    for (size_t i = 1; i <= n_half; i ++){
        ak_halfs = lya_pixel_width*vec[2*i-2].first/2.0;
        sum += vec[2*i-2].second*sin(ak_halfs)/ak_halfs*cos(order*ak_halfs);
        
        ak_halfs = lya_pixel_width*vec[2*i-1].first/2.0;
        sum += 4.0*vec[2*i-1].second*sin(ak_halfs)/ak_halfs*cos(order*ak_halfs);
        
        ak_halfs = lya_pixel_width*vec[2*i].first/2.0;
        sum += vec[2*i].second*sin(ak_halfs)/ak_halfs*cos(order*ak_halfs);
    }
    sum *= h/6.0/acos(-1.0);
    
    return sum;
}




