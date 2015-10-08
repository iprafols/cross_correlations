/**
 lya_spectrum.h
 Purpose: This file defines the class LyaSpectrum. This class contains the variables necessary to store a lyman-alpha spectrum. This class is a specialization of the cAstroObject class
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 09/18/2014
 
 */

#ifndef _LyaSpectrum_h
#define _LyaSpectrum_h



// libraries needed
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <valarray>
#include <vector>


#include <CCfits>
////////

// classes needed
#include "astro_object.h"
#include "interpolation_map.h"
#include "lya_pixel.h"
#include "sphere_point.h"
////////

// functions needed
////////

#include "defines.h"

class LyaSpectrum: public AstroObject{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs a "bad data" LyaSpectrum
    LyaSpectrum(double bad_data);

    // constructs object and initializes its variables
    LyaSpectrum(const std::string& filename, const double& lya_wl, const bool radians = true);
    LyaSpectrum(const double& ra, const double& dec, const int& plate, const int& fiber, const int& mjd, const double& z, const std::valarray<double>& lobs, std::valarray<double>& delta, std::valarray<double>& weight, const double& lya_wl, const bool radians = true);
    
    // -------------------------------------------------------------
    // access methods

    // access method for spectrum_
    std::vector<LyaPixel> spectrum() const {return spectrum_;}
    LyaPixel spectrum(size_t i) const;
    
    // -------------------------------------------------------------
    // other methods
    
    // set the distance to object
    void SetDistance(const InterpolationMap& redshift_distance_map);
    
    // return the size of spectrum_;
    size_t SpectrumSize() const {return spectrum_.size();}
    
    
    
private:
    // lyman-alpha spectrum
    std::vector<LyaPixel> spectrum_;
    
};    


#endif
