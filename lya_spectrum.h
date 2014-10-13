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
#include <cstring>
#include <string>
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


class LyaSpectrum: public AstroObject{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs object and initializes its variables
    LyaSpectrum(const std::string& filename, const double& lya_wl, const bool radians = true);
    
    // -------------------------------------------------------------
    // access methods

    // access method for spectrum_
    std::vector<LyaPixel> spectrum() const {return spectrum_;}
    LyaPixel spectrum(size_t i) const {return spectrum_[i];}
    
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
