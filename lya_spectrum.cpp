/**
 lya_spectrum.cpp
 Purpose: This files contains the body for the functions defined in lya_spectrum.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 09/18/2014
 */

#include "lya_spectrum.h"

LyaSpectrum::LyaSpectrum(double bad_data){
    /**
     EXPLANATION:
     Cosntructs a bad_data AstroObject instance
     
     INPUTS:
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaSpectrum
     LyaSpectrum
     
     FUNCITONS USED:
     NONE
     */
    if (bad_data != _BAD_DATA_){
        std::cout << "Error while initializing a LyaSpectrum 'bad data' instance" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    z_ = _BAD_DATA_;
    dist_ = _BAD_DATA_;
}

LyaSpectrum::LyaSpectrum(const std::string& filename, const double& lya_wl, const bool radians){
    /**
     EXPLANATION:
     Cosntructs a LyaSpectrum instance
     
     INPUTS:
     filename - a string containing the spectrum's fits file name
     lya_wl - restframe lyman alpha wavelength
     radians - a boolean specifying if angles are given in radians (true) or deg (false)
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaSpectrum
     LyaPixel
     SpherePoint
     
     FUNCITONS USED:
     NONE
     */

    // setting the catalog columns to be read
    std::vector<std::string> fields(3);
    fields[0] = "LOGLAM";
    fields[1] = "FOREST";
    fields[2] = "WEIGHT";
    

    // construct fits object
    std::auto_ptr<CCfits::FITS> pInfile; 
    
    try{
        
        pInfile = std::auto_ptr<CCfits::FITS>(new CCfits::FITS(filename,CCfits::Read,1,true,fields));
        
    } catch(CCfits::FITS::CantOpen x) {
        
        throw "Error: in LyaSpectrum::LyaSpectrum : Couldn't open catalog file: " + filename;
    }
    CCfits::ExtHDU& data = pInfile->extension(1);
    
    // number of lines in the file
    long NAXIS2 = data.axis(1);
    size_t nobj = NAXIS2;
    spectrum_.reserve(nobj);
    
    // this will store the information
    std::valarray<double> loglam, forest, weight;
    
    // reading data
    data.column(fields[0]).read(loglam, 1, nobj); // logarithm of the wavelength value
    data.column(fields[1]).read(forest, 1, nobj); // normalized flux in the ly-alpha forest
    data.column(fields[2]).read(weight, 1, nobj); // weight
    
    for (int i=0;i<nobj;i++){
        // create LyaPixel
        LyaPixel object(loglam[i], lya_wl, forest[i], weight[i]);
            
        // adding object to spectrum_
        spectrum_.push_back(object);
    }
    
    // reading header
    data.readAllKeys();
    
    // extract right ascension
    std::string sra("RA");
    double ra;
    data.keyWord(sra).value(ra);
     
    // extract declination
    std::string sdec("DEC");
    double dec;
    data.keyWord(sdec).value(dec);
     
    // set angular position
    if (not radians){
        ra *= acos(-1)/180.0;
        dec *= acos(-1)/180.0;
    }
    SpherePoint angle(ra, dec);
    angle_ = angle;
    
    // extract redshift
    std::string sz("Z");
    data.keyWord(sdec).value(z_);
    
    // extract plate, mjd and fiber numbers
    std::string spmf("PMF");
    std::string pmf;
    data.keyWord(spmf).value(pmf);
    
    plate_ = atoi(strtok((char*)pmf.c_str(),"-"));
    mjd_ = atoi(strtok((char*)pmf.c_str(),"-"));
    fiber_ = atoi(strtok((char*)pmf.c_str(),"-"));

}

LyaSpectrum::LyaSpectrum(const double& ra, const double& dec, const int& plate, const int& fiber, const int& mjd, const double& z, const std::vector<double>& lobs, std::vector<double>& delta, std::vector<double>& weight, const double& lya_wl, const bool radians){
    /**
     EXPLANATION:
     Cosntructs a LyaSpectrum instance
     
     INPUTS:
     ra - right ascention
     dec - declination
     plate - plate number
     fiber - fiber number
     mjd - modified julian day of observation
     z - quasar redshift
     lobs - observed wavelengths
     delta - measured overdensities (deltas)
     weight - weights
     lya_wl - restframe lyman alpha wavelength
     radians - a boolean specifying if angles are given in radians (true) or deg (false)
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaSpectrum
     LyaPixel
     SpherePoint
     
     FUNCITONS USED:
     NONE
     */
    
    // set angular position
    ra_local = ra;
    dec_local = dec;
    if (not radians){
        ra_local *= acos(-1)/180.0;
        dec_local*= acos(-1)/180.0;
    }
    SpherePoint angle(ra_local, dec_local);
    angle_ = angle;
    plate_ = plate;
    fiber_ = fiber;
    mjdi_ = mjd;
    z_ = z;
    
    for (int i=0; i < lobs.size(); i++){
        // create LyaPixel
        LyaPixel object(lobs[i], lya_wl, delta[i], weight[i], false);
        
        // adding object to spectrum_
        spectrum_.push_back(object);
    }
    
}

LyaPixel LyaSpectrum::spectrum(size_t i) const {
    /**
     EXPLANATION:
     Access function for spectrum_
     
     INPUTS:
     i - index of the selected spectrum_ element
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaSpectrum
     LyaPixel
     
     FUNCITONS USED:
     NONE
     */
    if (i < spectrum_.size()){
        return spectrum_[i];
    }
    else{
        LyaPixel lp(_BAD_DATA_);
        return lp;
    }
    
}

void LyaSpectrum::SetDistance(const InterpolationMap& redshift_distance_map){
    /**
     EXPLANATION:
     Sets the distance to object
     
     INPUTS:
     redshif_distance_map - a InterpolationMap instance with the redshift-distance relation
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     
     FUNCITONS USED:
     NONE
     */
    
    
    dist_ = NAN;
    
    for (size_t i = 0; i < spectrum_.size(); i++){
        spectrum_[i].SetDistance(redshift_distance_map);
    }
    
}