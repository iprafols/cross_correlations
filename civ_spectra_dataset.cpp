/**
 civ_spectra_dataset.cpp
 Purpose: This files contains the body for the functions defined in civ_spectra_dataset.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 08/08/2014
 */

#include "civ_spectra_dataset.h"

CIVSpectraDataset::CIVSpectraDataset(const Input& input){
    /**
     EXPLANATION:
     Cosntructs a CIVSpectraDataset instance and loads a catalog of LyaSpectrum objects into it
     
     INPUTS:
     input - object of type cInput
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     CIVSpectraDataset
     Input
     
     FUNCITONS USED:
     NONE
     */
    
    flag_verbose_civ_spectra_dataset_ = input.flag_verbose_civ_spectra_dataset();
    name_ = input.dataset2_name();
    Load(input.dataset2(), "", input.lya_wl());
    
}

void CIVSpectraDataset::Load(const std::string& lya_spectra_catalog, const std::string& lya_spectra_dir, const double& lya_wl){
    /**
     EXPLANATION:
     Loads the object dataset from a catalog file
     
     INPUTS:
     lya_spectra_catalog - a string with the name of the list of fits files containing the spectra
     lya_spectra_dir - a string with the directory where the ly-a spectrum files are stored
     lya_wl - a double with the wavelength of the lyman-alpha line (in Angstroms)
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaSpectrum
     LyaSpectraDataset
     
     FUNCITONS USED:
     NONE
     */
    
    if (flag_verbose_civ_spectra_dataset_ >= 1){
        std::cout << "Loading civ spectra dataset" << std::endl;
    }

    // setting the catalog columns to be read
    std::vector<std::string> fields(9);
    fields[0] = "LOBS";
    fields[1] = "DELTA";
    fields[2] = "WEIGHT";
    fields[3] = "RA";
    fields[4] = "DEC";
    fields[5] = "Z";
    fields[6] = "PLATE";
    fields[7] = "MJD";
    fields[8] = "FIBER_ID";
    
    // construct fits object
    std::auto_ptr<CCfits::FITS> pInfile;
    
    try{
        
        pInfile = std::auto_ptr<CCfits::FITS>(new CCfits::FITS(civ_filename,CCfits::Read,1,true,fields));
        
    } catch(CCfits::FITS::CantOpen x) {
        
        throw "Error: in QuasarDataset::Load : Couldn't open catalog file: " + civ_filename;
    }
    CCfits::ExtHDU& data = pInfile->extension(1);
    
    // number of lines in the file
    long NAXIS2 = data.axis(1);
    size_t nobj = NAXIS2;
    
    // this will store the information
    std::valarray<int> plate, mjd, fiber;
    std::valarray<double> ra, dec, z;
    std::vector<std::valarray<double> > lobs, delta, weight;
    
    // reading data
    data.column(fields[0]).readArrays(lobs, 1, nobj); // observed wavelength
    data.column(fields[1]).readArrays(delta, 1, nobj); // measured overdensities
    data.column(fields[2]).readArrays(weight, 1, nobj); // weights
    data.column(fields[3]).read(ra, 1, nobj); // ra
    data.column(fields[4]).read(dec, 1, nobj); // dec
    data.column(fields[5]).read(z, 1, nobj); // z
    data.column(fields[6]).read(plate, 1, nobj); // plate number
    data.column(fields[7]).read(mjd, 1, nobj); // mjd
    data.column(fields[8]).read(fiber, 1, nobj); // fiber number
    
    // setting size to zero and creating entries in list_ map
    size_ = 0;
    for (int i = 0; i < nobj; i++){
        
        if (list_.find(plate[i]) == list_.end()){
            std::vector<LyaSpectrum> v;
            list_[plate[i]] = v;
            num_objects_in_plate_[plate[i]] = 0;
        }
        
    }
    
    // adding objects to list_ and to list_by_plate
    for (int i = 0; i < nobj; i++){
        if ( not (ra[i] == 0.0 and dec[i] == 0.0)){
            
            // create LyaSpectrum
            LyaSpectrum object(ra[i], dec[i], plate[i], fiber[i], mjd[i], z[i], lobs[i], delta[i], weight[i], false);
            
            // adding object to list_
            (*list_.find(plate[i])).second.push_back(object);
            
            // updating size_
            size_ ++;
            
            // updating number_of_objects_in_plate
            (*num_objects_in_plate_.find(plate[i])).second ++;
            
            if (flag_verbose_civ_spectra_dataset_ >= 3 or (flag_verbose_civ_spectra_dataset_ >= 2 and size_ == size_/1000*1000)){
                std::cout << "Loaded " << size_ << " civ spectra" << std::endl;
            }
            
            
        }
    }
    
    if (flag_verbose_civ_spectra_dataset_ >= 1){
        std::cout << "Loaded " << size_ << " civ spectra" << std::endl;
    }

}