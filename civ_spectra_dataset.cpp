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
    std::vector<std::string> fields_hdu1(3);
    fields_hdu1[0] = "LOBS";
    fields_hdu1[1] = "DELTA";
    fields_hdu1[2] = "WEIGHT";
    std::vector<std::string> fields_hdu2(6);
    fields_hdu2[0] = "RA";
    fields_hdu2[1] = "DEC";
    fields_hdu2[2] = "Z";
    fields_hdu2[3] = "PLATE";
    fields_hdu2[4] = "MJD";
    fields_hdu2[5] = "FIBER_D";
    
    // construct fits object
    std::auto_ptr<CCfits::FITS> pInfile;
    
    try{
        pInfile = std::auto_ptr<CCfits::FITS>(new CCfits::FITS(lya_spectra_catalog,CCfits::Read,true));
        
    } catch(CCfits::FITS::CantOpen x) {
        
        throw "Error: in QuasarDataset::Load : Couldn't open catalog file: " + lya_spectra_catalog;
    }
    std::cout << "loading hdu 1\n";
    CCfits::ExtHDU& data_hdu1 = pInfile->extension(1);
    std::cout << "loading hdu 1\n";
    CCfits::ExtHDU& data_hdu2 = pInfile->extension(2);
    
    // number of lines in the file
    long NAXIS2_HDU1 = data_hdu1.axis(1);
    size_t nobj_hdu1 = NAXIS2_HDU1;
    long NAXIS2_HDU2 = data_hdu2.axis(1);
    size_t nobj_hdu2 = NAXIS2_HDU2;

    
    // this will store the information
    std::valarray<int> plate, mjd, fiber;
    std::valarray<double> ra, dec, z, lobs, delta, weight;
    
    // reading data
    std::cout << "reading data from hdu 1\n";
    data_hdu1.column(fields_hdu1[0]).read(lobs,1,nobj_hdu1); // observed wavelength
    data_hdu1.column(fields_hdu1[1]).read(delta,1,nobj_hdu1); // measured overdensity (delta)
    data_hdu1.column(fields_hdu1[2]).read(weight,1,nobj_hdu1); // weight
    std::cout << "reading data from hdu 1\n";
    data_hdu2.column(fields_hdu2[0]).read(ra,1,nobj_hdu2); // ra (in degrees)
    data_hdu2.column(fields_hdu2[1]).read(dec,1,nobj_hdu2); // dec (in degrees)
    data_hdu2.column(fields_hdu2[2]).read(z,1,nobj_hdu2); // z
    data_hdu2.column(fields_hdu2[3]).read(plate,1,nobj_hdu2); // plate
    data_hdu2.column(fields_hdu2[4]).read(mjd,1,nobj_hdu2); // mjd
    data_hdu2.column(fields_hdu2[5]).read(fiber,1,nobj_hdu2); // fiber
    
    // setting size to zero and creating entries in list_ map
    size_ = 0;
    for (int i = 0; i < nobj_hdu2; i++){
        
        if (list_.find(plate[i]) == list_.end()){
            std::vector<LyaSpectrum> v;
            list_[plate[i]] = v;
            num_objects_in_plate_[plate[i]] = 0;
        }
        
    }
    
    // adding objects to list_ and to list_by_plate
    size_t forest_size = nobj_hdu2/nobj_hdu1;
    std::vector<double> lobs_q(forest_size, 0.0), delta_q(forest_size, 0.0), weight_q(forest_size, 0.0);
    for (int i=0;i<nobj_hdu2;i++){
        if ( not (ra[i] == 0.0 and dec[i] == 0.0)){
            
            // create LyaSpectrum
            for (int j=0; j < forest_size; j++){
                lobs_q[j] = lobs[i + j];
                delta_q[j] = delta[i + j];
                weight_q[j] = weight[i + j];
            }
            LyaSpectrum object(ra[i], dec[i], plate[i], fiber[i], mjd[i], z[i], lobs_q, delta_q, weight_q, false);
            
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