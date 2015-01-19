/**
 quasar_dataset.cpp
 Purpose: This files contains the body for the functions defined in quasar_dataset.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 11/25/2014
 */

#include "quasar_dataset.h"


QuasarDataset::QuasarDataset(const Input& input){
    /**
     EXPLANATION:
     Cosntructs a QuasarDataset instance and loads a catalog of AstroObjects into it
     
     INPUTS:
     input - object of type Input
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     QuasarDataset
     Input
     
     FUNCITONS USED:
     NONE
     */

    flag_verbose_quasar_dataset_ = input.flag_verbose_quasar_dataset();
    name_ = input.dataset1_name();
    Load(input.z_min(), input.z_max(), input.dataset1());
    
}

void QuasarDataset::Load(const double& z_min, const double& z_max, const std::string& dataset1){
    /**
     EXPLANATION:
     Loads the object dataset from a catalog file
     
     INPUTS:
     z_min - minimum redshift to accept AstroObject into the dataset
     z_max - maximum redshift to accept AstroObject into the dataset
     dataset1 - name of the catalog file
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     AstroObjectDataset
     
     FUNCITONS USED:
     NONE
     */
    if (flag_verbose_quasar_dataset_ >= 1){
        std::cout << "Loading quasar dataset" << std::endl;
    }

    // setting the catalog columns to be read
    std::vector<std::string> fields(8);
    fields[0] = "PLATE";
    fields[1] = "MJD";
    fields[2] = "FIBERID";
    fields[3] = "RA";
    fields[4] = "DEC";
    fields[5] = "Z_VI";
    fields[6] = "PSFMAG";
    //fields[7] = "SPECPRIMARY";
    
    // construct fits object
    std::auto_ptr<CCfits::FITS> pInfile; 
    
    try{
        
        pInfile = std::auto_ptr<CCfits::FITS>(new CCfits::FITS(dataset1,CCfits::Read,1,true,fields));
        
    } catch(CCfits::FITS::CantOpen x) {
        
        throw "Error: in QuasarDataset::Load : Couldn't open catalog file: " + dataset1;
    }
    CCfits::ExtHDU& data = pInfile->extension(1);
    
    // number of lines in the file
    long NAXIS2 = data.axis(1);
    size_t nobj = NAXIS2;
    
    // this will store the information
    std::valarray<int> plate,mjd,fiber;
    std::valarray<double> ra,dec,z;
    std::vector<valarray<double> > psf_mag;
    std::valarray<bool> boss_target1_flag,specprimary_flag;
    
    // reading data
    data.column(fields[0]).read(plate,1,nobj); // plate
    data.column(fields[1]).read(mjd,1,nobj); // mjd
    data.column(fields[2]).read(fiber,1,nobj); // fiber
    data.column(fields[3]).read(ra,1,nobj); // ra
    data.column(fields[4]).read(dec,1,nobj); // dec
    data.column(fields[5]).read(z,1,nobj); // z
    data.column(fields[6]).readArrays(psf_mag,1,nobj); // magnitudes
    //data.column(fields[7]).read(specprimary_flag,1,nobj); // specprimary
    
    // setting size to zero and creating entries in list_ map
    size_ = 0;
    for (int i=0;i<nobj;i++){
        // this is to avoid repeated objects
        //if (specprimary_flag[i] > 0){
            
            if (list_.find(plate[i]) == list_.end()){
                std::vector<AstroObject> v;
                list_[plate[i]] = v;
                num_objects_in_plate_[plate[i]] = 0;
            }
            
        //}
    }
    
    // adding objects to list_ and to list_by_plate
    for (int i=0;i<nobj;i++){
        // this is to avoid repeated objects
        //if ((specprimary_flag[i] > 0) and (z_min <= z[i] and z[i] <= z_max) and not (ra[i] == 0.0 and dec[i] == 0.0)){
        if ((z_min <= z[i] and z[i] <= z_max) and not (ra[i] == 0.0 and dec[i] == 0.0)){
            
            // create AstroObject
            AstroObject object(ra[i], dec[i], plate[i], fiber[i], mjd[i], z[i], false);
            
            // adding object to list_
            (*list_.find(plate[i])).second.push_back(object);
            
            // updating size_
            size_ ++;
            
            // updating number_of_objects_in_plate
            (*num_objects_in_plate_.find(plate[i])).second ++;
            
            if (flag_verbose_quasar_dataset_ >= 3 or (flag_verbose_quasar_dataset_ >= 2 and size_ == size_/1000*1000)){
                std::cout << "Loaded " << size_ << " quasars" << std::endl;
            }

            
        }
    }
    if (flag_verbose_quasar_dataset_ >= 1){
        std::cout << "Loaded " << size_ << " quasars" << std::endl;
    }
}
