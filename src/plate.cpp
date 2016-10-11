/**
 plate.cpp
 Purpose: This files contains the body for the functions defined in plate.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 06/26/2014
 */

#include "plate.h"

Plate::Plate(const int& plate_number){
    /**
     EXPLANATION:
     Cosntructs a Plate instance and initializes all its variables assigning zero values
     
     INPUTS:
     plate_number - an integer containing the plate number
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Plate
     
     FUNCITONS USED:
     NONE
     */
    
    // set plate number
    plate_number_ = plate_number;
    
    // set number of averaged objects to 1
    number_of_objects_ = 0.0;
    
    // set position angle
    SpherePoint angle(0.0, 0.0);
    angle_ = angle;
}

Plate::Plate(const std::string& filename){
    /**
     EXPLANATION:
     Cosntructs a Plate instance and initializes all its variables
     
     INPUTS:
     filename - string containing the full path of the fits file containing the object
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Plate
     
     FUNCITONS USED:
     NONE
     */
    
    // constructs fits object
    std::auto_ptr<CCfits::FITS> pInfile;//std::unique_ptr<FITS> pInfile;
    
    try{
        
        pInfile = std::auto_ptr<CCfits::FITS>(new CCfits::FITS(filename, CCfits::Read));//std::unique_ptr<FITS>(new FITS(filename, Read));
        
    } catch(CCfits::FITS::CantOpen x) {
        
        throw "Error : In Plate::Plate : " + filename + " failed";
    }

    // define a reference for clarity
    CCfits::ExtHDU& data = (*pInfile).extension(1);
    
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
    
    // extract plate number
    std::string spmf("PMF");
    std::string pmf;
    data.keyWord(spmf).value(pmf);
    
    // set plate number
    plate_number_ = atoi(strtok((char*)pmf.c_str(),"-"));
    
    // set number of averaged objects to 1
    number_of_objects_ = 1.0;
    
    // set position angle
    SpherePoint angle(ra, dec);
    angle_ = angle;
}

bool Plate::IsNeighbour(const Plate& plate,const double& neighbours_max_distance){
    /**
     EXPLANATION:
     Returns true if the given plate is a neighbour plate and false otherwise
     
     INPUTS:
     plate - plate to check neighbourhood with
     dist - maximum angular separation for a couple of plates to be considered neighbours (in radians)
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Plate
     
     FUNCITONS USED:
     NONE
     */
    
    return angle_.AngularDistance(plate.angle()) <= neighbours_max_distance;
}

void Plate::Normalize(){
    /**
     EXPLANATION:
     Normalizes angle_ by dividing it by number_of_objects_. Then it sets number_of_objects_ to 0 and updates sin_dec_ and cos_dec_
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Plate
     
     FUNCITONS USED:
     NONE
     */
    
    // check that the plate has not already been normalized
    if (number_of_objects_ == _NORM_){
        std::cout << "Warning: in Plate::Normalize : Plate has already been normalized. Ignoring..." << std::endl;
        return;
    }
    angle_ /= number_of_objects_;
    number_of_objects_ = _NORM_;
}

void Plate::UpdateRADECValues(const Plate& plate){
    /**
     EXPLANATION:
     Adds the ra_ and dec_ values of another cPlate object and increases number_of_objects_ by the corresponding value
     
     INPUTS:
     plate - plate instace with the same plate number
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Plate
     
     FUNCITONS USED:
     NONE
     */
    
    // check that the plate numbers are indeed the same
    if (plate_number_ != plate.plate_number()){
        std::cout << "Error: in Plate::UpdateRADECValues : plates numbers are not the same. skipping addition of ra and dec values" << std::endl;
        return;
    }
    angle_ += plate.angle();
    number_of_objects_ += plate.number_of_objects();
}

void Plate::UpdateRADECValues(const SpherePoint& angle){
    /**
     EXPLANATION:
     Adds the specified ra and dec values and increases number_of_objects_ by 1
     
     INPUTS:
     angle - SpherePoint object containing the ra and dec values to add
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Plate
     
     FUNCITONS USED:
     NONE
     */

    angle_ += angle;
    number_of_objects_ += 1.0;
}

