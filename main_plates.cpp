/**
 main_plates_neighbours.cpp
 Purpose: Create a list of all the used plates and their neighbouring plates
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 06/17/2014
 */

// libraries used
#include <iostream>
#include <iterator>
#include <fstream>
#include <map>
#include <string>
#include <time.h>
#include <vector>
////////

// classes used
#include "global_variables.h"
#include "plate.h"
#include "plate_neighbours.h"
////////

// functions used
////////

#include "typedefs.h"

int main(){
    /**
     EXPLANATION:
     Create a list of all the used plates and their neighbouring plates
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     cGlobalVariables
     cPlate
     cPlateNeighbours
     
     FUNCITONS USED:
     NONE
     */
    
    // load time control variables
    time_t start_time,end_time;
    time(&start_time);
    
    std::cout << "Initializing variables" << std::endl;
    // load global variables
    const GlobalVariables kGlobalVariables;
    const std::string kLyaSpectraDir = kGlobalVariables.lya_spectra_dir();
    const double kNeighboursMaxDistance = kGlobalVariables.neighbours_max_distance();
    // load empty plate map
    PlatesMapSimple<Plate>::map plate_list;
    
    std::cout << "Reading spectra catalog" << std::endl;
    std::ifstream catalog(kGlobalVariables.lya_spectra_catalog().c_str());
    // open catalog
    int number_of_read_spectra = 0;
    if (catalog.is_open()){
        std::string file("");        
        while (getline(catalog,file)){
            // load plate information
            Plate plate(kLyaSpectraDir+file);
            PlatesMapSimple<Plate>::map::iterator it = plate_list.find(plate.plate_number());
            // if plate is still not in the list, include it
            if (it == plate_list.end()){
                plate_list[plate.plate_number()] = plate;
            }
            // otherwise update the mean right ascension and declination
            else{
                (*it).second.UpdateRADECValues(plate);
            }
            
            number_of_read_spectra ++;
            if (number_of_read_spectra/100*100 == number_of_read_spectra){
                std::cout << "read " << number_of_read_spectra << " spectra" << std::endl;
            }
            
        }
        catalog.close();
    }
    else{
        std::cout << "Error: could not read spectra catalog" << std::endl;
        return 0;
    }
    
    // show contents
    /*for (PlatesMapSimple<Plate>::map::iterator it = plate_list.begin(); it != plate_list.end(); it ++){
        std::cout << "plate " << (*it).first << " contains " << (*it).second.number_of_objects() << " objects; ra = " << (*it).second.angle().ra() << "; dec = " << (*it).second.angle().dec() << std::endl; 
    }*/

    // normalize contents
    std::cout << "Normalizing mean right ascension and declination for the different plates" << std::endl;
    for (PlatesMapSimple<Plate>::map::iterator it = plate_list.begin(); it != plate_list.end(); it ++){
        (*it).second.Normalize();
    }
    
    // show contents
    /*for (PlatesMapSimple<Plate>::map::iterator it = plate_list.begin(); it != plate_list.end(); it ++){
        std::cout << "plate " << (*it).first << " : ra = " << (*it).second.angle().ra() << "; dec = " << (*it).second.angle().dec() << std::endl; 
    }*/
    
    // load empty neighbours map
    std::cout << "Creating neighbours map" << std::endl;
    PlateNeighbours neighbours_list;
    
    // every plate is considered a neighbour of itself
    for (PlatesMapSimple<Plate>::map::iterator it = plate_list.begin(); it != plate_list.end(); it ++){
        //std::cout << "inserting plate " << (*it).first;
        neighbours_list.AddPlate((*it).second);
        //std::cout << "; neigbours list size is " << neighbours_list.plates().size() << std::endl;
    }
    
    // show contents
    /*PlatesMapVector<Plate>::map neighbours_plates = neighbours_list.plates();
    for (PlatesMapVector<Plate>::map::iterator it = neighbours_plates.begin(); it != neighbours_plates.end(); it ++){
        std::cout << "found plate " << (*it).first << std::endl; 
    }*/

    
    // look for neighbours
    std::cout << "Looking for neighbours" << std::endl;
    for (PlatesMapSimple<Plate>::map::iterator it = plate_list.begin(); it != plate_list.end(); it ++){ // loop over plates -> plate1
        for (PlatesMapSimple<Plate>::map::iterator it2 = it; it2 != plate_list.end(); it2 ++){ // loop over plates -> plate2
            if (it == it2){
                continue;
            }
            
            //std::cout << "checking plates " << (*it).first << " and " << (*it2).first << std::endl;
            // if they are neigbours
            if ((*it).second.IsNeighbour((*it2).second,kNeighboursMaxDistance)){
                neighbours_list.AddNeighbours((*it).second,(*it2).second);
            }
            
        }
    }
    
    // print results
    /*PlatesMapVector<Plate>::map neighbours_plates2 = neighbours_list.plates();
    for (PlatesMapVector<Plate>::map::iterator it = neighbours_plates2.begin(); it != neighbours_plates2.end(); it ++){
        
        std::cout << "plate " << (*it).first << " has the following neighbours:" << std::endl; 
        std::vector<int> v = neighbours_list.GetNeighboursList((*it).first);
        
        for (size_t i = 0;i < v.size(); i++){
            std::cout << v[i] << " ";
        }
        std::cout << std::endl;
    }*/
    
    // save results
    neighbours_list.Save(kGlobalVariables.plate_neighbours());
    
    // display time required to run the program
    std::cout << "End of program" << std::endl;
    time(&end_time);
    double time_spent = difftime(end_time, start_time);
    std::cout << "The program lasted " << time_spent << " seconds. This corresponds to " << time_spent/60.0 << " minutes or " << time_spent/3600.0 << " hours" << std::endl;

    return 0;
    
    


}