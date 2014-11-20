/**
 function_compute_plate_neighbours.cpp
 Purpose: Create a list of all the used plates and their neighbouring plates
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 06/17/2014
 */

// libraries used
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <fstream>
#include <map>
#include <string>
#include <time.h>
#include <vector>
////////

// classes used
#include "input.h"
#include "plate.h"
#include "plate_neighbours.h"
////////

// functions used
////////

#include "typedefs.h"

void ComputePlateNeighbours(const Input& input){
    /**
     EXPLANATION:
     Create a list of all the used plates and their neighbouring plates
     
     INPUTS:
     input - a Input instance
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Input
     Plate
     PlateNeighbours
     
     FUNCITONS USED:
     NONE
     */
    
    // load time control variables
    time_t start_time,end_time;
    time(&start_time);
    
    std::cout << "Computing plate neighbours list" << std::endl;
    // load variables
    const std::string kLyaSpectraDir = input.lya_spectra_dir();
    const double kNeighboursMaxDistance = input.neighbours_max_distance();
    // load empty plate map
    PlatesMapSimple<Plate>::map plate_list;
    
    std::cout << "Reading spectra catalog" << std::endl;
    std::ifstream catalog(input.dataset2().c_str());
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
        std::exit;
    }
    
    // normalize contents
    std::cout << "Normalizing mean right ascension and declination for the different plates" << std::endl;
    for (PlatesMapSimple<Plate>::map::iterator it = plate_list.begin(); it != plate_list.end(); it ++){
        (*it).second.Normalize();
    }
    
    // load empty neighbours map
    std::cout << "Creating neighbours map" << std::endl;
    PlateNeighbours neighbours_list;
    
    // every plate is considered a neighbour of itself
    for (PlatesMapSimple<Plate>::map::iterator it = plate_list.begin(); it != plate_list.end(); it ++){
        neighbours_list.AddPlate((*it).second);
    }
    
    // look for neighbours
    std::cout << "Looking for neighbours" << std::endl;
    for (PlatesMapSimple<Plate>::map::iterator it = plate_list.begin(); it != plate_list.end(); it ++){ // loop over plates -> plate1
        for (PlatesMapSimple<Plate>::map::iterator it2 = it; it2 != plate_list.end(); it2 ++){ // loop over plates -> plate2
            if (it == it2){
                continue;
            }
            
            // if they are neigbours
            if ((*it).second.IsNeighbour((*it2).second,kNeighboursMaxDistance)){
                neighbours_list.AddNeighbours((*it).second,(*it2).second);
            }
            
        }
    }
    
    // save results
    neighbours_list.Save(input.plate_neighbours());
    
    // display time required to run the program
    std::cout << "Plate neightbours list computed" << std::endl;
    time(&end_time);
    double time_spent = difftime(end_time, start_time);
    if (time_spent < 60.0){
        std::cout << "It took " << time_spent << " seconds to compute the neighbours list" << std::endl;
    }
    else if (time_spent < 3600.0){
        std::cout << "It took " << time_spent/60.0 << " minutes to compute the neighbours list" << std::endl;
    }
    else{
        std::cout << "It took " << time_spent/3600.0 << " hours to compute the neighbours list" << std::endl;
    }
    
}