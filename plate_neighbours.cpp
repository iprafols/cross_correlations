/**
 plate_neighbours.cpp
 Purpose: This files contains the body for the functions defined in plate_neighbours.cpp.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 06/26/2014
 */

#include "plate_neighbours.h"

PlateNeighbours::PlateNeighbours(const Input& input){
    /**
     EXPLANATION:
     Cosntructs a PlateNeighbours instance and initializes all its variables
     
     INPUTS:
     input - a Input instance
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     PlateNeighbours
     
     FUNCITONS USED:
     ComputePlateNeighbours
     */
    
    std::ifstream plates_file(input.plate_neighbours().c_str());
    if (plates_file.is_open()){
        // read plate neighbours
        ReadPlateNeighbours(plates_file);
        
        plates_file.close();
    }
    else{
        std::cout << "Error: Unable to open file: " << std::endl << input.plate_neighbours() << std::endl;
        std::cout << "This may be because the list is not created and the corresponding flag was not set." << std::endl << "Do you want to compute the plate neighbours list? (this may take a while) [y/n]" << std::endl;
        std::string ans = "";
        while (ans == ""){
            std::cin >> ans;
            if ((ans != "y") and (ans != "n")){
                ans = "";
                std::cout << "Please answer 'y' or 'n'" << std::endl;
            }
        }
        if (ans == "y"){
            ComputePlateNeighbours(input);
            std::ifstream plates_file2(input.plate_neighbours().c_str());
            if (plates_file2.is_open()){
                // read plate neighbours
                ReadPlateNeighbours(plates_file2);
                
                plates_file2.close();
            }
            else{
                std::cout << "Error: Unable to open file: " << std::endl << input.plate_neighbours() << std::endl;
                std::exit;
            }
        }
        else{
            std::exit;
        }
    }
    
}

void PlateNeighbours::AddNeighbours(const Plate& plate1, const Plate& plate2){
    /**
     EXPLANATION:
     The two given plates are considered neighbours, add them to the respective neighbours list
     
     INPUTS:
     plate1 - a Plate object 
     plate2 - another Plate object
     
     OUTPUTS:
     plates - a vector containg the list of used plates
     
     CLASSES USED:
     Plate
     PLateNeighbours
     
     FUNCITONS USED:
     ToStr
     */
    
    // set second plate as neighbour of the first
    std::map<int,std::vector<Plate> >::iterator it = plates_.find(plate1.plate_number());
    (*it).second.push_back(plate2);
    
    // set first plate as neighbour of the second
    it = plates_.find(plate2.plate_number());
    (*it).second.push_back(plate1);
    
}

void PlateNeighbours::AddPlate(const Plate& plate){
    /**
     EXPLANATION:
     If the given plate does not have an entry, then make a new entry
     
     INPUTS:
     NONE
     
     OUTPUTS:
     plates - a vector containg the list of used plates
     
     CLASSES USED:
     PLateNeighbours
     
     FUNCITONS USED:
     ToStr
     */
        
    // check whether or not the entry is already existing
    if (plates_.find(plate.plate_number()) != plates_.end()){
        std::cout << "The given plate is already included" << std::endl;
        return;
    }
    
    std::vector<Plate> v(1,plate);
    plates_[plate.plate_number()] = v;

}

std::vector<int> PlateNeighbours::GetNeighboursList(int plate) const{
    /**
     EXPLANATION:
     Returns a vector conatining the plate number of all the neighbours of a specified plate
     
     INPUTS:
     plate - plate number of the plate from which the neighbours are recovered
     
     OUTPUTS:
     neighbours - a vector containg the list of used plates
     
     CLASSES USED:
     PLateNeighbours
     
     FUNCITONS USED:
     NONE
     */
    
    if (plate == _NORM_){
        std::vector<int> neighbours;
        return neighbours;
    }
    else{
        // locate plate
        PlatesMapVector<Plate>::map::const_iterator map_it = plates_.find(plate);
        
        // if not found, return empty vector
        if (map_it == plates_.end()){
            std::vector<int> neighbours;
            return neighbours;
        }
        
        // create vector of neighbours with the correct size, fill it with zeros
        std::vector<int> neighbours((*map_it).second.size(),0);
        
        std::vector<int>::iterator vec_it1;
        std::vector<Plate>::const_iterator vec_it2;
        for (vec_it1 = neighbours.begin(), vec_it2 = (*map_it).second.begin(); vec_it1 != neighbours.end() and vec_it2 != (*map_it).second.end(); vec_it1 ++, vec_it2 ++){
            
            (*vec_it1) = (*vec_it2).plate_number();
        }
        
        return neighbours;
    }
    
}

std::vector<int> PlateNeighbours::GetPlatesList() const {
    /**
     EXPLANATION:
     Returns a vector conatining the plate number of all the used plates
     
     INPUTS:
     NONE
     
     OUTPUTS:
     plates - a vector containg the list of used plates
     
     CLASSES USED:
     PLateNeighbours
     
     FUNCITONS USED:
     NONE
     */
    
    std::vector<int> plates_list;
    plates_list.reserve(plates_.size());
    
    for (PlatesMapVector<Plate>::map::const_iterator it = plates_.begin(); it != plates_.end(); it ++){
        plates_list.push_back((*it).first);
    }
    
    return plates_list;
}

void PlateNeighbours::ReadPlateNeighbours(std::ifstream& plates_file){
    /**
     EXPLANATION:
     Reads the plate neighbours list
     
     INPUTS:
     plates_file - a ifstream instance with the file containing the neighbours list
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     PlateNeighbours
     
     FUNCITONS USED:
     NONE
     */
    std::string line;
    while (getline(plates_file,line)){
        if (line[0] != '#'){
            std::vector<std::string> cols = Split(line," ");
            
            // store plate number
            int plate_number = atoi(cols[0].c_str());
            
            // store neighbours in a vector
            std::vector<Plate> neighbours;
            neighbours.resize(cols.size()-1);
            
            for (size_t i = 0,j = 1; i < neighbours.size() and j < cols.size();i ++, j++){
                neighbours[i].set_plate_number(atoi(cols[j].c_str()));
            }
            
            // add plate
            if (plates_.find(plate_number) == plates_.end()){
                plates_[plate_number] = neighbours;
            }
            else{
                std::cout << "Attempted to create an entry for plate " << plate_number << " when the entry already exists. Check plate neighbours file for repeated lines" << std::endl;
            }
        }
        
    }
}

void PlateNeighbours::Save(const std::string& filename){
    /**
     EXPLANATION:
     Writes a list of plates and their neighbours into the specified file. Format is "plate neighbour1 neighbour2 ..."
     
     INPUTS:
     filename - name of the file in which to write the list
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     PLateNeighbours
     
     FUNCITONS USED:
     NONE
     */
    
    std::cout << "Saving list of plates and neighbours" << std::endl;
    
    std::ofstream file(filename.c_str(),std::ofstream::trunc); 
    if (file.is_open()){
        
        // write header in file
        file << "# plate neighbour1 neighbour2 ..." << std::endl;        
        
        // make sure the map iterator points at the beginnning of the map
        PlatesMapVector<Plate>::map::iterator plates_it = plates_.begin();
        
        while (plates_it != plates_.end()){
            // write plate number
            int plate = (*plates_it).first;
            file << (*plates_it).first << " ";
            
            
            for (std::vector<Plate>::iterator it = (*plates_it).second.begin(); it != (*plates_it).second.end(); it ++){ // loop over neighbours
                
                // write neighbour's plate number
                file << (*it).plate_number() << " ";
            }
            
            file << std::endl;
            
            plates_it ++;
        }
        
        file.close();
    }
    else{
        std::cout << "Unable to open file:" << std::endl << filename << std::endl;
    }
    
}
