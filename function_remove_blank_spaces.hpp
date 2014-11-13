/**
 function_remove_blank_spaces.cpp
 Purpose: This files contains the function RemoveBlankSpaces()
 
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 10/22/2014
 */

#ifndef _RemoveBlankSpaces_h
#define _RemoveBlankSpaces_h


// libraries used
#include <string>
////////

// classes used
////////

// functions used
////////

//-----------------------------------------------------
inline std::string RemoveBlankSpaces(const std::string& str){
    /* 
     EXPLANATION:
     Removes all the leading and ending blank spaces
     
     INPUTS:
     str - a string
     
     OUTPUTS:
     a string without the leading and ending blank spaces
     
     GLOBAL VARIABLES USED:
     NONE
     
     CLASSES USED:
     NONE
     
     FUNCITONS USED:
     NONE
     */
    
    // deleting leading blank spaces
    size_t lead_pos = str.find_first_not_of(" ");
    size_t end_pos = str.find_last_not_of(" ");

    return str.substr(lead_pos, end_pos-lead_pos+1);

}

#endif
