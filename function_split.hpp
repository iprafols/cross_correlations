/**
 function_split.cpp
 Purpose: This files contains the function Split()

 
 @author Ignasi Pérez-Ràfols
 @version 1.0 06/17/2014
 */

#ifndef _Split_h
#define _Split_h


// libraries used
#include <iostream>
#include <vector>
#include <string>
////////

// classes used
////////

// functions used
////////

//-----------------------------------------------------
inline std::vector<std::string> Split(const std::string& line,const std::string& delimiter){
    /* 
     EXPLANATION:
     splits a string into a vector of char*, which are sequences of contiguous characters separated by delimiter
     INPUTS:
     line - string to be splited
     delimiter - string denoting the separation criteria
     OUTPUTS:
     cols - char* vector
     GLOBAL VARIABLES USED:
     NONE
     CLASSES USED:
     NONE
     FUNCITONS USED:
     NONE
     */
    
    // declaring variables
    std::vector<std::string> cols;
    char* col;
    
    // splitting string
    col = strtok((char*)line.c_str(),delimiter.c_str());
    while (col){
        cols.push_back(col);
        col = strtok(NULL,delimiter.c_str());
    }
    
    return cols;
}

#endif
