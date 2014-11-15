/**
 typedefs.h
 Purpose: This file contains type definition instructions
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 09/23/2014
 
 */

#ifndef _typedefs_h
#define _typedefs_h

#include <map>
#include <vector>

template <typename T> 
struct PlatesMapVector{
    typedef std::map<int, std::vector<T> > map;
};


template <typename T>
struct PlatesMapSimple{
    typedef std::map<int, T> map;
};

typedef std::map<std::string, bool> InputFlag;



#endif
