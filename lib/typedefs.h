/**
 typedefs.h
 Purpose: This file contains type definition instructions
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 09/23/2014
 
 */

#ifndef _typedefs_h
#define _typedefs_h

#include <map>
#include <utility>
#include <string>
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

typedef std::map<std::pair<int,int>,double> CovMat;
typedef std::map<std::pair<int,int>,double> DistMat;
typedef std::vector<std::pair<std::string, double> > InputMetals;

#endif
