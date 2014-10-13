/**
 plots_object.h
 Purpose: This file defines the class PlotsObject. This class contains the functions required to build the ploting python scripts
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 09/25/2014
 
 */

#ifndef _Plots_h
#define _Plots_h

// libraries needed
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
////////

// classes needed
#include "correlation_results.h"
#include "dataset.h"
////////

// functions needed
////////

class PlotsObject{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs object and initializes its variables
    PlotsObject(std::string plots_dir);
    
    // -------------------------------------------------------------
    // access methods
    
    // access functions for plots_dir_
    std::string plots_dir() const {return plots_dir_;}
    
    // -------------------------------------------------------------
    // other methods
    
    // Plots the correlation function
    void PlotCrossCorrelation(const CorrelationResults& res, const bool update_script = false) const;
    
    // Plots the RA-DEC dispersion for the given objects
    void PlotRADECDispersion(Dataset& dataset, const bool update_script = false) const;
        
    // Plots the RA-DEC dispersion for the given objects
    void PlotZHistogram(Dataset& dataset, const bool update_script = false) const;
    
    
    
    
    
    
private:
    // path where the plots will be stored
    std::string plots_dir_;
    
};


#endif
