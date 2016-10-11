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
#include "input.h"
#include "dataset.h"
#include "lya_auto_interpolation_map.h"
////////

// functions needed
////////

class PlotsObject{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs object and initializes its variables
    PlotsObject(const Input& input);
    
    // -------------------------------------------------------------
    // access methods
    
    // access functions for plots_dir_
    std::string plots_dir() const {return plots_dir_;}
    
    // -------------------------------------------------------------
    // other methods
    
    // Plots the correlation function
    void PlotCrossCorrelation(const CorrelationResults& res, const Input& input, const bool update_script = false) const;
    
    // Plots the lya autocorrealtion
    void PlotLyaAuto(const std::vector<LyaAutoInterpolationMap>& lya_auto, const bool update_script = false) const;
 
        
    // Plots the RA-DEC dispersion for the given objects
    void PlotRADECDispersion(const Dataset& dataset, const bool update_script = false) const;
        
    // Plots the RA-DEC dispersion for the given objects
    void PlotZHistogram(const Dataset& dataset, const bool update_script = false) const;
    
    
    
    
    
    
private:
    // path where the plots will be stored
    std::string plots_dir_;
    
    // output base name
    std::string output_base_name_;
    
    
    // -------------------------------------------------------------
    // methods
    
    // Writes a makefile for plotting
    void MakePlottingMakefile() const;
};


#endif
