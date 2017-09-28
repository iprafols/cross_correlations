/**
 metal_grid.h
 Purpose: This file defines the class MetalGrid. This class contains the variables necessary to store the measurement of the metal contamination grid for a single grid
 
 @author Ignasi Pérez-Ràfols (iprafols@gmail.com)
 @version 1.0 28/09/2017
 
 */

#ifndef _MetalGrid_h
#define _MetalGrid_h

// libraries needed
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
////////

// classes needed
////////

// functions needed
////////

#include "defines.h"
#include "typedefs.h"

class MetalGrid{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs empty object
    MetalGrid(){};
    
    // constructs "bad data" MetalGrid
    MetalGrid(int bad_data);
    
    // constructs object and initializes its variables
    MetalGrid(const size_t& num_bins, const double& metal_wl, const std::string& metal_name, const size_t& flag_verbose_metal_grid);
    
    
    // -------------------------------------------------------------
    // access methods
    
    // access functions for mean_pi_
    std::vector<double> mean_pi() const {return mean_pi_;}
    double mean_pi(size_t index) const;
    
    // access functions for mean_sigma_
    std::vector<double> mean_sigma() const {return mean_sigma_;}
    double mean_sigma(size_t index) const;
    
    // access function for mean_z_
    std::vector<double> mean_z() const {return mean_z_;}
    double mean_z(size_t index) const;
    
    // access function for num_bins_
    size_t num_bins() const {return num_bins_;}
    
    // access functions for weight_
    std::vector<double> weight() const {return weight_;}
    double weight(size_t index) const;
    
    // -------------------------------------------------------------
    // set methods
    
    // set mean_pi_
    void set_mean_pi(size_t index, double value);
    
    // set mean_sigma_
    void set_mean_sigma(size_t index, double value);
    
    // set mean_z
    void set_mean_z(size_t index, double value);
    
    // set weight_
    void set_weight(size_t index, double value);
    
    // -------------------------------------------------------------
    // other methods
    
    // Adds pair contribution to the specified bin
    void AddPair(const size_t& k_index, const double& weight, const double& pi, const double& sigma, const double& z);
    
    // Normalizes the values in the grid
    void Normalize();
    
    // -------------------------------------------------------------
    // operator overload
    
    // adition
    void operator+= (const MetalGrid& other);
    
    // subtraction
    MetalGrid operator- (const MetalGrid& other);
    
    // multiplication
    MetalGrid operator* (const MetalGrid& other);


private:
    // verbose flag
    size_t flag_verbose_metal_grid_;
        
    // mean value of parallel separation in bin
    std::vector<double> mean_pi_;
    
    // mean value of perpendicular separation in bin
    std::vector<double> mean_sigma_;
    
    // mean redshift
    std::vector<double> mean_z_;
    
    // weight
    std::vector<double> weight_;
    
    // number of bins
    size_t num_bins_;
    
    // name of the metal contaminant to compute grids from
    std::string metal_name_;
    
    // name of the metal contaminant to compute grids from
    double metal_wl_;
};


#endif
