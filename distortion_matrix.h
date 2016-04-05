/**
 distortion_matrix.h
 Purpose: This file defines the class DistortionMatrix. This class contains the variables necessary to compute and store the distortion matrix of the measurement.
 
 @author Ignasi Pérez-Ràfols (iprafols@icc.ub.edu)
 @version 1.0 04/04/2016
 
 */

#ifndef _DistortionMatrix_h
#define _DistortionMatrix_h

// libraries needed
#include <string>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>
#include "omp.h"

#include <iomanip>
////////

// classes needed
#include "distortion_plate.h"
#include "input.h"
#include "plate_neighbours.h"
#include "pair_dataset.h"
#include "spectra_dataset.h"
////////

// functions needed
////////

#include "typedefs.h"
#include "defines.h"

class DistortionMatrix{
    
public:
    // -------------------------------------------------------------
    // constructors
    
    // constructs object and initializes its variables
    DistortionMatrix(const Input& input, const PlateNeighbours& kPlateNeighbours);
    
    
    // -------------------------------------------------------------
    // access methods
    
    // access function for dist_mat_
    CovMat dist_mat() const {return dist_mat_;}
    double dist_mat(size_t i, size_t j) const;
    
    // access function for flag_verbose_distortion_matrix_
    size_t flag_verbose_distortion_matrix() const {return flag_verbose_distortion_matrix_;}
    
    // access function for num_bins_
    size_t num_bins() const {return num_bins_;}

    // access function for output_base_name_
    std::string output_base_name() const {return output_base_name_;}

    
    
    // -------------------------------------------------------------
    // other methods
    
    // compute distortion matrix
    void ComputeDistMat(const AstroObjectDataset& object_list, const SpectraDataset& spectra_list, const Input& input, const PlateNeighbours& kPlateNeighbours);

    
    
    
private:
    // distortion matrix
    CovMat dist_mat_;
    
    // vector containing the distortion matrix computed with different threads
    std::vector<DistortionPlate> distortion_threads_;
    
    // distortion_matrix verbose flag
    size_t flag_verbose_distortion_matrix_;
    
    // normalized distortion matrix
    DistortionPlate normalized_dist_mat_;
    
    // number of bins
    size_t num_bins_;
    
    // output base name
    std::string output_base_name_;
    
    // list of plates
    std::vector<int> plates_list_;
    
    // number of plates that have to be skipped
    int skip_plates_;
    
    // -------------------------------------------------------------
    // methods
    
    // normalize the distortion matrix
    void NormalizeDistMat();
    
    // save the distortion matrix
    void SaveDistMat();
    
};



#endif