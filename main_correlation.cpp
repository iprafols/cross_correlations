/**
 main_correlation.cpp
 Purpose: Compute the cross correlation of the Lyman-alpha forest and quasars. Future versions should compute the cross-correlation of two species in general
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 06/17/2014
 */

// libraries used
#include <iostream>
#include <time.h>
////////

// classes used
#include "astro_object.h"
#include "astro_object_dataset.h"
#include "correlation_plate.h"
#include "correlation_results.h"
#include "covariance_matrix.h"
#include "input.h"
#include "lya_spectra_dataset.h"
#include "plate_neighbours.h"
#include "plots_object.h"
////////

// functions used
#include "defines.h"
////////

#include "typedefs.h"


int main(int argc, char *argv[]){
    /**
     EXPLANATION:
     Compute the cross correlation of the Lyman-alpha forest and quasars
          
     INPUTS:
     input_file[optional] - a file containing the input settings
     
     OUTPUTS:
     NONE
          
     CLASSES USED:
     AstroObjectDataset
     Input
     PlateNeighbours
     LyaSpectraDataset
     
     FUNCITONS USED:
     ComputePlateNeighbours
     */
    
    // load time control variables
    time_t start_time,end_time;
    time(&start_time);

    // load global variables and plot object
    std::cout << "Initializing variables" << std::endl;
    string input_filename = "";
    if (argc > 1){
        input_filename += argv[1];
    }
    Input input(input_filename);
    const PlotsObject kPlots(input);

    // check whether or not the plate list needs to be computed
    if (input.flag_compute_plate_neighbours()){
        ComputePlateNeighbours(input);
    }
    // load plate list
    std::cout << "Loading plate list" << std::endl;
    const PlateNeighbours kPlateNeighbours(input);
        
    // initialize cross-correlation results
    CorrelationResults results(input, kPlateNeighbours);
    
    {
        // load quasar dataset
        std::cout << "Loading quasar dataset" << std::endl;
        AstroObjectDataset object_list(input);
        std::cout << "Loaded " << object_list.size() << " quasars" << std::endl;
        std::cout << "Plotting quasar dataset information" << std::endl;
        kPlots.PlotRADECDispersion(object_list, true);
        kPlots.PlotZHistogram(object_list, true);
        
        // load spectra dataset
        std::cout << "Loading spectra dataset" << std::endl;
        LyaSpectraDataset spectra_list(input);
        std::cout << "Loaded " << spectra_list.size() << " spectra" << std::endl;
        std::cout << "Plotting spectra dataset information" << std::endl;
        kPlots.PlotRADECDispersion(spectra_list, true);
        kPlots.PlotZHistogram(spectra_list, true);
        
        // compute distances to the objects (in Mpc/h)
        {
            std::cout << "Computing distances (in Mpc/h) to objects" << std::endl;
            InterpolationMap redshift_distance_map(input);
            object_list.SetDistances(redshift_distance_map);
            spectra_list.SetDistances(redshift_distance_map);
        }
        
        // compute the cross-correlation
        std::cout << "Computing the cross-correlation" << std::endl;
        //CorrelationResults results(input, kPlateNeighbours);
        results.ComputeCrossCorrelation(object_list, spectra_list, input);
        std::cout << "Plotting cross-correlation" << std::endl;
        kPlots.PlotCrossCorrelation(results, true);
    }
    /*// load quasar dataset
    std::cout << "Loading quasar dataset" << std::endl;
    AstroObjectDataset object_list(input);
    std::cout << "Loaded " << object_list.size() << " quasars" << std::endl;
    std::cout << "Plotting quasar dataset information" << std::endl;
    kPlots.PlotRADECDispersion(object_list, true);
    kPlots.PlotZHistogram(object_list, true);
    
    // load spectra dataset
    std::cout << "Loading spectra dataset" << std::endl;
    LyaSpectraDataset spectra_list(input);
    std::cout << "Loaded " << spectra_list.size() << " spectra" << std::endl;
    std::cout << "Plotting spectra dataset information" << std::endl;
    kPlots.PlotRADECDispersion(spectra_list, true);
    kPlots.PlotZHistogram(spectra_list, true);
    
    // compute distances to the objects (in Mpc/h)
    {
        std::cout << "Computing distances (in Mpc/h) to objects" << std::endl;
        InterpolationMap redshift_distance_map(input);
        object_list.SetDistances(redshift_distance_map);
        spectra_list.SetDistances(redshift_distance_map);
    }
    
    // compute the cross-correlation
    std::cout << "Computing the cross-correlation" << std::endl;
    //CorrelationResults results(input, kPlateNeighbours);
    results.ComputeCrossCorrelation(object_list, spectra_list, input);
    std::cout << "Plotting cross-correlation" << std::endl;
    kPlots.PlotCrossCorrelation(results, true);
        
    // remove datasets from memory
    */
    
    // compute covariance matrix
    std::cout << "Computing covariance matrix" << std::endl;
    CovarianceMatrix cov_mat(input);
    if (input.flag_compute_bootstrap()){
        cov_mat.ComputeBootstrapCovMat(results.bootstrap());
    }
    
    {
        /*
         load objects where the results are stored
        --> cCovMatrixResults
         
         for all cross-correlation bins:
            
            load pixel list 1
            --> cExtendedForestPixelDataset
                cExtendedForestPixelDataset <-- Dataset
         
                --> vector<cExtendedForestPixel> 
                    cExtendedForestPixel <-- ForestPixel + ra, dec
         
            for all cross-correlation bins AFTER this one:
            
                * load pixel list 2
                --> cExtendedForestPixelDataset
         
                * do N times:
                    + randomly select item from pixel_list1 
                    + retrieve plate neighbours of the selected item's plate
                    --> cPlateNeighbours::GetNeighboursList
                    + randomly select item from the subset of pixel_list2 comprised of all the neighbouring plates
                    + add to covariance matrix
            

         */
    }
    
    
    // save the covariance matrix results
    //--> cCovMatrixResults::Save()
    
    // run the plotting scripts
    
    
    // display time required to run the program
    std::cout << "End of program" << std::endl;
    time(&end_time);
    double time_spent = difftime(end_time, start_time);
    if (time_spent < 60.0){
        std::cout << "Program lasted " << time_spent << " seconds" << std::endl;
    }
    else if (time_spent < 3600.0){
        std::cout << "Program lasted " << time_spent/60.0 << " minutes" << std::endl;
    }
    else{
        std::cout << "Program lasted " << time_spent/3600.0 << " hours" << std::endl;
    }
    
    return 0;
}