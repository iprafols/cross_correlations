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
#include "global_variables.h"
#include "lya_spectra_dataset.h"
#include "plate_neighbours.h"
#include "plots_object.h"
////////

// functions used
////////



int main(){
    /**
     EXPLANATION:
     Compute the cross correlation of the Lyman-alpha forest and quasars
          
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
          
     CLASSES USED:
     AstroObjectDataset
     GlobalVariables
     PlateNeighbours
     LyaSpectraDataset
     
     FUNCITONS USED:
     NONE
     */
    
    // load time control variables
    time_t start_time,end_time;
    time(&start_time);

    // load global variables and plot object
    std::cout << "Initializing variables" << std::endl;
    const GlobalVariables kGlobalVariables;
    const PlotsObject kPlots(kGlobalVariables.plots());

    // load plate list
    std::cout << "Loading plate list" << std::endl;
    const PlateNeighbours kPlateNeighbours(kGlobalVariables.plate_neighbours());
        
    // load quasar dataset
    std::cout << "Loading quasar dataset" << std::endl;
    AstroObjectDataset object_list(kGlobalVariables);
    std::cout << "Loaded " << object_list.size() << " quasars" << std::endl;
    std::cout << "Plotting quasar dataset information" << std::endl;
    kPlots.PlotRADECDispersion(object_list, true);
    kPlots.PlotZHistogram(object_list, true);
    
    // load spectra dataset
    std::cout << "Loading spectra dataset" << std::endl;
    LyaSpectraDataset spectra_list(kGlobalVariables);
    std::cout << "Loaded " << spectra_list.size() << " spectra" << std::endl;
    std::cout << "Plotting spectra dataset information" << std::endl;
    kPlots.PlotRADECDispersion(spectra_list, true);
    kPlots.PlotZHistogram(spectra_list, true);
    
    // compute distances to the objects (in Mpc/h)
    {
        std::cout << "Computing distances (in Mpc/h) to objects" << std::endl;
        InterpolationMap redshift_distance_map(kGlobalVariables);
        object_list.SetDistances(redshift_distance_map);
        spectra_list.SetDistances(redshift_distance_map);
    }
    
    // compute the cross-correlation
    std::cout << "Creating object in where cross-correlation results will be stored" << std::endl;
        CorrelationResults results(kGlobalVariables, kPlateNeighbours);
    std::cout << "Computing the cross-correlation" << std::endl;
    results.ComputeCrossCorrelation(object_list, spectra_list, kGlobalVariables);
    std::cout << "Plotting cross-correlation" << std::endl;
    kPlots.PlotCrossCorrelation(results, true);
        
    
    // save the cross-correlation results
    // --> cCorrelationResults::Save()
    
    // remove datasets from memory
    
    
    // compute covariance matrix
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
    
    // create the python scripts for plotting
    
    
    // display time required to run the program
    std::cout << "End of program" << std::endl;
    time(&end_time);
    double time_spent = difftime(end_time, start_time);
    std::cout << "The program lasted " << time_spent << " seconds. This corresponds to " << time_spent/60.0 << " minutes or " << time_spent/3600.0 << " hours" << std::endl;
    
    return 0;
}