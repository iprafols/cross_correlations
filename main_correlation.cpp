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
#include "baofit_setup.h"
#include "correlation_plate.h"
#include "correlation_results.h"
#include "covariance_matrix.h"
#include "dla_dataset.h"
#include "input.h"
#include "lya_spectra_dataset.h"
#include "plate_neighbours.h"
#include "plots_object.h"
#include "quasar_dataset.h"
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
    size_t flag_verbose_main = input.flag_verbose_main();
    const PlotsObject kPlots(input);

    // check whether or not the plate list needs to be computed
    if (input.flag_compute_plate_neighbours()){
        ComputePlateNeighbours(input);
    }
    // load plate list
    const PlateNeighbours kPlateNeighbours(input);
        
    // initialize cross-correlation results
    CorrelationResults results(input, kPlateNeighbours);
    
    if (input.flag_compute_cross_correlation() or input.flag_load_only() or input.flag_plot_catalog_info()){
        // load quasar dataset
        std::auto_ptr<AstroObjectDataset> object_list;
        if (input.dataset1_type() == "quasar"){
            object_list.reset(new QuasarDataset(input));
        }
        else if (input.dataset1_type() == "dla"){
            object_list.reset(new DLADataset(input));
        }
        else{
            std::cout << "The selected type for dataset1 is not enabled. Current options are: " << input.dataset1_type_options() << std::endl;
            return 0;
        }
        if (input.flag_plot_catalog_info()){
            if (flag_verbose_main >= 2){
                std::cout << "Plotting dataset1 information" << std::endl;
            }
            kPlots.PlotRADECDispersion(*object_list, true);
            kPlots.PlotZHistogram(*object_list, true);
        }
        
        // load spectra dataset
        LyaSpectraDataset spectra_list(input);
        if (input.flag_plot_catalog_info()){
            if (flag_verbose_main >= 2){
                std::cout << "Plotting spectra dataset information" << std::endl;
            }
            kPlots.PlotRADECDispersion(spectra_list, true);
            kPlots.PlotZHistogram(spectra_list, true);
        }
        
        // check load_only flag and end program if set
        if (input.flag_load_only()){
            if (flag_verbose_main >= 1){
                std::cout << "End of program" << std::endl;
            }
            time(&end_time);
            double time_spent = difftime(end_time, start_time);
            if (flag_verbose_main >= 1){
                if (time_spent < 60.0){
                    std::cout << "Program lasted " << time_spent << " seconds" << std::endl;
                }
                else if (time_spent < 3600.0){
                    std::cout << "Program lasted " << time_spent/60.0 << " minutes" << std::endl;
                }
                else{
                    std::cout << "Program lasted " << time_spent/3600.0 << " hours" << std::endl;
                }
            }
            return 1;
        }
        
        // compute the cross-correlation
        if (input.flag_compute_cross_correlation()){
            {
                // compute distances to the objects (in Mpc/h)
                if (flag_verbose_main >= 1){
                    std::cout << "Computing distances (in Mpc/h) to objects" << std::endl;
                }
                InterpolationMap redshift_distance_map(input);
                (*object_list).SetDistances(redshift_distance_map);
                spectra_list.SetDistances(redshift_distance_map);
            }
            // compute cross-correlations
            results.ComputeCrossCorrelation(*object_list, spectra_list, input);
            if (flag_verbose_main >= 2){
                std::cout << "Plotting cross-correlation" << std::endl;
            }
            kPlots.PlotCrossCorrelation(results, input, true);
        }
    }
    
    if (input.flag_compute_covariance()){
        // compute covariance matrix
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


    }
    
    // setup BAOFIT configuration and run fitting program
    BaofitSetup baofit(input);
    if (input.flag_set_baofit()){
        if (input.flag_compute_bootstrap()){
            baofit.Set(input,true);
        }
        // baofit.Set(input); // full covariance matrix needs to be computed first
    }
    if (input.flag_run_baofit()){
        if (input.flag_compute_bootstrap()){
            baofit.Run(input, true);
        }
        // baofit.Run(input); // full covariance matrix needs to be computed first
    }
    if (input.flag_set_baofit_best_fit()){
        if (input.flag_compute_bootstrap()){
            baofit.SetBestFit(input,true);
        }
        // baofit.SetBestFit(input); // full covariance matrix needs to be computed first
    }
    if (input.flag_run_baofit_best_fit()){
        if (input.flag_compute_bootstrap()){
            baofit.RunBestFit(input, true);
        }
        // baofit.RunBestFit(input); // full covariance matrix needs to be computed first        
    }
    
    // make the plots
    
    
    // display time required to run the program
    if (flag_verbose_main >= 1){
        std::cout << "End of program" << std::endl;
    }
    time(&end_time);
    double time_spent = difftime(end_time, start_time);
    if (flag_verbose_main >= 1){
        if (time_spent < 60.0){
            std::cout << "Program lasted " << time_spent << " seconds" << std::endl;
        }
        else if (time_spent < 3600.0){
            std::cout << "Program lasted " << time_spent/60.0 << " minutes" << std::endl;
        }
        else{
            std::cout << "Program lasted " << time_spent/3600.0 << " hours" << std::endl;
        }
    }
    return 1;
}