/**
 main_correlation.cpp
 Purpose: Compute the average of projected deltas as a function of redshift
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 06/17/2014
 */

// libraries used
#include <iostream>
#include <time.h>
#include "omp.h"
////////

// classes used
#include "astro_object.h"
#include "astro_object_dataset.h"
#include "civ_spectra_dataset.h"
#include "correlation_plate.h"
#include "correlation_results.h"
#include "covariance_matrix.h"
#include "dla_dataset.h"
#include "distortion_matrix.h"
#include "input.h"
#include "lya_spectra_dataset.h"
#include "lya_auto_interpolation_map.h"
#include "plate_neighbours.h"
#include "plots_object.h"
#include "quasar_dataset.h"
#include "spectra_dataset.h"
#include "z_dist_interpolation_map.h"
////////

// functions used
#include "defines.h"
////////

#include "typedefs.h"



int main(int argc, char *argv[]){
    /**
     EXPLANATION:
     This function computes the average of the (projected) deltas as a function of redshift
     Note that this program should work with typical ini file developed to pass to correlation.run but that not all options are used here.
     This program should be run for each set of the paramenters z_min_interpolation, z_max_interpolation, num_points_interpolation, 
     or pixel_separation
          
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
    
    // load spectra dataset
    std::auto_ptr<SpectraDataset> spectra_list;
    if (input.dataset2_type() == "lya"){
        spectra_list.reset(new LyaSpectraDataset(input, true));
    }
    else if (input.dataset2_type() == "civ"){
        spectra_list.reset(new CIVSpectraDataset(input));
    }
    else{
        std::cout << "Error : The selected type for dataset2 is not enabled. Current options are: " << input.dataset2_type_options() << std::endl;
        return 1;
    }
    if (input.flag_plot_catalog_info()){
        if (flag_verbose_main >= 2){
            std::cout << "Plotting spectra dataset information" << std::endl;
        }
        kPlots.PlotRADECDispersion(*spectra_list, true);
        kPlots.PlotZHistogram(*spectra_list, true);
        if (flag_verbose_main >= 2){
            std::cout << "done" << std::endl;
        }
    }
    
    // define vectors to store the average projected delta as a function of redshift
    size_t num_points_interpolation = input.num_points_interpolation();
    std::vector<double> z_mean_delta(num_points_interpolation, 0.0);
    double z_step = (input.z_max_interpolation() - input.z_min_interpolation())/float(z_mean_delta.size());
    for (size_t index = 0; index < z_mean_delta.size(); index++){
        z_mean_delta[index] = input.z_min_interpolation() + z_step*float(index);
    }
    std::vector<double> mean_delta(num_points_interpolation, 0.0);
    std::vector<double> total_weight(num_points_interpolation, 0.0);
    
    // define copies of mean_delta and total_weight for each of the threads
    int num_threads = atoi(std::getenv("OMP_NUM_THREADS"));
    std::vector<std::vector<double> > mean_delta_thread;
    mean_delta_thread.resize(num_threads, mean_delta);
    std::vector<std::vector<double> > total_weight_thread;
    total_weight_thread.resize(num_threads, total_weight);
    
    // load list of plates
    PlatesMapVector<LyaSpectrum>::map plates_list = (*spectra_list).list();
    
    // compute the average of projected deltas
    std::cout << "Computing the average of projected deltas" << std::endl;
    // loop over plates
    for (PlatesMapVector<LyaSpectrum>::map::iterator it = plates_list.begin(); it != plates_list.end(); it ++){
        
        std::vector<LyaSpectrum> spec_list = (*it).second;
        // loop over spectra
        #pragma omp parallel for schedule(dynamic)
        for (size_t spectrum_number = 0; spectrum_number < spec_list.size(); spectrum_number ++){
            // get the thread number
            int thread_num = omp_get_thread_num();
            double weight;
            
            std::vector<LyaPixel> spectrum = spec_list[spectrum_number].spectrum();
            // loop over pixel
            for (size_t pixel_number = 0; pixel_number < spectrum.size(); pixel_number++){
                // locate the bin in z the pixel belongs to
                int z_index = int((spectrum[pixel_number].z() - input.z_min_interpolation())/z_step);
                if ((z_index < 0) or (z_index > z_mean_delta.size())){
                    #pragma omp critical (cerr)
                    {
                        std::cerr << "Bad index, z_index = " << z_index << std::endl;
                    }
                }
                else{
                    weight = spectrum[pixel_number].weight();
                    mean_delta_thread[thread_num][z_index] += spectrum[pixel_number].delta()*weight;
                    total_weight_thread[thread_num][z_index] += weight;
                }
            }
        }
    }
    
    // combine the measurements from the different threads
    for (size_t index = 0; index < num_points_interpolation; index ++){
        for (size_t thread_num = 0; thread_num < mean_delta_thread.size(); thread_num ++){
            mean_delta[index] += mean_delta_thread[thread_num][index];
            total_weight[index] += total_weight_thread[thread_num][index];
        }
    }
    
    // normalize the correlation
    for (int index =  0; index < num_points_interpolation; index ++){
        // check that the weight is not zero
        if (total_weight[index] == 0.0){
            std::cerr << "Warning : In compute_projection_correction : Index " << index;
            std::cerr << " shows zero weight. Consider reducing num_points_interpolation." << std::endl;
        }
        else{
            mean_delta[index] /= total_weight[index];
        }
    }
    
    // save results
    std::string filename = input.lya_projection_correction();
    std::cout << "Saving results to " << filename << std::endl;
    std::ofstream file;
    file.open(filename.c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
    if (file.is_open()){
        // print header
        file << "z mean_delta" << std::endl;
        // print results
        for (size_t index = 0; index < z_mean_delta.size(); index ++){
            file << z_mean_delta[index] << " " << mean_delta[index] <<  std::endl;
        }
        
        file.close();
    }
    else{
        std::cerr << "Error : In compute_projection_correction : Unable to open file:" << std::endl << filename << std::endl;
    }
    
    
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
    return 0;
}