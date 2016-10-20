/**
 main_correlation.cpp
 Purpose: Compute the variance, the contribution to the covariance matrix as a function of redshift and pixel separation
 
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
     This function computes the 1D correlation of the Lyman-alpha mean transmission as a function of redshift and pixel separation
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
        spectra_list.reset(new LyaSpectraDataset(input));
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
    
    // define vectors to store the correlation as a function of redshift and pixel separation
    size_t num_points_interpolation = input.num_points_interpolation();
    size_t max_pixels_separation = input.pixels_separation() + 1;
    std::vector<double> z_correlation(num_points_interpolation, 0.0);
    double z_step = (input.z_max_interpolation() - input.z_min_interpolation())/float(z_correlation.size());
    for (size_t index = 0; index < z_correlation.size(); index++){
        z_correlation[index] = input.z_min_interpolation() + z_step*float(index);
    }
    std::vector<std::vector<double> > correlation;
    {
        std::vector<double> correlation_separation_N(num_points_interpolation, 0.0);
        correlation.resize(max_pixels_separation, correlation_separation_N);
    }
    std::vector<std::vector<double> > total_weight = correlation;
    
    // define copies of correlation and total_weight for each of the threads
    int num_threads = atoi(std::getenv("OMP_NUM_THREADS"));
    std::vector<std::vector<std::vector<double> > > correlation_thread;
    correlation_thread.resize(num_threads, correlation);
    std::vector<std::vector<std::vector<double> > > total_weight_thread;
    total_weight_thread.resize(num_threads, total_weight);
    
    // load list of plates
    PlatesMapVector<LyaSpectrum>::map plates_list = (*spectra_list).list();
    
    // compute the correlation for different pixel separations
    std::cout << "Computing the correlation for different pixel separations" << std::endl;
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
                
                // loop over pixel separation
                for (size_t pixel_separation = 0; pixel_separation < max_pixels_separation; pixel_separation ++){
                    // check that we are not reaching the end of the spectrum
                    if (pixel_separation + pixel_number < spectrum.size()){
                        // locate the bin in z the pixel belongs to
                        double mean_z = (spectrum[pixel_number].z() + spectrum[pixel_number + pixel_separation].z())/2.0;
                        int z_index = int((mean_z - input.z_min_interpolation())/z_step);
                        if ((z_index < 0) or (z_index > z_correlation.size())){
                            #pragma omp critical (cerr)
                            {
                                std::cerr << "Bad index, z_index = " << z_index << std::endl;
                            }
                        }
                        else{
                            weight = spectrum[pixel_number].weight()*spectrum[pixel_number + pixel_separation].weight();
                            correlation_thread[thread_num][pixel_separation][z_index] += spectrum[pixel_number].delta()*spectrum[pixel_number + pixel_separation].delta()*weight;
                            total_weight_thread[thread_num][pixel_separation][z_index] += weight;
                        }
                        
                        // TODO: remove this test
                        /*if (mean_z >= 3.426){
                            #pragma omp critical (cout)
                            {
                                std::cout << "In thread " << thread_num << ": mean_z = " << mean_z;
                                std::cout << "delta1 = " << spectrum[pixel_number].delta();
                                std::cout << ", delta2 = " << spectrum[pixel_number + pixel_separation].delta();
                                std::cout << ", weight = " << weight;
                                std::cout << ", corr[" << pixel_separation << "] = " << correlation_thread[thread_num][pixel_separation][z_index];
                                std::cout << ", weight[" << pixel_separation << "] = " << total_weight_thread[thread_num][pixel_separation][z_index];
                                std::cout << std::endl;
                            }
                        }*/
                        // end of test
                    }
                }
            }
        }
    }
    
    // combine the measurements from the different threads
    for (size_t pixel_separation =  0; pixel_separation < max_pixels_separation; pixel_separation ++){
        for (size_t index = 0; index < num_points_interpolation; index ++){
            for (size_t thread_num = 0; thread_num < correlation_thread.size(); thread_num ++){
                correlation[pixel_separation][index] += correlation_thread[thread_num][pixel_separation][index];
                total_weight[pixel_separation][index] += total_weight_thread[thread_num][pixel_separation][index];
            }
        }
    }
    
    // normalize the correlation
    for (size_t pixel_separation =  0; pixel_separation < max_pixels_separation; pixel_separation ++){
        for (int index =  0; index < num_points_interpolation; index ++){
            // check that the weight is not zero
            if (total_weight[pixel_separation][index] == 0.0){
                std::cerr << "Warning : In compute_lya_1d : For pixel separation " << pixel_separation << " index " << index;
                std::cerr << " shows zero weight. Consider reducing the num_points_interpolation or pixels_separation." << std::endl;
            }
            else{
                correlation[pixel_separation][index] /= total_weight[pixel_separation][index];
            }
        }
    }
    
    // save results
    std::string filename = input.lya_auto_correlation();
    std::cout << "Saving results to " << filename << std::endl;
    std::ofstream file;
    file.open(filename.c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
    if (file.is_open()){
        // print header
        file << "z";
        for (size_t pixel_separation = 0; pixel_separation < max_pixels_separation; pixel_separation ++){
            file << " xi(" << pixel_separation << ")";
        }
        file << std::endl;
        // print results
        for (size_t index = 0; index < z_correlation.size(); index ++){
            file << z_correlation[index];
            for (size_t pixel_separation = 0; pixel_separation < max_pixels_separation; pixel_separation ++){
                file << " " << correlation[pixel_separation][index];
            }
            file <<  std::endl;
        }
        
        file.close();
    }
    else{
        std::cerr << "Error : In compute_lya_1d : Unable to open file:" << std::endl << filename << std::endl;
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