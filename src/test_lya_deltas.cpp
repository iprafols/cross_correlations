/**
 main_correlation.cpp
 Purpose: Compute the cross correlation of the Lyman-alpha forest and quasars. Future versions should compute the cross-correlation of two species in general
 
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
     This function computes the variance of the Lyman-alpha mean transmission as a function of wavelength and compares it with its estimation
     from the weights given in Niolas Busca files. Note that this program should work with typical ini file developed to pass to 
     correlation.run but that not all options are used here
          
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
    
    // define vectors to store the variance as a function of redshift
    size_t num_bins_z = 500;
    std::vector<double> z_list(num_bins_z, 0.0);
    double z_step = (input.z_max() - input.z_min())/float(z_list.size());
    for (size_t index = 0; index < z_list.size(); index++){
        z_list[index] = input.z_min() + z_step*float(index);
    }
    std::vector<double> variance(num_bins_z, 0.0);
    std::vector<double> correlation_first_neighbours(num_bins_z, 0.0);
    std::vector<double> total_weight(num_bins_z, 0.0);
    
    // define varaibles to compute the variance estimate
    std::vector<double> vairance_estimate(num_bins_z, 0.0);
    double gamma_halfs = 3.8/2.0;
    double one_plus_z0 = 1.0 + 2.25;
    double one_plus_z0_to_the_half_gamma = pow(1.0 + 2.25, gamma_halfs);
    
    // define copies of mean, variance, vairance_estimate and total_weight for each of the threads
    int num_threads = atoi(std::getenv("OMP_NUM_THREADS"));
    std::vector<std::vector<double> > variance_thread;
    variance_thread.resize(num_threads, variance);
    std::vector<std::vector<double> > correlation_first_neighbours_thread;
    correlation_first_neighbours_thread.resize(num_threads, correlation_first_neighbours);
    std::vector<std::vector<double> > vairance_estimate_thread;
    vairance_estimate_thread.resize(num_threads, vairance_estimate);
    std::vector<std::vector<double> > total_weight_thread;
    total_weight_thread.resize(num_threads, total_weight);
    
    // load list of plates
    PlatesMapVector<LyaSpectrum>::map plates_list = (*spectra_list).list();
    
    
    // compute the variance
    std::cout << "Computing variance of delta" << std::endl;
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
                int z_index = int((spectrum[pixel_number].z() - input.z_min())/z_step);
                if ((z_index < 0) or (z_index > z_list.size())){
                    #pragma omp critical (cerr)
                    {
                        std::cerr << "Bad index, z_index = " << z_index << std::endl;
                    }
                }
                // add to average
                else{
                    weight = spectrum[pixel_number].weight()*spectrum[pixel_number].weight();
                    variance_thread[thread_num][z_index] += spectrum[pixel_number].delta()*spectrum[pixel_number].delta()*weight;
                    vairance_estimate_thread[thread_num][z_index] += pow(1.0+spectrum[pixel_number].z(), gamma_halfs)*spectrum[pixel_number].weight()/one_plus_z0_to_the_half_gamma;
                    total_weight_thread[thread_num][z_index] += weight;
                }
            }
        }
    }
    
    // combine the measurements from the different threads and reset total_weight_thread
    for (int index =  0; index < variance.size(); index ++){
        for (size_t thread_num = 0; thread_num < variance_thread.size(); thread_num ++){
            variance[index] += variance_thread[thread_num][index];
            vairance_estimate[index] += vairance_estimate_thread[thread_num][index];
            total_weight[index] += total_weight_thread[thread_num][index];
            total_weight_thread[thread_num][index] = 0.0;
        }
    }
    
    // normalize the variance and reset total_weight
    for (int index =  0; index < variance.size(); index ++){
        variance[index] /= total_weight[index];
        vairance_estimate[index] = vairance_estimate[index]/total_weight[index];
        total_weight[index] = 0.0;
    }
    
    // finally, compute the correlation with the first neighbour
    std::cout << "Computing correlation of delta with the first neighbour" << std::endl;
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
            for (size_t pixel_number = 0; pixel_number < spectrum.size() - 1; pixel_number++){
                
                // locate the bin in z the pixel belongs to
                int z_index = int((spectrum[pixel_number].z() - input.z_min())/z_step);
                int z_index_first_neighbour = int((spectrum[pixel_number + 1].z() - input.z_min())/z_step);
                if ((z_index < 0) or (z_index > z_list.size()) or (z_index_first_neighbour < 0) or (z_index_first_neighbour > z_list.size())){
                    #pragma omp critical (cerr)
                    {
                        std::cerr << "Bad index, z_index = " << z_index << std::endl;
                    }
                }
                // add to average
                else{
                    weight = spectrum[pixel_number].weight()*spectrum[pixel_number + 1].weight();
                    correlation_first_neighbours_thread[thread_num][z_index] += spectrum[pixel_number].delta()*spectrum[pixel_number + 1].delta()*weight;
                    total_weight_thread[thread_num][z_index] += weight;
                    
                }
            }
        }
    }
    
    // combine the measurements from the different threads
    for (int index =  0; index < correlation_first_neighbours.size(); index ++){
        for (size_t thread_num = 0; thread_num < correlation_first_neighbours_thread.size(); thread_num ++){
            correlation_first_neighbours[index] += correlation_first_neighbours_thread[thread_num][index];
            total_weight[index] += total_weight_thread[thread_num][index];
        }
    }
    
    // normalize the variance
    for (int index =  0; index < correlation_first_neighbours.size(); index ++){
        correlation_first_neighbours[index] /= total_weight[index];
    }
    
    // now compute the xi^1D estimate from Palanque-Delabrouille
    LyaAutoInterpolationMap lya_auto = LyaAutoInterpolationMap(input, 1);

    
    // save results
    std::cout << "Saving results" << std::endl;
    std::string filename("lya_deltas_variance.dat");
    std::ofstream file;
    file.open(filename.c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
    if (file.is_open()){
        file << "z variance variance_estimate corr_first_neighbours xi1DNatalie" << std::endl;
        for (int index = 0; index < z_list.size(); index ++){
            file << z_list[index] << " " << variance[index] << " " << vairance_estimate[index] << " " << correlation_first_neighbours[index] << " " << lya_auto.LinearInterpolation(z_list[index]) <<  std::endl;
        }
        
        file.close();
    }
    else{
        std::cout << "Error : In test_lya_deltas : Unable to open file:" << std::endl << filename << std::endl;
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