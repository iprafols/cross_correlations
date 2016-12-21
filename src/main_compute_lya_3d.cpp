/**
 main_correlation.cpp
 Purpose: Compute the contribution to the covariance matrix as a function of parallel and perpendicular distance for
 pairs-of-pairs in different lines-of-sight
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 06/17/2014
 */

// libraries used
#include <iostream>
#include <time.h>
#include <vector>
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
     This function computes the 3D correlation of the Lyman-alpha mean transmission as a function of parallel and perpendicular separation
     Note that this program should work with typical ini file developed to pass to correlation.run but that not all options are used here.
     This program should be run for each set of the paramenters max_pi_auto, step_pi_auto, num_pi_bins_auto, max_sigma_auto, step_sigma_auto,
     and num_sigma_bins_auto
          
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
    // load list of plates
    const PlateNeighbours kPlateNeighbours(input);
    std::vector<int> plates_list = kPlateNeighbours.GetPlatesList();
    
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
    
    // load binning settings
    double max_pi = input.max_pi_auto();
    double max_sigma = input.max_sigma_auto();
    double step_pi = input.step_pi_auto();
    double step_sigma = input.step_sigma_auto();
    double num_pi_bins = input.num_pi_bins_auto();
    double num_sigma_bins = input.num_sigma_bins_auto();
    double num_bins = input.num_bins_auto();
    
    // define vectors to store the correlation as a function of parallel and perpendicular separation
    std::vector<double> correlation, total_weight;
    correlation.resize(num_bins, 0.0);
    total_weight.resize(num_bins, 0.0);

    // define copies of correlation and total_weight for each of the threads
    int num_threads = atoi(std::getenv("OMP_NUM_THREADS"));
    std::vector<std::vector<double> > correlation_thread;
    std::vector<std::vector<double> > total_weight_thread;
    
    
    // compute the correlation as a function of parallel and perpendicular separation
    std::cout << "Computing the correlation as a function of parallel and perpendicular separation" << std::endl;
    // loop over plates
    #pragma omp parallel for schedule(dynamic)
    for (size_t plate_number = 0; plate_number < plates_list.size(); plate_number ++){
        // get the thread number
        int thread_num = omp_get_thread_num();
        double weight;
        
        // get the plate neighbours
        std::vector<int> plate_neighbours = kPlateNeighbours.GetNeighboursList(plates_list[plate_number]);
        
        // loop over spectra
        size_t number_of_spectra = (*spectra_list).num_objects_in_plate(plates_list[plate_number]);
        for (size_t spectrum_number = 0; spectrum_number < number_of_spectra; spectrum_number ++){
            
            LyaSpectrum lya_spectrum = (*spectra_list).list(plates_list[plate_number], spectrum_number);
            // checking that the object was loaded successfully
            if (lya_spectrum.dist() == _BAD_DATA_){
                if (flag_verbose_main >= 2){
                    #pragma omp critical (cout)
                    {
                        std::cout << "_BAD_DATA_ LyaSpectrum found. Ignoring..." << std::endl;
                    }
                }
                continue;
            }
            std::vector<LyaPixel> spectrum = lya_spectrum.spectrum();
        
        
            // loop over plates 2
            for (size_t plate_number2 = plate_number; plate_number2 < plates_list.size(); plate_number2 ++){
                
                // check that this plate is in the neighbours list
                bool is_neighbour = false;
                for (size_t index =  0; index < plate_neighbours.size(); index ++){
                    if (plate_neighbours[index] == plates_list[plate_number2]){
                        is_neighbour = true;
                        break;
                    }
                }
                if (not is_neighbour){
                    continue;
                }
                
                // loop over spectra 2
                size_t number_of_spectra2 = (*spectra_list).num_objects_in_plate(plates_list[plate_number2]);
                for (size_t spectrum_number2 = 0; spectrum_number2 < number_of_spectra; spectrum_number2 ++){
                    
                    // check that the two spectra are different
                    if (plate_number == plate_number2 and spectrum_number == spectrum_number2){
                        continue;
                    }
                    
                    LyaSpectrum lya_spectrum2 = (*spectra_list).list(plates_list[plate_number2], spectrum_number2);
                    // checking that the object was loaded successfully
                    if (lya_spectrum2.dist() == _BAD_DATA_){
                        if (flag_verbose_main >= 2){
                            #pragma omp critical (cout)
                            {
                                std::cout << "_BAD_DATA_ LyaSpectrum found. Ignoring..." << std::endl;
                            }
                        }
                        continue;
                    }
                    std::vector<LyaPixel> spectrum2 = lya_spectrum2.spectrum();
            
                    // compute angular separation
                    double cos_theta = lya_spectrum.angle().CosAngularDistance(lya_spectrum2.angle());
                    
                    double sigma_aux;
                    if (cos_theta == 1.0){
                        sigma_aux = 0.0;
                    }
                    else{
                        sigma_aux = 2.0*(1.0-cos_theta); // auxiliar variable to compute sigma values
                    }
                    
                    // check if the spectra are too far apart
                    double pair_min_sigma = sqrt(sigma_aux*spectrum[0].dist()*spectrum2[0].dist()); // minimum distance obtained for lowest redshift pixels
                    if (pair_min_sigma > max_sigma){ // if the minimum value for sigma (r_perp) is too large, the whole spectra is discarded
                        if (flag_verbose_main >= 2){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Spectrum rejected: min_sigma is too large" << std::endl;
                            }
                        }
                        continue;
                    }
                    
                    // compute the minimum distance between pairs of pixels
                    double pair_min_pi;
                    // case 1: first spectrum is at lower redshift than the second
                    if (spectrum[0].dist() <= spectrum2[0].dist()){
                        pair_min_pi = spectrum[0].dist() - spectrum2.back().dist();
                    }
                    // case 2: second spectrum is at lower redshift than the first
                    else{
                        pair_min_pi = spectrum2[0].dist() - spectrum.back().dist();
                    }
                    if (pair_min_pi > max_pi){ // if minimum value for pi (r_par) is too large, the whole spectra is discarded
                        if (flag_verbose_main >= 2){
                            #pragma omp critical (cout)
                            {
                                std::cout << "Spectrum rejected: min_pi is too small" << std::endl;
                            }
                        }
                        continue;
                    }
                    
                    // loop over pixel
                    for (size_t pixel_number = 0; pixel_number < spectrum.size(); pixel_number++){
                        
                        // cehck that the weight is not zero
                        if (spectrum[pixel_number].weight() == 0.0){
                            continue;
                        }
                        
                        // loop over pixel2
                        for (size_t pixel_number2 = 0; pixel_number2 < spectrum2.size(); pixel_number2++){
                            
                            
                            // cehck that the weight is not zero
                            if (spectrum2[pixel_number2].weight() == 0.0){
                                continue;
                            }
                            
                            // compute pi and sigma
                            double sigma = sqrt(sigma_aux*spectrum[pixel_number].dist()*spectrum2[pixel_number2].dist());
                            if (sigma > max_sigma){ // if sigma (r_perp) is too large, discard pixel
                                if (flag_verbose_main >= 3){
                                    #pragma omp critical (cout)
                                    {
                                        std::cout << "Pixel rejected: sigma is too large" << std::endl;
                                        std::cout << "sigma = " << sigma << std::endl;
                                    }
                                }
                                continue;
                            }
                            
                            double pi = spectrum[pixel_number].dist()-spectrum2[pixel_number2].dist();
                            if (pi < 0.0){
                                pi = -pi;
                            }
                            if ((pi > max_pi) or (pi <- max_pi)){ // if pi (r_par) is too large or too small, discard pixel
                                if (flag_verbose_main >= 3){
                                    #pragma omp critical (cout)
                                    {
                                        std::cout << "Pixel rejected: abs(pi) is too large" << std::endl;
                                        std::cout << "pi = " << pi << std::endl;
                                    }
                                }
                                continue;
                            }
                            
                            // locate pi pixel (i_index)
                            int i_index = int(pi/step_pi);
                            if (i_index < 0){
                                if (flag_verbose_main >= 2){
                                    #pragma omp critical (cout)
                                    {
                                        std::cout << "Pixel rejected: Bad indexing: i_index = " << i_index << std::endl;
                                    }
                                }
                                continue;
                            }
                            
                            // locate sigma pixel (j_index)
                            int j_index = int(sigma/step_sigma);
                            if (j_index < 0){
                                if (flag_verbose_main >= 2){
                                    #pragma omp critical (cout)
                                    {
                                        std::cout << "Pixel rejected: Bad indexing: j_index = " << j_index << std::endl;
                                    }
                                }
                                continue;
                            }
                            
                            // locate xi pixel (k)
                            int k_index = i_index*num_sigma_bins+j_index;
                            if (k_index < 0){
                                if (flag_verbose_main >= 2){
                                    #pragma omp critical (cout)
                                    {
                                        std::cout << "Pixel rejected: Bad indexing: k_index = " << k_index << std::endl;
                                    }
                                }
                                continue;
                                
                            }

                            // add pair contribution
                            weight = spectrum[pixel_number].weight()*spectrum2[pixel_number2].weight();
                            correlation_thread[thread_num][k_index] += spectrum[pixel_number].delta()*spectrum2[pixel_number2].delta()*weight;
                            total_weight_thread[thread_num][k_index] += weight;
                            
                        }
                    }
                }
            }
        }
    }

    // combine the measurements from the different threads
    for (size_t index = 0; index < num_bins; index ++){
        for (size_t thread_num = 0; thread_num < correlation_thread.size(); thread_num ++){
            correlation[index] += correlation_thread[thread_num][index];
            total_weight[index] += total_weight_thread[thread_num][index];
        }
    }
    
    // normalize the correlation
    for (int index =  0; index < num_bins; index ++){
        // check that the weight is not zero
        if (total_weight[index] == 0.0){
            std::cerr << "Warning : In compute_lya_1d : Index " << index;
            std::cerr << " shows zero weight. Consider increasing the bin size" << std::endl;
        }
        else{
            correlation[index] /= total_weight[index];
        }
    }
    
    // save results
    std::string filename = input.lya_auto_correlation_3d();
    std::cout << "Saving results to " << filename << std::endl;
    std::ofstream file;
    file.open(filename.c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
    if (file.is_open()){
        for (int index =  0; index < num_bins; index ++){
            file << index << " " << correlation[index] << std::endl;
        }
        
        file.close();
    }
    else{
        std::cerr << "Error : In compute_lya_3d : Unable to open file:" << std::endl << filename << std::endl;
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