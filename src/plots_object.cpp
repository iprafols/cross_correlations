/**
 plots.cpp
 Purpose: This files contains the body for the functions defined in plots.h
 
 @author Ignasi Pérez-Ràfols
 @version 1.0 09/25/2014
 */

#include "plots_object.h"

PlotsObject::PlotsObject(const Input& input){
    /**
     EXPLANATION:
     Cosntructs a Plots instance
     
     INPUTS:
     input - a Input instance
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Plots
     
     FUNCITONS USED:
     NONE
     */
    plots_dir_ = input.output() + input.plots();
    output_base_name_ = input.output_base_name();
    MakePlottingMakefile();
}

void PlotsObject::MakePlottingMakefile() const{
    /*
     EXPLANATION:
     Writes a makefile for plotting
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Plots
     
     FUNCITONS USED:
     NONE
     */
    
    // set the name of the makefile
    std::string filename = plots_dir_ + "makefile";
    
    std::ofstream file;
    file.open(filename.c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
    if (file.is_open()){
        file << "# lines starting with '#' are comments" << std::endl;
        file << "# usage:    make 'target'" << std::endl;
        file << "# 'target' specifies the task to be performed. Default is all" << std::endl;
        file << std::endl;        
        file << "#" << std::endl;
        file << "# VARIABLES:" << std::endl;
        file << "#" << std::endl;
        file << std::endl;
        file << "# this variable sets the compiler to use" << std::endl;
        file << "CC = g++" << std::endl;
        file << "# this variable contains the list of sources" << std::endl;
        file << "SOURCES = $(wildcard *.py)" << std::endl;
        file << "# this variable contains the list of object files" << std::endl;
        file << "OBJECTS = $(patsubst %.py,%.pyc,$(SOURCES))" << std::endl;
        file << std::endl;
        file << "#" << std::endl;
        file << "# TARGETS:" << std::endl;
        file << "#" << std::endl;
        file << std::endl;
        file << "all: $(OBJECTS)" << std::endl;
        file << std::endl;
        file << "%.pyc: %.py" << std::endl;
        file << "\t" << "python $<" << std::endl;
        file.close();
    }
    else{
        std::cout << "Error :  In PlotsObject::MakePlottingMakefile : Unable to open file:" << std::endl << filename << ".py" << std::endl;
    }    
}

void PlotsObject::PlotCrossCorrelation(const CorrelationResults& res, const Input& input, const bool update_script) const{
    /*
     EXPLANATION:
     Plots the correlation function
     
     INPUTS:
     res - object where the results are stored
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     Plots
     
     FUNCITONS USED:
     NONE
    */
    
    // set the name of the results and script file
    std::string filename = "correlation_measurements";
        
    // checking if the script has to be rewriten
    if (update_script){
        // opening file
        std::ofstream script;
        script.open((plots_dir_ + filename + ".py").c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
        if (script.is_open()){
            script << "\"\"\"" << std::endl;
            script << "EXPLANATION:" << std::endl;
            script << "    Plots the measured correlation function" << std::endl;
            script << "\"\"\"" << std::endl;
            script << std::endl;
            script << "import sys" << std::endl;
            script << "sys.path.append('/Users/iprafols/software/development/cross_correlations/lib/')" << std::endl;
            script << std::endl;
            script << "import correlation_process as cp" << std::endl;
            script << "import numpy as np" << std::endl;

            script << std::endl;
            script << "def main():" << std::endl;
            script << "    filename = '../" << input.output_base_name() << "'" << std::endl;
            script << std::endl;
            script << "    # uncomment this block to plot models, change model names to the appropiate names" << std::endl;
            script << "    '''" << std::endl;
            script << "    rmin = 10" << std::endl;
            script << "    model_filename = '../fit/DLA_DR12_bin2Mpc_80Mpc_sampleA_rmin10_0Mpc_rmax90_0Mpc_residuals.dat'" << std::endl;
            script << std::endl;
            script << "    rmin2 = 5" << std::endl;
            script << "    model_filename2 = '../fit/DLA_DR12_bin2Mpc_80Mpc_sampleA_rmin5_0Mpc_rmax90_0Mpc_residuals.dat'" << std::endl;
            script << "    '''" << std::endl;
            script << std::endl;
            script << "    # define pi and sigma arrays" << std::endl;
            script << "    pi  = np.arange(-" << input.max_pi() << ", " << input.max_pi() + input.step_pi() << ", " << input.step_pi() << ", dtype=float)" << std::endl;
            script << "    sigma = np.arange(0, " << input.max_sigma() + input.step_sigma() << ", " << input.step_sigma() << ", dtype=float)" << std::endl;
            script << std::endl;
            script << "    # define a more convenient array to plot rebinned plots" << std::endl;
            script << "    pi_new = np.array([-80, -72, -64, -56, -48, -40, -32, -28, -20, -16, -12, -8, -4, 0," << std::endl;
            script << "                       4, 8, 12, 16, 20, 28, 32, 40, 48, 56, 64, 72, 80], dtype=float)" << std::endl;
            script << "    sigma_new = np.array([0, 2, 4, 8, 16, 24, 32, 40, 80], dtype=float)" << std::endl;
            script << std::endl;
            script << "    # load data" << std::endl;
            script << "    print 'loading data from {}'.format(filename)" << std::endl;
            script << "    data = cp.CorrData(filename, pi, sigma)" << std::endl;
            script << std::endl;
            script << "    # uncomment this block to plot models" << std::endl;
            script << "    '''" << std::endl;
            script << "    # load models" << std::endl;
            script << "    print 'loading model from {}'.format(model_filename)" << std::endl;
            script << "    model = cp.CorrModel(model_filename)" << std::endl;
            script << "    print 'loading model from {}'.format(model_filename2)" << std::endl;
            script << "    model2 = cp.CorrModel(model_filename2)" << std::endl;
            script << "    '''" << std::endl;
            script << std::endl;
            script << "    # rebin data" << std::endl;
            script << "    print 'rebinning data from {}'.format(filename)" << std::endl;
            script << "    transformation_matrix = data.computeTransformationMatrix(pi_new, sigma_new, keep_matrix=True)" << std::endl;
            script << "    data.rebinData(transformation_matrix, keep_results=True, rebin_grid=True)" << std::endl;
            script << std::endl;
            script << "    # uncomment this block to plot models" << std::endl;
            script << "    '''" << std::endl;
            script << "    # rebin models" << std::endl;
            script << "    print 'rebinning model from {}'.format(model_filename)" << std::endl;
            script << "    model.rebinModel(pi_new, sigma_new, keep_results=True)" << std::endl;
            script << "    print 'rebinning model from {}'.format(model_filename2)" << std::endl;
            script << "    model2.rebinModel(pi_new, sigma_new, keep_results=True)" << std::endl;
            script << std::endl;
            script << "    # to plot models uncomment the first set of plotting functions and comment the second set; otherwise comment the first set and uncomment the second set" << std::endl;
            script << "    # plot input data" << std::endl;
            script << "    #cp.plot(data, './not_rebinned/', fmt_list='k.', sigma_bins=sigma, plot_rebinned_list=False, model_list = [model, model2], fmt_model_list = ['k-', 'r--'], plot_model_rebinned_list=False, base_fig_name='cross_correlation_not_rebinned', rmin_list=[rmin, rmin2])" << std::endl;
            script << "    cp.plot(data, './not_rebinned/', fmt_list='k.', sigma_bins=sigma, plot_rebinned_list=False, base_fig_name='cross_correlation_not_rebinned'" << std::endl;
            script << std::endl;
            script << "    # plot rebinned data" << std::endl;
            script << "    #cp.plot(data, './rebinned/', fmt_list='k.', sigma_bins=sigma_new, plot_rebinned_list=True, model_list = [model, model2], fmt_model_list = ['k-', 'r--'], plot_model_rebinned_list=True, rmin_list=[rmin, rmin2], single_plot=True)" << std::endl;
            script << "    cp.plot(data, './rebinned/', fmt_list='k.', sigma_bins=sigma_new, plot_rebinned_list=True, single_plot=True)" << std::endl;
            script << std::endl;
            script << "    # plot contours" << std::endl;
            script << "    #cp.plot(data, "./", model_list=model2, plot_rebinned_list=False, plot_model_rebinned_list=False, contour=True)" << std::endl;
            script << std::endl;
            script << std::endl;
            script << "if __name__ == '__main__':" << std::endl;
            script << "    main()" << std::endl;
            script << "" << std::endl;
            script.close();
        }
        else{
            std::cout << "Error : In PlotsObject::PlotCrossCorrelation : Unable to open file:" << std::endl << filename << ".py" << std::endl;
        }
    }
    

}

void PlotsObject::PlotLyaAuto(const std::vector<LyaAutoInterpolationMap>& lya_auto, const bool update_script) const{
    /**
     EXPLANATION:
     Plots the lyman-alpha autocorrelation
     
     INPUTS:
     NONE
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     LyaAutoInterpolationMap
     PlotsObject
     
     FUNCITONS USED:
     NONE
     */
    
    // set the name of the results and script file
    std::string filename =  "lya_autocorrelation";
    
    // open results file
    std::ofstream result;
    result.open((plots_dir_ + filename + ".dat").c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
    if (result.is_open()){
        
        result << "# data required to redo the '" << filename << ".eps' plot." << std::endl;
        result << "# n z lya_auto" << std::endl;
        
        // load the LyaInterpolationMap
        for (size_t i = 0; i < lya_auto.size(); i ++){
            std::map<double,double> interpolation_map = lya_auto[i].interpolation_map();
            
            for (std::map<double,double>::const_iterator it = interpolation_map.begin(); it != interpolation_map.end(); it ++){
                result << i << " " << (*it).first << " " << (*it).second << std::endl;
            }
        }
        
        result.close();
    }
    else{
        std::cout << "Error : In PlotsObject::PlotLyaAuto : Unable to open file:" << std::endl << filename << ".dat" << std::endl;
    }
    
    // checking if the script has to be rewriten
    if (update_script){
        // open script file
        std::ofstream script;
        script.open((plots_dir_ + filename + ".py").c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
        if (script.is_open()){
            script << "import numpy as np" << std::endl;
            script << "import math" << std::endl;
            script << "from math import acos" << std::endl;
            script << "import matplotlib.pyplot as plt" << std::endl;
            script << "import matplotlib.colors" << std::endl;
            script << "from matplotlib.colors import colorConverter" << std::endl;
            script << std::endl;
            script << "import matplotlib.ticker as ax" << std::endl;
            script << "from matplotlib.ticker import MultipleLocator, FormatStrFormatter, NullFormatter, ScalarFormatter" << std::endl;
            script << "\"\"\"" << std::endl;
            script << "EXPLANATION:" << std::endl;
            script << "    Plots the lya autocorrelation" << std::endl;
            script << "\"\"\"" << std::endl;
            script << "# loading variables" << std::endl;
            script << "filename = '" << filename << ".dat'" << std::endl;
            script << "data = np.genfromtxt(filename, names = True, skip_header = 1)" << std::endl;
            script << std::endl;
            script << "max_n = np.amax(data['n'])" << std::endl;
            script << std::endl;
            script << "# plotting lya autocorrelation vs z" << std::endl;
            script << "fig = plt.figure(figsize=(18, 9))" << std::endl;
            script << "ax = fig.add_subplot(1, 1, 1)" << std::endl;
            script << "ax.tick_params(axis='both', which='major', labelsize=35, length=6, width=2)" << std::endl;
            script << "ax.tick_params(axis='both', which='minor', labelsize=35, length=4, width=1)" << std::endl;
            script << "ax.set_xlabel(r'$z$', fontsize=35)" << std::endl;
            script << "ax.set_ylabel(r'$\\xi$', fontsize=35)" << std::endl;
            script << "for aux_n in range(0, int(max_n)):" << std::endl;
            script << "    data_z = [z for (z, n) in zip(data['z'], data['n']) if n == aux_n]" << std::endl;
            script << "    data_xi = [xi for (xi, n) in zip(data['lya_auto'], data['n']) if n == aux_n]" << std::endl;
            script << "    ax.plot(data_z ,data_xi, '-', label='n = '+str(aux_n))" << std::endl;
            script << "ax.legend(loc=0, numpoints=1, prop={'size':20}, frameon=False)" << std::endl;
            script << "plt.tight_layout()" << std::endl;
            script << "fig.savefig('" << filename << "_vs_z.eps')" << std::endl;
            script << std::endl;
            script << "# plotting lya autocorrelation vs n" << std::endl;
            script << "fig2 = plt.figure(figsize=(18, 9))" << std::endl;
            script << "ax2 = fig2.add_subplot(1, 1, 1)" << std::endl;
            script << "ax2.tick_params(axis='both', which='major', labelsize=35, length=6, width=2)" << std::endl;
            script << "ax2.tick_params(axis='both', which='minor', labelsize=35, length=4, width=1)" << std::endl;
            script << "ax2.set_xlabel(r'$n$', fontsize=35)" << std::endl;
            script << "ax2.set_ylabel(r'$\\xi$', fontsize=35)" << std::endl;
            script << "for aux_z in data_z:" << std::endl;
            script << "    data_n = [n for (z, n) in zip(data['z'], data['n']) if z == aux_z]" << std::endl;
            script << "    data_xi = [xi for (xi, z) in zip(data['lya_auto'], data['z']) if z == aux_z]" << std::endl;
            script << "    ax2.plot(data_n ,data_xi, '-', label='z = '+str(aux_z))" << std::endl;
            script << "ax2.legend(loc=0, numpoints=1, prop={'size':20}, frameon=False)" << std::endl;
            script << "plt.tight_layout()" << std::endl;
            script << "fig2.savefig('" << filename << "_vs_n.eps')" << std::endl;
            script.close();
        }
        else{
            std::cout << "Error : In PlotsObject::PlotLyaAuto : Unable to open file:" << std::endl << filename << ".py" << std::endl;
        }
    }
    
}

void PlotsObject::PlotRADECDispersion(const Dataset& dataset, const bool update_script) const{
    /**
     EXPLANATION:
     Plots the RA-DEC dispersion for the given objects
     
     INPUTS:
     dataset - a Dataset instance from which to plot
     update_script - a boolean; if true, python script is updated, if false, it is not - defaut = False
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     Plots
     
     FUNCITONS USED:
     NONE
     */
    
    // set the name of the results and script file
    std::string filename =  dataset.name() + "_RA_DEC_dispersion";
    
    // open results file
    std::ofstream result;
    result.open((plots_dir_ + filename + ".dat").c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
    if (result.is_open()){
        
        result << "# data required to redo the '" << filename << ".png' plot. All angles are in radians" << std::endl;
        result << "# RA DEC" << std::endl;
        
        // load the AstroObject pointers
        dataset.GiveRADEC(result);                
        result.close();
    }
    else{
        std::cout << "Error : In PlotsObject::PlotRADECDispersion : Unable to open file:" << std::endl << filename << ".dat" << std::endl;
    }
    
    // checking if the script has to be rewriten
    if (update_script){
        // open script file
        std::ofstream script;
        script.open((plots_dir_ + filename + ".py").c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
        if (script.is_open()){
            script << "import numpy as np" << std::endl;
            script << "import math" << std::endl;
            script << "from math import acos" << std::endl;
            script << "import matplotlib.pyplot as plt" << std::endl;
            script << "import matplotlib.colors" << std::endl;
            script << "from matplotlib.colors import colorConverter" << std::endl;
            script << std::endl;
            script << "import matplotlib.ticker as ax" << std::endl;
            script << "from matplotlib.ticker import MultipleLocator, FormatStrFormatter, NullFormatter, ScalarFormatter" << std::endl;
            script << "\"\"\"" << std::endl;
            script << "EXPLANATION:" << std::endl;
            script << "    Plots the RA - DEC dispersion of the " << dataset.name() << " sample" << std::endl;
            script << "\"\"\"" << std::endl;
            script << "# loading variables" << std::endl;
            script << "filename = '" << filename << ".dat'" << std::endl;
            script << "data = np.genfromtxt(filename, names = True, skip_header = 1)" << std::endl;
            script << std::endl;
            script << "for i in range(0,len(data['RA'])):" << std::endl;
            script << "    if data['RA'][i] < 1.0:" << std::endl;
            script << "        data['RA'][i] += 2.0*acos(-1.0)" << std::endl;
            script << std::endl;
            script << "# plotting RA/DEC dispersion" << std::endl;
            script << "fig = plt.figure(figsize=(18,9))" << std::endl;
            script << "ax = fig.add_subplot(1,1,1)" << std::endl;
            script << "ax.tick_params(axis='both',which='major',labelsize=35,length=6,width=2)" << std::endl;
            script << "ax.tick_params(axis='both',which='minor',labelsize=35,length=4,width=1)" << std::endl;
            script << "#ax.set_xticklabels([' ',r'$0$',r'$50$',r'$100$',r'$150$',r'$200$',r'$250$',r'$300$'])" << std::endl;
            script << "#ax.set_yticklabels([' ',r'$0$',r'$10$',r'$20$',r'$30$',r'$40$',r'$50$',r'$60$',r'$70$'])" << std::endl;
            script << "ax.set_xlabel('J2000 RA (rad)',fontsize=35)" << std::endl;
            script << "ax.set_ylabel('J2000 DEC (rad)',fontsize=35)" << std::endl;
            script << "ax.set_xlim(1,8)" << std::endl;
            script << "ax.set_ylim(-0.2,1.2)" << std::endl;
            script << "ax.plot(data['RA'],data['DEC'],'k.',alpha=0.5)" << std::endl;
            script << "fig.savefig('" << filename << ".png')" << std::endl;
            
            script.close();
        }
        else{
            std::cout << "Error : In PlotsObject::PlotRADECDispersion : Unable to open file:" << std::endl << filename << ".py" << std::endl;
        }
    }
    
}

void PlotsObject::PlotZHistogram(const Dataset& dataset, const bool update_script) const{
    /**
     EXPLANATION:
     Plots the redshift histogram for the given objects
     
     INPUTS:
     dataset - a Dataset instance from which to plot
     update_script - a boolean; if true, python script is updated, if false, it is not - defaut = False
     
     OUTPUTS:
     NONE
     
     CLASSES USED:
     AstroObject
     Plots
     
     FUNCITONS USED:
     NONE
     */
    
    // set the name of the results and script file
    std::string filename =  dataset.name() + "_z_histogram";
    
    // open results file
    std::ofstream result;
    result.open((plots_dir_ + filename + ".dat").c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
    if (result.is_open()){
        
        result << "# data required to redo the '" << filename << ".png' plot" << std::endl;
        result << "# z" << std::endl;
        
        // load the AstroObject pointers
        dataset.GiveZ(result);                
        result.close();
    }
    else{
        std::cout << "Error : In PlotsObject::PlotZHistogram : Unable to open file:" << std::endl << filename << ".dat" << std::endl;
    }
    
    // checking if the script has to be rewriten
    if (update_script){
        // open script file
        std::ofstream script;
        script.open((plots_dir_ + filename + ".py").c_str(),std::ofstream::trunc); // opens the file erasing the previous contents
        if (script.is_open()){
            script << "import numpy as np" << std::endl;
            script << "import matplotlib.pyplot as plt" << std::endl;
            script << "import matplotlib.colors" << std::endl;
            script << "from matplotlib.colors import colorConverter" << std::endl;
            script << std::endl;
            script << "import matplotlib.ticker as ax" << std::endl;
            script << "from matplotlib.ticker import MultipleLocator, FormatStrFormatter, NullFormatter, ScalarFormatter" << std::endl;
            script << "\"\"\"" << std::endl;
            script << "EXPLANATION:" << std::endl;
            script << "    Plots the redshift histogram of the " << dataset.name() << " sample" << std::endl;
            script << "\"\"\"" << std::endl;
            script << "# loading variabes" << std::endl;
            script << "filename = '" << filename << ".dat'" << std::endl;
            script << "data = np.genfromtxt(filename)" << std::endl;
            script << std::endl;
            script << "# ploting redshift histogram" << std::endl;
            script << "fig = plt.figure(figsize=(18,9))" << std::endl;
            script << "ax = fig.add_subplot(1,1,1)" << std::endl;
            script << "ax.tick_params(axis='both',which='major',labelsize=35,length=6,width=2)" << std::endl;
            script << "ax.tick_params(axis='both',which='minor',labelsize=35,length=4,width=1)" << std::endl;
            script << "#ax.set_xticklabels(['$2.0$',r'$2.2$',r'$2.4$',r'$2.6$',r'$2.8$',r'$3.0$',r'$3.2$',r'$3.4$',r'$3.6$'])" << std::endl;
            script << "#ax.set_yticklabels([' ',r'$1000$',r'$2000$',r'$3000$',r'$4000$',r'$5000$',r'$6000$',r'$7000$',r'$8000$'])" << std::endl;
            script << "ax.set_xlabel('z',fontsize=35)" << std::endl;
            script << "ax.set_ylabel('number of objects',fontsize=20)" << std::endl;
            script << "ax.hist(data,50,histtype='step',color='black')" << std::endl;
            script << "fig.savefig('" << filename << ".png')" << std::endl;
            
            script.close();
        }
        else{
            std::cout << "Error : In PlotsObject::PlotZHistogram : Unable to open file:" << std::endl << filename << ".py" << std::endl;
        }

    }
}




