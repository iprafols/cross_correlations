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
    plots_dir_ = input.plots();
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
        file << "	python $<" << std::endl;
        file.close();
    }
    else{
        std::cout << "Unable to open file:" << std::endl << filename << ".py" << std::endl;
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
            script << "import numpy as np" << std::endl;
            script << "import matplotlib.pyplot as plt" << std::endl;
            script << "import matplotlib.colors" << std::endl;
            script << "from matplotlib.colors import colorConverter" << std::endl;
            script << std::endl;
            script << "import matplotlib.ticker as ax" << std::endl;
            script << "from matplotlib.ticker import MultipleLocator, FormatStrFormatter, NullFormatter, ScalarFormatter" << std::endl;
            script << std::endl;
            script << "import math" << std::endl;
            script << "from math import sqrt" << std::endl;
            script << "\"\"\"" << std::endl;
            script << "EXPLANATION:" << std::endl;
            script << "    Plots the measured correlation function" << std::endl;
            script << "\"\"\"" << std::endl;
            script << "# loading variables" << std::endl;
            script << "max_sigma = " << input.max_sigma() << std::endl;
            script << "sigma_step = " << input.step_sigma() << std::endl;
            script << "num_sigma_bins = " << input.num_sigma_bins() << std::endl;
            script << "filename = '../" << output_base_name_ << ".full.data'" << std::endl;
            script << "cov_filename = '../" << output_base_name_ << ".bootstrap.diag.cov'" << std::endl;
            script << "data = np.genfromtxt(filename, names = True)" << std::endl;
            script << "cov = np.genfromtxt(cov_filename, usecols = (1,2))" << std::endl;
            script << "data.sort(order = 'bin_index')" << std::endl;
            script << "cov.sort(axis = 0)" << std::endl;
            script << std::endl;
            script << "# plotting" << std::endl;
            script << "for j in range (0, num_sigma_bins):" << std::endl;
            script << "    pi = [value for i,value in enumerate(data['mean_pi']) if (i % num_sigma_bins == j)]" << std::endl;
            script << "    xi = [value for i,value in enumerate(data['xi']) if (i % num_sigma_bins == j)]" << std::endl;
            script << "    dxi = [sqrt(value) for i,value in enumerate(cov[:,1]) if (i % num_sigma_bins == j)]" << std::endl;
            script << "    fig = plt.figure(figsize=(14,7))" << std::endl;
            script << "    ax = fig.add_subplot(1,1,1)" << std::endl;
            script << "    ax.set_xlabel(r'$\\pi\\left(h^{-1}Mpc\\right)$', fontsize = 20)" << std::endl;
            script << "    ax.set_ylabel(r'$\\xi\\left(\\pi, \\sigma\\right)$', fontsize = 20)" << std::endl;
            script << "    ax.errorbar(pi,xi,yerr = dxi, fmt = 'b.')" << std::endl;
            script << "    ax.text(0.05,0.05,r'$\\sigma/\\left(h^{-1}Mpc\\right)\\in\\left[' + str(j*sigma_step) + ', '+str((j+1)*sigma_step )+ r'\\right) $',transform = ax.transAxes, fontsize = 20)" << std::endl;
            script << "    fig.savefig('correlation_measurements_sigma_bin_' + str(j) + '.eps')" << std::endl;
            script << "    plt.close(fig)" << std::endl;
            
            script.close();
        }
        else{
            std::cout << "Unable to open file:" << std::endl << filename << ".py" << std::endl;
        }
    }
    

}

void PlotsObject::PlotRADECDispersion(Dataset& dataset, const bool update_script) const{
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
        std::cout << "Unable to open file:" << std::endl << filename << ".dat" << std::endl;
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
            std::cout << "Unable to open file:" << std::endl << filename << ".py" << std::endl;
        }
    }
    
}

void PlotsObject::PlotZHistogram(Dataset& dataset, const bool update_script) const{
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
        std::cout << "Unable to open file:" << std::endl << filename << ".dat" << std::endl;
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
            std::cout << "Unable to open file:" << std::endl << filename << ".py" << std::endl;
        }

    }
}




