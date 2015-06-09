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
            script << "import numpy as np" << std::endl;
            script << "import matplotlib.pyplot as plt" << std::endl;
            script << "import matplotlib.colors" << std::endl;
            script << "from matplotlib.colors import colorConverter" << std::endl;
            script << std::endl;
            script << "import matplotlib.ticker as ax" << std::endl;
            script << "from matplotlib.ticker import MultipleLocator, FormatStrFormatter, NullFormatter, ScalarFormatter" << std::endl;
            script << std::endl;
            script << "\"\"\"" << std::endl;
            script << "EXPLANATION:" << std::endl;
            script << "    Plots the measured correlation function" << std::endl;
            script << "\"\"\"" << std::endl;
            script << "# loading variables" << std::endl;
            script << "max_sigma = " << input.max_sigma() << std::endl;
            script << "sigma_step = " << input.step_sigma() << std::endl;
            script << "num_sigma_bins = " << input.num_sigma_bins() << std::endl;
            script << "num_pi_bins = " << input.num_pi_bins() << std::endl;
            script << std::endl;
            script << "filename = '../" << input.fit() << output_base_name_ << "_residuals.dat'" << std::endl;
            script << "data = np.genfromtxt(filename)" << std::endl;
            script << "filename = '../" << input.best_fit() << output_base_name_ << "_residuals.dat'" << std::endl;
            script << "data_all = np.genfromtxt(filename)" << std::endl;
            script << "filename = '../" << input.best_fit() << output_base_name_ << "_altresiduals.dat'" << std::endl;
            script << "data_iso = np.genfromtxt(filename)" << std::endl;
                        script << std::endl;
            script << "# plotting sigma bins" << std::endl;
            script << "for j in range (0, num_sigma_bins):" << std::endl;
            script << "    pi = [value for i,value in zip(data_all[:,0],data_all[:,1]) if (i % num_sigma_bins == j)]" << std::endl;
            script << "    xi = [value for i,value in zip(data_all[:,0],data_all[:,8]) if (i % num_sigma_bins == j)]" << std::endl;
            script << "    dxi = [value for i,value in zip(data_all[:,0],data_all[:,9]) if (i % num_sigma_bins == j)]" << std::endl;
            script << "    model_xi = [value for i,value in zip(data_all[:,0],data_all[:,7]) if (i % num_sigma_bins == j)]" << std::endl;
            script << std::endl;
            script << "    pi_fit_pos = [value for i,value in zip(data[:,0],data[:,1]) if (i % num_sigma_bins == j) and value >= -sigma_step]" << std::endl;
            script << "    model_xi_fit_pos = [value for i,value,r_par in zip(data[:,0],data[:,7],data[:,1]) if (i % num_sigma_bins == j) and r_par >= -sigma_step]" << std::endl;
            script << "    pi_fit_neg = [value for i,value in zip(data[:,0],data[:,1]) if (i % num_sigma_bins == j) and value < 0]" << std::endl;
            script << "    model_xi_fit_neg = [value for i,value,r_par in zip(data[:,0],data[:,7],data[:,1]) if (i % num_sigma_bins == j) and r_par < 0]" << std::endl;
            script << std::endl;
            script << "    fig = plt.figure(figsize=(14,7))" << std::endl;
            script << "    ax = fig.add_subplot(1,1,1)" << std::endl;
            script << "    ax.set_xlabel(r'$\\pi\\left(h^{-1}Mpc\\right)$', fontsize = 20)" << std::endl;
            script << "    ax.set_ylabel(r'$\\xi\\left(\\pi, \\sigma\\right)$', fontsize = 20)" << std::endl;
            script << "    ax.errorbar(pi,xi,yerr = dxi, fmt = 'b.', label = 'data')" << std::endl;
            script << "    ax.plot(pi, model_xi, 'r--', label = 'model prediction')" << std::endl;
            script << "    ax.plot(pi_fit_pos, model_xi_fit_pos, 'b-', label = 'best-fit model')" << std::endl;
            script << "    ax.plot(pi_fit_neg, model_xi_fit_neg, 'b-')" << std::endl;
            script << "    ax.text(0.05,0.05,r'$\\sigma/\\left(h^{-1}Mpc\\right)\\in\\left[' + str(j*sigma_step) + ', '+str((j+1)*sigma_step )+ r'\\right) $',transform = ax.transAxes, fontsize = 20)" << std::endl;
            script << "    ax.legend(numpoints = 1, loc = 4, prop = {'size':20}, frameon = False)" << std::endl;
            script << "    fig.savefig('correlation_measurements_sigma_bin_' + str(j) + '.eps')" << std::endl;
            script << "    plt.close(fig)" << std::endl;
            script << std::endl;
            script << "    # plotting 2D contours" << std::endl;
            script << "    pi = [value for i,value in zip(data_all[:,0],data_all[:,1]) if (i % num_sigma_bins == 0)]" << std::endl;
            script << "    sigma = [value for i,value in zip(data_all[:,0],data_all[:,2]) if (int(i) / int(num_sigma_bins) == 0)]" << std::endl;
            script << std::endl;  
            script << "    r2xi = np.reshape(data_all[:,8]*data_all[:,4]*data_all[:,4],(num_pi_bins, num_sigma_bins))" << std::endl;
            script << "    dxi = np.reshape(data_all[:,9]*data_all[:,4]*data_all[:,4],(num_pi_bins, num_sigma_bins))" << std::endl;
            script << "    model_r2xi = np.reshape(data_all[:,7]*data_all[:,4]*data_all[:,4],(num_pi_bins, num_sigma_bins))" << std::endl;
            script << "    model_r2xi_iso = np.reshape(data_iso[:,7]*data_iso[:,4]*data_iso[:,4],(num_pi_bins, num_sigma_bins))" << std::endl;
            script << "    model_r2xi_broaddand = model_r2xi - model_r2xi_iso" << std::endl;
            script << "    residual = (r2xi-model_r2xi)/dxi*(r2xi-model_r2xi)/dxi" << std::endl;
            script << std::endl;
            script << "    # fiitng range" << std::endl;
            script << "    rmin = 40" << std::endl;
            script << "    rmax = 180" << std::endl;
            script << "    inner_circle_sigma = np.arange(0,rmin,0.001)" << std::endl;
            script << "    inner_circle_pi = np.sqrt(rmin*rmin-inner_circle_sigma*inner_circle_sigma)" << std::endl;
            script << "    outer_circle_sigma = np.arange(0,rmax,0.001)" << std::endl;
            script << "    outer_circle_pi = np.sqrt(rmax*rmax-outer_circle_sigma*outer_circle_sigma)" << std::endl;
            script << std::endl;
            script << "    # contour settings" << std::endl;
            script << "    cmap = plt.cm.get_cmap('RdYlBu')" << std::endl;
            script << "    vmin = min(np.amin(r2xi),np.amin(model_r2xi))" << std::endl;
            script << "    vmax = max(np.amax(r2xi),np.amax(model_r2xi))" << std::endl;
            script << "    vmin = -60.0" << std::endl;
            script << "    vmax = 60.0#np.amax(model_r2xi)" << std::endl;
            script << "    num_colors = 20" << std::endl;
            script << std::endl;
            script << "    levels = np.arange(vmin,vmax,(vmax-vmin)/num_colors)" << std::endl;
            script << "    residuals_levels = np.arange(0,vmax/6,vmax/6/num_colors)" << std::endl;
            script << std::endl;
            script << "    # plotting" << std::endl;
            script << "    fig = plt.figure(figsize = (28,7))" << std::endl;
            script << std::endl;
            script << "    ax = fig.add_subplot(1,5,1)" << std::endl;
            script << "    ax.set_xlim(np.amin(sigma),np.amax(sigma))" << std::endl;
            script << "    ax.set_ylim(np.amin(pi),np.amax(sigma))" << std::endl;
            script << "    ax.set_xlabel(r'$\\sigma\\left(h^{-1}Mpc\\right)$', fontsize = 20)" << std::endl;
            script << "    ax.set_ylabel(r'$\\pi\\left(h^{-1}Mpc\\right)$', fontsize = 20)" << std::endl;
            script << "    ax.set_title(r'$r^{2}\\xi\\left(h^{-1}Mpc^{2}\\right)$', fontsize = 20)" << std::endl;
            script << "    cs = ax.contourf(sigma, pi, r2xi, levels, cmap = cmap, vmin = vmin, vmax = vmax)" << std::endl;
            script << "    fig.colorbar(cs, ax = ax, shrink = 0.9)" << std::endl;
            script << std::endl;
            script << "    ax2 = fig.add_subplot(1,5,2)" << std::endl;
            script << "    ax2.set_xlim(np.amin(sigma),np.amax(sigma))" << std::endl;
            script << "    ax2.set_ylim(np.amin(pi),np.amax(sigma))" << std::endl;
            script << "    ax2.set_xlabel(r'$\\sigma\\left(h^{-1}Mpc\\right)$', fontsize = 20)" << std::endl;
            script << "    ax2.set_title(r'$r^{2}\\xi_{\\rm model}\\left(h^{-1}Mpc^{2}\\right)$', fontsize = 20)" << std::endl;
            script << "    cs2 = ax2.contourf(sigma, pi, model_r2xi, levels, cmap = cmap, vmin = vmin, vmax = vmax)" << std::endl;
            script << "    fig.colorbar(cs2, ax = ax2, shrink = 0.9)" << std::endl;
            script << std::endl;
            script << "    ax3 = fig.add_subplot(1,5,3)" << std::endl;
            script << "    ax3.set_xlim(np.amin(sigma),np.amax(sigma))" << std::endl;
            script << "    ax3.set_ylim(np.amin(pi),np.amax(sigma))" << std::endl;
            script << "    ax3.set_xlabel(r'$\\sigma\\left(h^{-1}Mpc\\right)$', fontsize = 20)" << std::endl;
            script << "    ax3.set_title(r'$r^{2}\\xi_{\\rm cosmo}\\left(h^{-1}Mpc^{2}\\right)$', fontsize = 20)" << std::endl;
            script << "    cs3 = ax3.contourf(sigma, pi, model_r2xi_iso, levels, cmap = cmap, vmin = vmin, vmax = vmax)" << std::endl;
            script << "    fig.colorbar(cs3, ax = ax3, shrink = 0.9)" << std::endl;
            script << std::endl;
            script << "    ax4 = fig.add_subplot(1,5,4)" << std::endl;
            script << "    ax4.set_xlim(np.amin(sigma),np.amax(sigma))" << std::endl;
            script << "    ax4.set_ylim(np.amin(pi),np.amax(sigma))" << std::endl;
            script << "    ax4.set_xlabel(r'$\\sigma\\left(h^{-1}Mpc\\right)$', fontsize = 20)" << std::endl;
            script << "    ax4.set_title(r'$r^{2}\\xi_{\\rm broadband}\\left(h^{-1}Mpc^{2}\\right)$', fontsize = 20)" << std::endl;
            script << "    cs4 = ax4.contourf(sigma, pi, model_r2xi_broaddand, levels, cmap = cmap, vmin = vmin, vmax = vmax)" << std::endl;
            script << "    fig.colorbar(cs4, ax = ax4, shrink = 0.9)" << std::endl;
            script << std::endl;
            script << "    ax5 = fig.add_subplot(1,5,5)" << std::endl;
            script << "    ax5.set_xlim(np.amin(sigma),np.amax(sigma))" << std::endl;
            script << "    ax5.set_ylim(np.amin(pi),np.amax(sigma))" << std::endl;
            script << "    ax5.set_xlabel(r'$\\sigma\\left(h^{-1}Mpc\\right)$', fontsize = 20)" << std::endl;
            script << "    ax5.set_title(r'$\\left(\\xi-\\xi_{\\rm model}\\right)^{2}/C$', fontsize = 20)" << std::endl;
            script << "    cs5 = ax5.contourf(sigma, pi, residual, residuals_levels, cmap = cmap, vmin = 0.0, vmax = vmax/6, extend = 'neither')" << std::endl;
            script << "    fig.colorbar(cs5, ax = ax5, shrink = 0.9)" << std::endl;
            script << std::endl;
            script << "    # plotting fiting range" << std::endl;
            script << "    ax.plot(inner_circle_sigma,inner_circle_pi,'k-')" << std::endl;
            script << "    ax.plot(inner_circle_sigma,-inner_circle_pi,'k-')" << std::endl;
            script << "    ax.plot(outer_circle_sigma,outer_circle_pi,'k-')" << std::endl;
            script << "    ax.plot(outer_circle_sigma,-outer_circle_pi,'k-')" << std::endl;
            script << "    ax2.plot(inner_circle_sigma,inner_circle_pi,'k-')" << std::endl;
            script << "    ax2.plot(inner_circle_sigma,-inner_circle_pi,'k-')" << std::endl;
            script << "    ax2.plot(outer_circle_sigma,outer_circle_pi,'k-')" << std::endl;
            script << "    ax2.plot(outer_circle_sigma,-outer_circle_pi,'k-')" << std::endl;
            script << "    ax3.plot(inner_circle_sigma,inner_circle_pi,'k-')" << std::endl;
            script << "    ax3.plot(inner_circle_sigma,-inner_circle_pi,'k-')" << std::endl;
            script << "    ax3.plot(outer_circle_sigma,outer_circle_pi,'k-')" << std::endl;
            script << "    ax3.plot(outer_circle_sigma,-outer_circle_pi,'k-')" << std::endl;
            script << "    ax4.plot(inner_circle_sigma,inner_circle_pi,'k-')" << std::endl;
            script << "    ax4.plot(inner_circle_sigma,-inner_circle_pi,'k-')" << std::endl;
            script << "    ax4.plot(outer_circle_sigma,outer_circle_pi,'k-')" << std::endl;
            script << "    ax4.plot(outer_circle_sigma,-outer_circle_pi,'k-')" << std::endl;
            script << "    ax5.plot(inner_circle_sigma,inner_circle_pi,'k-')" << std::endl;
            script << "    ax5.plot(inner_circle_sigma,-inner_circle_pi,'k-')" << std::endl;
            script << "    ax5.plot(outer_circle_sigma,outer_circle_pi,'k-')" << std::endl;
            script << "    ax5.plot(outer_circle_sigma,-outer_circle_pi,'k-')" << std::endl;
            script << std::endl;
            script << "    fig.tight_layout()" << std::endl;
            script << "    fig.savefig('correlation_measurements.eps')" << std::endl;

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
            std::map<double,double> interpolation_map_ = lya_auto[i].interpolation_map();
            
            for (std::map<double,double>::const_iterator it = interpolation_map_.begin(); it != interpolation_map.end(); it ++){
                result << i << " " << (*it).first << " " << (*it).second;
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
            script << "# plotting lya autocorrelation" << std::endl;
            script << "fig = plt.figure(figsize=(18,9))" << std::endl;
            script << "ax = fig.add_subplot(1,1,1)" << std::endl;
            script << "ax.tick_params(axis='both',which='major',labelsize=35,length=6,width=2)" << std::endl;
            script << "ax.tick_params(axis='both',which='minor',labelsize=35,length=4,width=1)" << std::endl;
            script << "ax.set_xlabel(r'$z$',fontsize=35)" << std::endl;
            script << "ax.set_ylabel(r'$\xi$',fontsize=35)" << std::endl;
            script << "for aux_n in range(0, max_n):"
            script << "    data_z = [z for (z,n) in zip(data['z'],data['n']) if n == aux_n]"
            script << "    data_xi = [xi for (xi,n) in zip(data['lya_auto'],data['n']) if n == aux_n]"
            script << "    ax.plot(data_z ,data_xi, '-', label='n = '+srt(n))" << std::endl;
            script << "fig.savefig('" << filename << ".eps')" << std::endl;
            
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




