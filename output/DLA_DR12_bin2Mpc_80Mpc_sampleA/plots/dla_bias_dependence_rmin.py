from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.colors import colorConverter
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter

import math

import matplotlib.ticker as ax
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, NullFormatter, ScalarFormatter

"""
    EXPLANATION:
    Plots the dla bias against z
    """


def main():
    
    rmin = 5.0
    rmax = 90.0
    rstep = 1.35#10
    
    
    #r_interval_limits = range(rmin, rmax + rstep, rstep)
    r_interval_limits = []
    aux = rmin
    while aux < rmax:
        r_interval_limits.append(aux)
        aux *= rstep
    r_interval_limits.append(aux)

    r_mid_interval = [(item + r_interval_limits[index + 1])/2.0 for index, item in enumerate(r_interval_limits[:-1])]
    r_half_interval = [(r_interval_limits[index + 1] - item)/2.0 for index, item in enumerate(r_interval_limits[:-1])]
    
    # bias of the dla for different bins in r
    dla_r = []
    dla_r_error = []
    chi2_r = []
    aux_r = rmin
    for i in range(1, 81):
        filename = "../fit_r_min/DLA_DR12_bin2Mpc_80Mpc_sampleA_rmin{}Mpc_rmax{}Mpc.log".format(aux_r, aux_r*rstep)
        print "reading file {}".format(filename)
        f = open (filename, 'r')
        dla = -1
        dla_error = -1
        chi2 = -1
        aux = 0
        for line in f.readlines():
            if aux < 10:
                aux += 1
                continue
            cols = line.split()
            if len(cols) == 0:
                continue
            if cols[0] == "bias2":
                dla = float(cols[2])
                dla_error = float(cols[4])
                break
            if cols[0] == "Fit":
                chi2 = "{} {}".format(cols[6],cols[8][0:-1])
        f.close()

        dla_r.append(dla)
            dla_r_error.append(dla_error)
        chi2_r.append(chi2)
        
        print "{rmin} & {rmax} & ${bias}\pm{error}$ & {chi2} \\\\".format(rmin=aux_r, rmax=aux_r*rstep, bias=dla, error=dla_error, chi2=chi2)
        
        aux_r *= rstep
    
    # bias of the overall analysis
    filename = "../fit/DLA_DR12_bin2Mpc_80Mpc_sampleA_rmin{}Mpc_rmax{}Mpc.log".format(5.0, 90.0)
    print "reading file {}".format(filename)
    f = open (filename, 'r')
    global_dla = -1
    global_dla_error = -1
    chi2 = -1
    aux = 0
    for line in f.readlines():
        if aux < 10:
            aux += 1
            continue
        cols = line.split()
        if len(cols) == 0:
            continue
        if cols[0] == "bias2":
            global_dla = float(cols[2])
            global_dla_error = float(cols[4])
            break
        if cols[0] == "Fit":
            chi2 = "{} {}".format(cols[6],cols[8][0:-1])
    f.close()
    
    gs = gridspec.GridSpec(1, 1)
    gs.update(hspace=0.02)
    fontsize = 28
    labelsize = 20
    
    fig = plt.figure(figsize=(9, 7))
    fig.subplots_adjust(bottom=0.18, left=0.18)
    
    ax = fig.add_subplot(gs[0])
    ax.set_xscale("log")
    ax.set_xlim(4, 110)
    ax.tick_params(axis="both", which="both", direction="in", pad=10)
    ax.tick_params(axis="both", labelsize=labelsize)
    ax.tick_params(axis="x", which="both", length=6)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    #ax.ticklabel_format(axis="x", style="plain")
    ax.set_ylabel(r"$b_{d}$", fontsize=fontsize)
    ax.set_xlabel(r"$r\left({\rm h^{-1}Mpc}\right)$", fontsize=fontsize)
    ax.errorbar(r_mid_interval, dla_r, yerr=dla_r_error, xerr=r_half_interval, fmt='k.')
    ax.plot(ax.get_xlim(), [global_dla]*2, 'k--')
    ax.plot(ax.get_xlim(), [global_dla + global_dla_error]*2, 'k:')
    ax.plot(ax.get_xlim(), [global_dla - global_dla_error]*2, 'k:')
    
    xlim = ax.get_xlim()
    ax.set_xlim(rmin-1, xlim[1])
    #ylim = ax.get_ylim()
    #if ylim[0] < 0.0:
    #    ax.set_ylim(0, ylim[1])
    
    #ax.legend(numpoints = 1, loc = 4, prop = {'size':20}, frameon = False)
    
    fig.savefig('dla_bias_dependence_rmin.eps')
    plt.close(fig)


main()
