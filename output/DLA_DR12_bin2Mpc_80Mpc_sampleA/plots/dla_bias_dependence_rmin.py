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
    
    # bias of the dla for different bins in r_min
    dla_r = []
    dla_r_error = []
    chi2_r = []
    rmin_list = np.arange(1, 72)
    for rmin in rmin_list:
        filename = "../fit_rmin/DLA_DR12_bin2Mpc_80Mpc_sampleA_rmin{:.1f}Mpc_rmax{:.1f}Mpc.log".format(rmin, 90.0)
        #print "reading file {}".format(filename)
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
                chi2 = "{:.4f} {}".format(float(cols[6]),cols[8][0:-1])
        f.close()

        dla_r.append(dla)
        dla_r_error.append(dla_error)
        chi2_r.append(chi2)
        
        print "{rmin:d} & ${bias:.2f}\pm{error:.2f}$ & {chi2} \\\\".format(rmin=rmin, bias=dla, error=dla_error, chi2=chi2)

    for index in range(1, len(dla_r) - 1):
        if dla_r[index] == -1.0:
            upper = index + 1
            while dla_r[upper] == -1.0 and upper < len(dla_r):
                upper += 1
            lower = index - 1
            while dla_r[lower] == -1.0 and lower < len(dla_r):
                upper -= 1
            if upper >= len(dla_r) or lower < 0:
                print "Warning, interpolation not possible, omitting point"
                continue
            
            dla_r[index] = (dla_r[lower] - dla_r[upper])/(rmin_list[lower] - rmin_list[upper])* (rmin_list[index] - rmin_list[lower]) + dla_r[lower]

        if dla_r_error[index] == -1.0 or dla_r_error[index] == 0.1:
            upper = index + 1
            while (dla_r_error[upper] == -1.0 or dla_r_error[upper] == 0.1) and upper < len(dla_r):
                upper += 1
            lower = index - 1
            while (dla_r_error[upper] == -1.0 or dla_r_error[upper] == 0.1) and lower < len(dla_r):
                upper -= 1
            if upper >= len(dla_r) or lower < 0:
                print "Warning, interpolation not possible, omitting point"
                continue
            
            dla_r_error[index] = (dla_r_error[lower] - dla_r_error[upper])/(rmin_list[lower] - rmin_list[upper])* (rmin_list[index] - rmin_list[lower]) + dla_r_error[lower]

    dla_r = np.array(dla_r)
    dla_r_error = np.array(dla_r_error)
    
    
    # bias of the overall analysis with rmin 10Mpc/h
    filename = "../fit/DLA_DR12_bin2Mpc_80Mpc_sampleA_rmin{}Mpc_rmax{}Mpc.log".format(10.0, 90.0)
    #print "reading file {}".format(filename)
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

    # bias of the overall analysis with rmin 5Mpc/h
    filename = "../fit/DLA_DR12_bin2Mpc_80Mpc_sampleA_rmin{}Mpc_rmax{}Mpc.log".format(5.0, 90.0)
    #print "reading file {}".format(filename)
    f = open (filename, 'r')
    global_dla_5mpc = -1
    global_dla_error_5mpc = -1
    chi2_5mpc = -1
    aux = 0
    for line in f.readlines():
        if aux < 10:
            aux += 1
            continue
        cols = line.split()
        if len(cols) == 0:
            continue
        if cols[0] == "bias2":
            global_dla_5mpc = float(cols[2])
            global_dla_error_5mpc = float(cols[4])
            break
        if cols[0] == "Fit":
            chi2_5mpc = "{} {}".format(cols[6],cols[8][0:-1])
    f.close()


    gs = gridspec.GridSpec(1, 1)
    gs.update(hspace=0.02)
    fontsize = 28
    labelsize = 20
    
    fig = plt.figure(figsize=(9, 7))
    fig.subplots_adjust(bottom=0.18, left=0.18)
    
    ax = fig.add_subplot(gs[0])
    #ax.set_xscale("log")
    ax.set_xlim(np.min(rmin_list), np.max(rmin_list))
    ax.tick_params(axis="both", which="both", direction="in", pad=10)
    ax.tick_params(axis="both", labelsize=labelsize)
    ax.tick_params(axis="x", which="both", length=6)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    #ax.ticklabel_format(axis="x", style="plain")
    ax.set_ylabel(r"$b_{d}$", fontsize=fontsize)
    ax.set_xlabel(r"$r_{\rm min}\left({\rm h^{-1}Mpc}\right)$", fontsize=fontsize)
    ax.plot(rmin_list[np.where(dla_r != -1)], dla_r[np.where(dla_r != -1)], 'k-')
    ax.fill_between(rmin_list[np.where(dla_r != -1)], dla_r[np.where(dla_r != -1)] + dla_r_error[np.where(dla_r != -1)], dla_r[np.where(dla_r != -1)] - dla_r_error[np.where(dla_r != -1)], color='0.9')
    ax.plot(ax.get_xlim(), [global_dla]*2, 'r--', label=r"$r_{\rm min}=10{\rm h^{-1}Mpc}$")
    ax.plot(ax.get_xlim(), [global_dla + global_dla_error]*2, 'r:')
    ax.plot(ax.get_xlim(), [global_dla - global_dla_error]*2, 'r:')
    ax.plot(ax.get_xlim(), [global_dla_5mpc]*2, 'b--', label=r"$r_{\rm min}=5{\rm h^{-1}Mpc}$")
    ax.plot(ax.get_xlim(), [global_dla_5mpc + global_dla_error_5mpc]*2, 'b:')
    ax.plot(ax.get_xlim(), [global_dla_5mpc - global_dla_error_5mpc]*2, 'b:')

    ax.legend(numpoints = 1, loc = 2, prop = {'size':20}, frameon = False)
    
    fig.savefig('dla_bias_dependence_rmin.eps')
    plt.close(fig)


main()
