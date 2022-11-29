
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import scipy 
from scipy.stats.stats import pearsonr
from matplotlib.ticker import FormatStrFormatter

# %matplotlib osx
"""
rep_scatterplot.py:  make scatterplots from Calling Cards data
options:
'-f1','--file1',help='signifcant hits file 1',required=False,default="../../organized_cbf1tye7_data/FL_Plasmid_Cbf1_Tye7/sig_promoters/sig_prom_CBF1_FL_P_Rep1.gnashy.txt"
'-f2','--file2',help='significant hits file 2',required=False,default="../../organized_cbf1tye7_data/FL_Plasmid_Cbf1_Tye7/sig_promoters/sig_prom_CBF1_FL_P_Rep2.gnashy.txt"
'-xn','--xname',help='name for x axis',required=False,default = "Cbf1 Rep1 TPH"
'-yn','--yname',help='name for y axis',required=False,default = "Cbf1 Rep2 TPH"
"""


def scatterplot_cc(plot_frame,xname = "X",yname = "Y"):
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout':True})
    fig = plt.figure(figsize = (9,6))
    ax1 = fig.add_subplot(111)
    #plt.tight_layout
    #define experimental hop cutoff

    props = dict(alpha=1, edgecolors='none' )
    #find line of best fit
    xd = plot_frame[xname]
    yd = plot_frame[yname]
    par = np.polyfit(xd, yd, 1, full=True)
    slope=par[0][0]
    intercept=par[0][1]
    xl = [min(xd), max(xd)]
    yl = [slope*xx + intercept  for xx in xl]
    variance = np.var(yd)
    residuals = np.var([(slope*xx + intercept - yy)  for xx,yy in zip(xd,yd)])
    Rsqr = np.round(1-residuals/variance, decimals=2)
    R = np.sqrt(Rsqr)
    print(R)
    R_test = pearsonr(xd,yd)
    print(R_test)
    plt.text(.13,.9,'r = %0.2f'% R, fontsize=24,ha="center",va="center",transform=ax1.transAxes)
    handles = []
    
    handles.append(ax1.scatter(plot_frame[xname], plot_frame[yname], c='red', s=50, marker='o', **props))
    handles.append(ax1.plot(xl, yl, '-k'))

    #ax1.set_ylim([0,3500])
    #ax1.set_xlim([0,3000])
    ax1.grid(False)
    for item in ([ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() +ax1.get_yticklabels()):
        item.set_fontsize(18)
    for label in ax1.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    #ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax1.set_xlabel(xname,fontsize=24)
    ax1.set_ylabel(yname,fontsize=24)
    plt.show()
    
def scatterplot2_cc(plot_frame,xname = "X",yname = "Y"):
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout':True})
    fig = plt.figure(figsize = (9,6))
    ax1 = fig.add_subplot(111)
    #plt.tight_layout
    #define experimental hop cutoff

    props = dict(alpha=1, edgecolors='none' )
    #find line of best fit
    xd = plot_frame[xname]
    yd = plot_frame[yname]
    par = np.polyfit(xd, yd, 1, full=True)
    slope=par[0][0]
    intercept=par[0][1]
    xl = [min(xd), max(xd)]
    yl = [slope*xx + intercept  for xx in xl]
    variance = np.var(yd)
    residuals = np.var([(slope*xx + intercept - yy)  for xx,yy in zip(xd,yd)])
    Rsqr = np.round(1-residuals/variance, decimals=2)
    R = np.sqrt(Rsqr)
    print(R)
    R_test = pearsonr(xd,yd)
    print(R_test)
    plt.text(.13,.9,'r = %0.2f'% R, fontsize=24,ha="center",va="center",transform=ax1.transAxes)
    handles = []
    
    handles.append(ax1.scatter(plot_frame[xname], plot_frame[yname], c='red', s=50, marker='o', **props))
    handles.append(ax1.plot(xl, yl, '-k'))

    #ax1.set_ylim([0,3500])
    #ax1.set_xlim([0,3000])
    ax1.grid(False)
    for item in ([ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() +ax1.get_yticklabels()):
        item.set_fontsize(18)
    for label in ax1.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax1.set_xlabel(xname,fontsize=24)
    ax1.set_ylabel(yname,fontsize=24)
    plt.show()
    
    
def scatterplot_cc_site(plot_frame, plot_frame2, xname = "X",yname = "Y"):
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout':True})
    fig = plt.figure(figsize = (9,6))
    ax1 = fig.add_subplot(111)
    #plt.tight_layout
    #define experimental hop cutoff

    props = dict(alpha=1, edgecolors='none' )
    #find line of best fit
    xd = plot_frame[xname]
    yd = plot_frame[yname]
    par = np.polyfit(xd, yd, 1, full=True)
    slope=par[0][0]
    intercept=par[0][1]
    xl = [min(xd), max(xd)]
    yl = [slope*xx + intercept  for xx in xl]
    variance = np.var(yd)
    residuals = np.var([(slope*xx + intercept - yy)  for xx,yy in zip(xd,yd)])
    Rsqr = np.round(1-residuals/variance, decimals=2)
    R = np.sqrt(Rsqr)
    print(R)
    R_test =scipy.stats.pearsonr(xd,yd)
    print(R_test)
    plt.text(.13,.9,'r = %0.2f'% R, fontsize=24,ha="center",va="center",transform=ax1.transAxes)
    handles = []
    
    #c = ["b" if y > 5 else "r" for y in y]
    
#     plot_frame.plot.scatter(x = 'xname', 
#                    y = 'yname', 
    c = ['red' if x else 'darkmagenta' for x in (plot_frame.sites > 0)] #,
#                    ax = ax1
#                   )


    handles.append(ax1.scatter(plot_frame[xname], plot_frame[yname], c=c, s=50, marker='o', **props))
    handles.append(ax1.plot(xl, yl, '-k'))
    #handles.append(ax1.scatter(plot_frame2.Occupancy, plot_frame2.TPH, color="orange"))



    #ax1.set_ylim([0,3500])
    #ax1.set_xlim([0,3000])
    ax1.grid(False)
    for item in ([ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() +ax1.get_yticklabels()):
        item.set_fontsize(18)
    for label in ax1.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    #ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax1.set_xlabel("Observed Binding (TPH)",fontsize=24)
    ax1.set_ylabel("Predicted Occupancy",fontsize=24)
    
    

# loop through to annotate multiple datapoints
#     for i in range(plot_frame2.shape[0]):
#         plt.annotate(plot_frame2.names.tolist()[i], (df_strong_gal.Predicted_Occupancy.tolist()[i], df_strong_gal.Observed_Binding_TPH.tolist()[i]))
#     print(df_strong_gal.shape[0])
    plt.tight_layout()


    plt.show()
    