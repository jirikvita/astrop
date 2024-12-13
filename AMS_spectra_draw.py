#!/usr/bin/python
# Tue 9 Nov 2021

# more tutorials:
# https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/subplots_demo.html
# https://matplotlib.org/2.0.2/api/lines_api.html
# 


from math import sqrt, pow, log, exp
import os, sys, getopt
import scipy.stats as stats

import numpy as np


##########################################

from scipy.optimize import curve_fit

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
# https://stackoverflow.com/questions/19165259/python-numpy-scipy-curve-fitting
def fit_func(x, a, b):
    return a*np.power(x,b)

##########################################

import matplotlib.pyplot as plt

from collections import OrderedDict

##########################################

def ReadData(fname = 'file.txt', npow = 3):
    infile = open(fname, 'r')
    x = []
    y = []
    ex = []
    ey = [] # stat errors
    esy = [] # syst errors
    iline = -1
    for xline in infile.readlines():
        line = xline[:-1]
        iline = iline + 1
        if line[0] == '#' or iline == 0:
            continue

        tokens = line.split(',')
        R1 = float(tokens[0])
        R1,R2,N,pb,pb_stat,pb_syst,pb_over_p,pb_over_p_stat,pb_over_p_syst = float(tokens[0]),float(tokens[1]),float(tokens[2]),float(tokens[3]),float(tokens[4]),float(tokens[5]),float(tokens[6]),float(tokens[7]),float(tokens[8])
        R = 0.5 * (R2+R1)
        print((R, pb))
        x.append(R)
        # try: y.append(1.*pb_over_p)
        # was: 
        y.append(pb)
        ex.append(0.)
        ey.append(1.*pb_stat)
        esy.append(1.*pb_syst)

    infile.close()
    return x,y,ex,ey, esy

stuff = []

##########################################
# https://www.tutorialspoint.com/python/python_command_line_arguments.htm
def main(argv):
    #if len(sys.argv) > 1:
    #  foo = sys.argv[1]

    ### https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    ### https://pymotw.com/2/getopt/
    ### https://docs.python.org/3.1/library/getopt.html
    gBatch = False
    gTag=''
    print((argv[1:]))
    try:
        # options that require an argument should be followed by a colon (:).
        opts, args = getopt.getopt(argv[2:], 'hbt:', ['help','batch','tag='])

        print('Got options:')
        print(opts)
        print(args)
    except getopt.GetoptError:
        print('Parsing...')
        print ('Command line argument error!')
        print(('{:} [ -h -b --batch -tTag --tag="MyCoolTag"]]'.format(argv[0])))
        sys.exit(2)
    for opt,arg in opts:
        print(('Processing command line option {} {}'.format(opt,arg)))
        if opt == '-h':
            print(('{:} [ -h -b --batch -tTag --tag="MyCoolTag"]'.format(argv[0])))
            sys.exit()
        elif opt in ("-b", "--batch"):
            gBatch = True
        elif opt in ("-t", "--tag"):
            gTag = arg
            print(('OK, using user-defined histograms tag for output pngs {:}'.format(gTag,) ))

    print('*** Settings:')
    print(('tag={:}, batch={:}'.format(gTag, gBatch)))


    fname = 'ams_data/table-smi.csv'
    emin = 1.
    emax = 500.
   
    x,y,ex,ey, esy = ReadData(fname, 0)

    # https://stackoverflow.com/questions/773814/plot-logarithmic-axes-with-matplotlib-in-python
    # https://matplotlib.org/stable/api/markers_api.html
    
    fig, ax = plt.subplots()
    #ax.plot(x,y)
    #plt.yscale()
    #plt.yscale('log',base=2) 
    ax.semilogy(x, y, color='red', lw=2, marker='o', ls='solid')
    #ax.loglog(x, y, color='red', lw=2)
    i = 0
    ax.set(xlabel='R [GV]', ylabel='flux [..]', title='AMS antiproton flux data')
    #ax.set_yscale('log')

    # initial index for the beginning of the fit:
    i1 = int(len(x)/3)
    xx = np.array(x[i1:])
    yy = np.array(y[i1:])
    eyy = np.array(ey[i1:])
    esyy = np.array(esy[i1:])

    # up and down, stat:
    ax.semilogy(xx, yy + eyy, color='red', lw=1, ls='dashed')
    ax.semilogy(xx, yy - eyy, color='red', lw=1, ls='dashed')
    # stat+syst
    ax.semilogy(xx, yy + np.sqrt(np.power(eyy,2) + np.power(esyy, 2)), color='red', lw=1, ls='dotted')
    ax.semilogy(xx, yy - np.sqrt(np.power(eyy,2) + np.power(esyy, 2)), color='red', lw=1, ls='dotted')
    

    # fit parameters and fit covariance matrix:
    popt, pcov = curve_fit(fit_func, xx, yy, sigma = eyy)
    a, b = popt
    print((a,b, pcov))

    yfit = [ fit_func(xi, a, b) for xi in x[i1:]]
    yyfit = np.array(yfit)
    ax.semilogy(xx, yyfit, color='blue', lw=1, label='fit: a=%5.3f, b=%5.3f' % tuple(popt) )
    plt.legend()
    
    # https://www.statology.org/chi-square-goodness-of-fit-test-python/
    #chi2test = stats.chisquare(f_obs=yy, f_exp=yyfit)
    #print(chi2test)
    
    ax.grid()
    plt.savefig('ams_{}.png'.format(i))
    plt.savefig('ams_{}.pdf'.format(i))
    plt.show()


    
    
    return

###################################
###################################
###################################

if __name__ == "__main__":
    # execute only if run as a script"
    main(sys.argv)
    
###################################
###################################
###################################

