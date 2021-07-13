#!/usr/bin/python
# Tue 13 Jul 15:20:25 CEST 2021

from __future__ import print_function

#
#https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.121106
#


import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

from collections import OrderedDict

cans = []
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
    print(argv[1:])
    try:
        # options that require an argument should be followed by a colon (:).
        opts, args = getopt.getopt(argv[2:], 'hbt:', ['help','batch','tag='])

        print('Got options:')
        print(opts)
        print(args)
    except getopt.GetoptError:
        print('Parsing...')
        print ('Command line argument error!')
        print('{:} [ -h -b --batch -tTag --tag="MyCoolTag"]]'.format(argv[0]))
        sys.exit(2)
    for opt,arg in opts:
        print('Processing command line option {} {}'.format(opt,arg))
        if opt == '-h':
            print('{:} [ -h -b --batch -tTag --tag="MyCoolTag"]'.format(argv[0]))
            sys.exit()
        elif opt in ("-b", "--batch"):
            gBatch = True
        elif opt in ("-t", "--tag"):
            gTag = arg
            print('OK, using user-defined histograms tag for output pngs {:}'.format(gTag,) )

    if gBatch:
        ROOT.gROOT.SetBatch(1)

    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, gBatch))

    canname = 'can'
    can = ROOT.TCanvas(canname, canname)
    cans.append(can)
    #filename = 'foo.root'
    #rfile = ROOT.TFile(filename, 'read')
    #hname = 'histo_h'
    #h1 = rfile.Get(hname)
    #stuff.append(h1)

    emin = 3.e18
    emax = 1.e20
    fitform = 'x^3 * [0] * (x/10^18.5)^(-[3]) * (1 + (x/[2])^(1./[1]))^(([3]-[4])*[1]) * (1 + (x/[5])^(1./[1]))^(([4]-[6])*[1]) * (1 + (x/[7])^(1./[1]))^(([6]-[8])*[1])'
    #fitform = 'x^3 * [0] * (x/10^18.5)^(-[3]) * TMath::Power(1 + TMath::Power(x/[2],1./[1]),([3]-[4])*[1])'
    fitname = 'AugerFit2021'
    fit = ROOT.TF1(fitname, fitform, emin, emax)
    #fit.SetParameters(1.315e-18, 3.29)
    
    fit.SetNpx(1000)
    pars = OrderedDict()
    pars['J0'] = 1.315e-18 # [0]
    pars['omega'] = 0.05 # [1]
    pars['E12'] = 5.e18  # [2]
    pars['gamma1'] = 3.29 # [3]
    pars['gamma2'] = 2.51 # [4]
    pars['E23'] = 13.e18 # [5]
    pars['gamma3'] = 3.05 # [6]
    pars['E34'] = 46.e18 # [7]
    pars['gamma4'] = 5.1 # [8]
    
    ip = -1
    for label in pars:
        ip = ip + 1
        if ip > fit.GetNpar():
            break
        fit.SetParameter(ip, pars[label])
        fit.SetParName(ip, label)

    for ip in range(0,fit.GetNpar()):
        print('fit par{} {}: {}'.format(ip, fit.GetParName(ip), fit.GetParameter(ip)))
        
    fit.Draw()
    ROOT.gPad.SetLogy(1)
    ROOT.gPad.SetLogx(1)
    ROOT.gPad.SetGridy(1)
    ROOT.gPad.SetGridx(1)
    ROOT.gApplication.Run()
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

