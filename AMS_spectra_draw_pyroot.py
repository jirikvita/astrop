#!/usr/bin/python
# Tue 9 Nov 2021



#
# 
#


import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

from collections import OrderedDict


def ReadData(fname = 'file.txt', npow = 3):
    infile = open(fname, 'r')
    gr = ROOT.TGraphErrors()
    ip = 0
    iline = -1
    for xline in infile.readlines():
        line = xline[:-1]
        iline = iline + 1
        if line[0] == '#' or iline == 0:
            continue
        #  \pm stat dex [km-2 yr-1 sr-1 eV-1]

        tokens = line.split(',')
        R1,R2,N,pb,pb_stat,pb_syst,pb_over_p,pb_over_p_stat,pb_over_p_syst = float(tokens[0]),float(tokens[1]),float(tokens[2]),float(tokens[3]),float(tokens[4]),float(tokens[5]),float(tokens[6]),float(tokens[7]),float(tokens[8])
        Phi = 1.*pb
        R = 0.5 * (R2+R1)
        print(ip, R, Phi)
        gr.SetPoint(ip, R, Phi)
        gr.SetPointError(ip, 0., 0.)
        ip = ip + 1

    gr.SetName('grPhi{}'.format(npow))
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(1.)
    gr.SetMarkerColor(ROOT.kBlack)
    gr.SetLineColor(ROOT.kBlack)
    infile.close()
    return gr
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

    ROOT.gStyle.SetPadLeftMargin(0.20)

    fname = 'ams_data/table-smi.csv'
    cw, ch = 1200, 800
    
    canname = 'AMS'
    can = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
    cans.append(can)
    #filename = 'foo.root'
    #rfile = ROOT.TFile(filename, 'r')
    #hname = 'histo_h'
    #h1 = rfile.Get(hname)
    #stuff.append(h1)

    emin = 1.
    emax = 500.
   
    #ROOT.gPad.Update()
    ROOT.gPad.SetLogy(1)
    ROOT.gPad.SetLogx(1)
    ROOT.gPad.SetGridy(1)
    ROOT.gPad.SetGridx(1)
    ROOT.gPad.Update()

    spect = ReadData(fname, 0)
    spect.Draw('AP')
    spect.GetXaxis().SetTitle('R [GV]')
    spect.GetYaxis().SetTitle('Flux')

    ROOT.gPad.SetLogy(1)
    ROOT.gPad.SetLogx(1)
    ROOT.gPad.SetGridy(1)
    ROOT.gPad.SetGridx(1)
    ROOT.gPad.Update()

    for can in cans:
        can.Print(can.GetName() + '.pdf')
        can.Print(can.GetName() + '.png')
    
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

