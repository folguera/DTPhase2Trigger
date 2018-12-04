#!/usr/bin/python

#
# PyROOT study of DT local reco 
# using a Hough Transform
#
import os, re, ROOT, sys, pickle, time
from pprint import pprint
from math import *
from array import array
from DataFormats.FWLite import Events, Handle
import numpy as np


# import ROOT in batch mode
import sys
import argparse
import math

oldargv = sys.argv[:]
sys.argv = [ '-b-' ]

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
sys.argv = oldargv

ZPOS = [0., 1.3, 2.6, 3.9, 23.8, 25.1, 26.4, 27.7]  # in cm

YCELLSIZE = 13.  # in mm
XCELLSIZE = 40.  # in mm, 50 cells for MB1

XMAX = 1990 # in mm, for MB1
YMAX = 290  # in mm

NYBINS = 23 # each cell own bin, middle empty
NXBINS = 8000

#### other global variables: 
n_events_limit=10000

##def generateSignal(nmuons, eff):
##    ## x-y pairs
##    hits = {}
##        
##    line = ROOT.TF1("line","pol1",0,1990) # in mm
##    ## restrict phi bending to +/- pi/8
##
###    for y in YPOS:
###        if (>eff) continue:
##
##    ## now duplicate hit position. 
##    return hits 
##
##def generateBackground(bkglevel):
##    hits = {}
##
##    ## bkglevel indicates the occupancy of the noise, should go from 0-1
##    if bkglevel > 1:
##        bkglevel=1
##
##        ## evenly distribute hits from noise... of course this is just partial simulation
##    
##def applyRECOeff(hits):
##    
##def findAndDrawReco(histo):
def print_canvas(canvas, output_name_without_extention, path):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png"%(path,output_name_without_extention))
    canvas.Print("%s/%s.pdf"%(path,output_name_without_extention))
    canvas.Print("%s/%s.root"%(path,output_name_without_extention))

def getWirePosition(sl, l, wire): 
    ##Numbers for MB2
    wirepos_z = [0., 1.3, 2.6, 3.9, 99999., 99999., 99999., 99999., 23.8, 25.1, 26.4, 27.7]  # in mc
    wirepos_x = [0.0, -2.0961, 0.0048, -2.0904, 99999., 99999., 99999., 99999., -2.0614, -4.1849, -2.0793, -4.1791] # in cm (relative to wire 1)


    index = (sl-1)*4 + (l-1) # Index for a given SL and layer
    cellLength     = 4.2  #  The length of a cell
    cellSemiLength = 2.1  # length from the wire to the end of the cell
#    chamberLength  = 210  # The length of the whole chamber

    x = wirepos_x[index] + wire*cellLength
    return x, wirepos_z[index]

def getDigiPosition(sl,l,wire,dtime):
    vdrift = 0.000545  #cm/ns
    
    posInCell = dtime*vdrift
    (xwire, zwire) = getWirePosition(sl,l,wire)

    return xwire-posInCell,xwire+posInCell,zwire

def makeHoughTransform(occupancy):
    
    linespace = ROOT.TH2D("linespace","linespace", 3600, 0, 2*math.pi, 500, -250., 250.) 
    
    ## loop over the occupancy (left-right) 
    x = ROOT.Double(0.)
    z = ROOT.Double(0.)
    rho = 0.
    phi = 0.
    
    for i in range(occupancy.GetN()):
        occupancy.GetPoint(i,x,z)
        
        phi = 0.
        while phi<2*math.pi: 
            rho = x*math.cos(phi)+z*math.sin(phi)
            linespace.Fill(phi,rho)
            phi = phi + 0.001 
        
    return linespace
        
    
def run(inputfile):
    print "Processing " 
    
    occupancy = ROOT.TGraph()

    linespace = ROOT.TH2D()
    linespace1D = ROOT.TH1D()
    
    f = ROOT.TFile.Open(inputfile)
    
    ## Loop over the events: 
    nevents=0
    filled=False
    for event in f.Get("DTTree"):
        nevents += 1
        if n_events_limit and nevents>=n_events_limit: break
        if (nevents+1)%10000==0: print "Event: ",nevents+1
        
        nhits = 0
        if filled: break 
        for i in range(len(event.digi_wheel)):
            if (event.digi_wheel[i]!=0): continue
            if (event.digi_station[i]!=2): continue
            if (event.digi_sl[i]==2):  continue
            
            (x1,x2,z) = getDigiPosition(event.digi_sl[i],event.digi_layer[i],event.digi_wire[i],event.digi_time[i])                                  
            
            occupancy.SetPoint(nhits,x1,z) 
            occupancy.SetPoint(nhits+1,x2,z)
            nhits += 2
        
        ## if enough
        if nhits > 20: 
            filled=True
        else: 
            for hit in range(nhits):
                occupancy.RemovePoint(hit)
    

    ## Now Draw and save TGraph
    c1 = ROOT.TCanvas("c1", "Occupancy",700,700)
    c1.SetLeftMargin(0.15)
    
    occupancy.SetTitle("x-z hits")
    occupancy.SetMarkerColor(ROOT.kBlack)
    occupancy.SetMinimum(-2)
    occupancy.SetMaximum(30)
    occupancy.SetMarkerStyle(5)
    occupancy.SetMarkerSize(1)
    occupancy.GetXaxis().SetRange(0, 250)
    occupancy.GetXaxis().SetTitle("x coordinate (in cm)")
    occupancy.GetYaxis().SetTitle("z coordinate (in cm)")
    occupancy.Draw("AP")

    
    path="/afs/cern.ch/user/f/folguera/www/private/L1TPhase2/DTwithHT/"
    print_canvas(c1,"Occupancy_MB2_Wh0",path)
    
    ## now make HT 
    linespace = makeHoughTransform(occupancy)            
    c2 = ROOT.TCanvas("c2", "HT",700,700)
    c2.SetLeftMargin(0.15)
    
    linespace.SetTitle("Linespace")
    linespace.GetXaxis().SetTitle("phi (rad)")
    linespace.GetYaxis().SetTitle("rho (cm)")
    linespace.Draw("COLZ")
    print_canvas(c2,"HoughTransform_MB2_Wh0",path)

    ## Now try to find segments
    findSegments(linespace)    
    


#### ========= MAIN =======================
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(usage="MakeHoughTransform.py [options]",
                                     description="Dummy test on Hough transform",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    #parser.add_argument("--sim", action='store_true')
    parser.add_argument("--input", dest="inputfile", default="", help='Input filename', required=True)
    #parser.add_argument("-o","--ofolder",dest="output", default="plots/", help='folder name to store results')

    args = parser.parse_args() 
    
    inputfile = args.inputfile
    print inputfile
    
    run(inputfile)
    
    print "DONE"
