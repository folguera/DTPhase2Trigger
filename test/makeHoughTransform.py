#!/usr/bin/python

#
# PyROOT study of DT local reco
# using a Hough Transform
#
import os, re, ROOT, sys, pickle, time
from pprint import pprint
from math import *
from array import array
from ctypes import *
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

### other global variables:
n_events_limit=10000
doFanae = False





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



class TransformHelper:
    '''Scaffold for doing the Hough transform.'''
    def __init__(self, sourcefile = ""):
        self.inputfile = sourcefile
        self.path      = ""


    def print_canvas(self, canvas, output_name_without_extention):
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        canvas.Print("%s/%s.png"%(self.path,  output_name_without_extention))
        canvas.Print("%s/%s.pdf"%(self.path,  output_name_without_extention))
        canvas.Print("%s/%s.root"%(self.path, output_name_without_extention))


    def getHitPosition(self, sl,l):
        wirepos_z = [0., 1.3, 2.6, 3.9, 99999., 99999., 99999., 99999., 23.8, 25.1, 26.4, 27.7]  # in mc
        index = (int(sl)-1)*4 + (int(l)-1) # Index for a given SL and layer
        
        return wirepos_z[index]


    def getLRHitPosition(self, sl,l,pos,lr,wire):
        wirepos_x = [-95.0142, -97.1103, -95.0094, -97.1046, 113.4, 115.5, 113.4, 115.5,-97.0756,-99.1991,-97.0935,-99.1933 ]
        
        index = (int(sl)-1)*4 + (int(l)-1) # Index for a given SL and layer
        cellLength     = 4.2  #  The length of a cell
        cellSemiLength = 2.1  # length from the wire to the end of the cell
        chamberLength  = 210.  # The length of the whole chamber

        vdrift = 0.000545  #cm/ns
        relativeChamberPosition = (pos - wirepos_x[index]) + chamberLength + 1.35;
        relativeCellPosition    = math.fmod(relativeChamberPosition, cellLength);
        positionInCell          = relativeCellPosition - cellSemiLength;
        
        return pos, pos + math.pow(-1, lr) * 2 * abs(positionInCell)


    def getWirePosition(self, sl, l, wire):
        ##Numbers for MB2
        wirepos_z = [0., 1.3, 2.6, 3.9, 99999., 99999., 99999., 99999., 23.8, 25.1, 26.4, 27.7]  # in mc
        wirepos_x = [-95.0142, -97.1103, -95.0094, -97.1046, 113.4, 115.5, 113.4, 115.5,-97.0756,-99.1991,-97.0935,-99.1933 ]
    ##    wirepos_x = [0.0, -2.0961, 0.0048, -2.0904, 99999., 99999., 99999., 99999., -2.0614, -4.1849, -2.0793, -4.1791] # in cm (relative to wire 1)
        
        index = (int(sl)-1)*4 + (int(l)-1) # Index for a given SL and layer
        cellLength     = 4.2  #  The length of a cell
        cellSemiLength = 2.1  # length from the wire to the end of the cell
    #    chamberLength  = 210  # The length of the whole chamber
        
        x = wirepos_x[index] + wire*cellLength
        return x, wirepos_z[index]


    def getDigiPosition(self, sl,l,wire,dtime):
        vdrift = 0.000545  #cm/ns
        posInCell = dtime*vdrift
        (xwire, zwire) = getWirePosition(sl,l,wire)

        return xwire-posInCell,xwire+posInCell,zwire


    def makeHoughTransform(self, occupancy):
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


    def findBestSegment(self, linespace):
        ## this function goes back to point-space (for now only one)
        locx = c_int(0)
        locy = c_int(0)
        locz = c_int(0)
        maxbin = linespace.GetMaximumBin(locx,locy,locz)

        
        theta_err = linespace.GetXaxis().GetBinLowEdge(locx.value)
        rho_err   = linespace.GetYaxis().GetBinLowEdge(locy.value)
        theta_pos = linespace.GetXaxis().GetBinCenter(locx.value)
        rho_pos   = linespace.GetYaxis().GetBinCenter(locy.value)
        
        m = -1./math.tan(theta_pos)
        n = rho_pos/math.sin(theta_pos)
        print "\n> Maximum bin in the dual space. x:", locx, "y:", locy, "z:", locz
        print "> Maximum bin in the dual space in terms of rho and theta. rho:", rho_pos, "theta:", theta_pos
        print "> Associated line in the direct space. m:", m, "n:", n, "\n"
        return m,n


    def run(self):
        print "> Beginning to process..."
        
        occupancy   = ROOT.TGraph()
        linespace   = ROOT.TH2D()
        linespace1D = ROOT.TH1D()
        
        f = ROOT.TFile.Open(self.inputfile)
        
        ## Loop over the events:
        nevents=0
        filled=False
        for event in f.Get("DTTree"):
            nevents += 1
            if n_events_limit and nevents>=n_events_limit: break
            if (nevents+1)%10000==0: print "Event: ",nevents+1
            
            nhits = 0
            if filled: break 
            for i in range(len(event.dtsegm4D_wheel)):
                if (event.dtsegm4D_wheel[i]!=0): continue
                if (event.dtsegm4D_station[i]!=2): continue
                for idtsegments in range(event.dtsegm4D_phinhits[i]):
                    if (event.dtsegm4D_phi_hitsSuperLayer[i][idtsegments]==2): continue
    ###SF                 print "dtsegm4D_phi_hitsPos %f" %(event.dtsegm4D_phi_hitsPos[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsPosCh %f" %(event.dtsegm4D_phi_hitsPosCh[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsPosErr %f" %(event.dtsegm4D_phi_hitsPosErr[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsSide %f" %(event.dtsegm4D_phi_hitsSide[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsWire %f" %(event.dtsegm4D_phi_hitsWire[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsLayer %f" %(event.dtsegm4D_phi_hitsLayer[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsSuperLayer %f" %(event.dtsegm4D_phi_hitsSuperLayer[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsTime %f" %(event.dtsegm4D_phi_hitsTime[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsTimeCali %f" %(event.dtsegm4D_phi_hitsTimeCali[i][idtsegments])
                    
                    z = self.getHitPosition(event.dtsegm4D_phi_hitsSuperLayer[i][idtsegments],event.dtsegm4D_phi_hitsLayer[i][idtsegments])

                    (x1,x2) = self.getLRHitPosition(event.dtsegm4D_phi_hitsSuperLayer[i][idtsegments],event.dtsegm4D_phi_hitsLayer[i][idtsegments],
                                                    event.dtsegm4D_phi_hitsPos[i][idtsegments],event.dtsegm4D_phi_hitsSide[i][idtsegments],
                                                    event.dtsegm4D_phi_hitsWire[i][idtsegments])
                    
                    occupancy.SetPoint(nhits,x1,z)
                    occupancy.SetPoint(nhits+1,x2,z)
                    nhits += 2

    ##        for i in range(len(event.digi_wheel)):
    ##            if (event.digi_wheel[i]!=0): continue
    ##            if (event.digi_station[i]!=2): continue
    ##            if (event.digi_sl[i]==2):  continue
    ##            
    ##            (x1,x2,z) = getDigiPosition(event.digi_sl[i],event.digi_layer[i],event.digi_wire[i],event.digi_time[i])
    ##            
    ##            print "digi pos: %f or %f"    %(x1,x2)
    ##            print "digi z: %f"      %(z)
    ##            print "digi time: %f"   %(event.digi_time[i])

            ## if enough
            if nhits > 20:
                filled = True
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

        self.print_canvas(c1, "Occupancy_MB2_Wh0")
        
        ## now make HT
        linespace = self.makeHoughTransform(occupancy)
        c2 = ROOT.TCanvas("c2", "HT",700,700)
        c2.SetLeftMargin(0.15)
        
        linespace.SetTitle("Linespace")
        linespace.GetXaxis().SetTitle("#theta (rad)")
        linespace.GetYaxis().SetTitle("#rho (cm)")
        linespace.Draw("COLZ")
        self.print_canvas(c2,"HoughTransform_MB2_Wh0")

        ## Now try to find segments
        (m,n) = self.findBestSegment(linespace)
        
        segm = ROOT.TF1("seg","[0]*x+[1]",-150,150)
        segm.SetParameters(m,n)
        c1.cd()
        segm.SetLineColor(ROOT.kRed)
        segm.SetLineWidth(2)
        segm.Draw("SAME")
        self.print_canvas(c1,"Occupancy_MB2_Wh0_fit")
        f.Close()
        return



### ========= MAIN =======================
if __name__ == "__main__":
    print "> Beginning to study DTntuples!"
    import argparse
    parser = argparse.ArgumentParser(usage = "MakeHoughTransform.py [options]",
                                     description="Dummy test on Hough transform",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-i", "--input", dest = "inputfile", default = "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2018/Prompt/DTTree_SingleMuon_ZMuSkim315252-315510.root", help = 'Input filename', required = False)
    parser.add_argument("-u", "--uvieu", dest = "inFanae",   default = False,    action= "store_true", help = 'Drop the result plots in fanae', required = False)
    #parser.add_argument("--sim", action='store_true')
    #parser.add_argument("-o", "--output",dest = "output",    default = "plots/", help = 'Folder name to store results', required = False)

    args = parser.parse_args()
    inputfile = args.inputfile
    inFanae   = args.inFanae
    print "We will analyse", inputfile
    
    hlpr = TransformHelper(inputfile)
    if not inFanae: hlpr.path = "/afs/cern.ch/user/f/folguera/www/private/L1TPhase2/DTwithHT_withRecHits1D/"
    else:           hlpr.path = "/nfs/fanae/user/vrbouza/www/Proyectos/trigger_primitives/results/houghtrans/"
    hlpr.run()
    
    print "> Done!"

#  LocalWords:  dtsegm
