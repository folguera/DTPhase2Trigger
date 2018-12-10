#!/usr/bin/python

#
# PyROOT study of DT local reco
# using a Hough Transform
#
import os, re, sys, pickle, time, argparse, math, copy

from pprint import pprint
import ROOT as r
from array import array
from ctypes import *
from DataFormats.FWLite import Events, Handle
import numpy as np
import warnings as wr


oldargv = sys.argv[:]
sys.argv = [ '-b-' ]

r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(0)
sys.argv = oldargv

ZPOS = [0., 1.3, 2.6, 3.9, 23.8, 25.1, 26.4, 27.7]  # in cm
# FOR MB2:
wirepos_x = [-95.0142, -97.1103, -95.0094, -97.1046, 113.4, 115.5, 113.4, 115.5, -97.0756, -99.1991, -97.0935, -99.1933]
# wirepos_x = [0.0, -2.0961, 0.0048, -2.0904, 99999., 99999., 99999., 99999., -2.0614, -4.1849, -2.0793, -4.1791] # in cm (relative to wire 1)
wirepos_z = [0., 1.3, 2.6, 3.9, 99999., 99999., 99999., 99999., 23.8, 25.1, 26.4, 27.7]  # in cm

YCELLSIZE = 13.  # in mm
XCELLSIZE = 40.  # in mm, 50 cells for MB1

XMAX = 1990 # in mm, for MB1
YMAX = 290  # in mm

NYBINS = 23 # each cell own bin, middle empty
NXBINS = 8000

cellLength     = 4.2  #  The length of a cell
cellSemiLength = cellLength/2  # length from the wire to the end of the cell
cellQuarterLength = cellSemiLength/2  # length from the wire to the middle of space between the border and the wire
chamberLength  = 210.  # The length of the whole chamber

vdrift = 0.000545  #cm/ns

### other global variables:
n_events_limit = 20000
doFanae        = False


##def generateSignal(nmuons, eff):
##    ## x-y pairs
##    hits = {}
##        
##    line = r.TF1("line","pol1",0,1990) # in mm
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
        self.binsangle = 3600
        self.binsrho   = 500
        self.doVarious = False
        self.nhits     = 20


    def print_canvas(self, canvas, output_name_without_extention):
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        canvas.Print("%s/%s.png"%(self.path,  output_name_without_extention))
        canvas.Print("%s/%s.pdf"%(self.path,  output_name_without_extention))
        canvas.Print("%s/%s.root"%(self.path, output_name_without_extention))


    def getHitPosition(self, sl,l):
        index     = (int(sl)-1)*4 + (int(l)-1) # Index for a given SL and layer
        return wirepos_z[index]


    def getLRHitPosition(self, sl, l, pos, lr, wire):
        index = (int(sl)-1)*4 + (int(l)-1) # Index for a given SL and layer
        relativeChamberPosition = (pos - wirepos_x[index]) + chamberLength + 1.35;
        relativeCellPosition    = math.fmod(relativeChamberPosition, cellLength);
        positionInCell          = relativeCellPosition - cellSemiLength;
        return pos, pos + math.pow(-1, lr) * 2 * abs(positionInCell)


    def getWirePosition(self, sl, l, wire):
        index = (int(sl)-1)*4 + (int(l)-1) # Index for a given SL and layer
        x = wirepos_x[index] + wire*cellLength
        return x, wirepos_z[index]


    def getDigiPosition(self, sl,l,wire,dtime):
        posInCell = dtime*vdrift
        (xwire, zwire) = getWirePosition(sl,l,wire)
        return xwire-posInCell,xwire+posInCell,zwire


    def makeHoughTransform(self, occupancy):
        linespace = r.TH2D("linespace", "linespace", self.binsangle, 0, 2*math.pi, self.binsrho, -250., 250.)
        
        ## loop over the occupancy (left-right)
        x = r.Double(0.)
        z = r.Double(0.)
        rho = 0.
        phi = 0.
        
        for i in range(occupancy.GetN()):
            occupancy.GetPoint(i,x,z)
            
            phi = 0.
            while phi < 2*math.pi:
                rho = x*math.cos(phi)+z*math.sin(phi)
                linespace.Fill(phi,rho)
                phi = phi + 0.001
            
        return linespace


    def findLocalMaxima(self, histo2D):
        maxima = []; maximaerr = []; # maximaerr still not implemented
        spc    = r.TSpectrum2()
        
        spc.Search(histo2D, 5, "noMarkov", 0.55)
        #print histo2D.ShowPeaks(1, "noMarkov", 0.8)
        
        npks = spc.GetNPeaks()
        if npks == 0: raise RuntimeError("We have not found any peak.")
        else:         print "\n> Number of found peaks:", npks
        posx = spc.GetPositionX(); posy = spc.GetPositionY()
        posx.SetSize(npks); posy.SetSize(npks);
        lx = list(posx); ly = list(posy);
        for i in range(npks):
            xtemp = c_int(0); ytemp = c_int(0); ztemp = c_int(0);
            histo2D.GetBinXYZ(histo2D.FindBin(lx[i], ly[i], 0), xtemp, ytemp, ztemp)
            theta_err = histo2D.GetXaxis().GetBinLowEdge(xtemp.value)
            rho_err   = histo2D.GetYaxis().GetBinLowEdge(ytemp.value)
            theta_pos = histo2D.GetXaxis().GetBinCenter(xtemp.value)
            rho_pos   = histo2D.GetYaxis().GetBinCenter(ytemp.value)
            
            if (abs(theta_pos - math.pi/2) < 0.2) or (abs(theta_pos - math.pi*3/2) < 0.2): continue
            
            m = -1. / math.tan(theta_pos)
            n = rho_pos / math.sin(theta_pos)
            maxima.append((m, n))
            #maximaerr.append((m, n))
        
        print "> Number of accepted peaks:", len(maxima)
        return maxima, maximaerr


    def findBestSegment(self, linespace):
        ## this function goes back to point-space (for now only one)
        locx = c_int(0)
        locy = c_int(0)
        locz = c_int(0)
        linespace.GetMaximumBin(locx, locy, locz)

        
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


    def run(self, mb = 2, wh = 0, se = -1):
        print "> Beginning to process..."
        if (mb <= 0) or (mb > 4):  mb = -1
        if (wh < -2) or (wh > +2): wh = -3
        if (se <= 0) or (se > 14): se = -1
        txtwh = str(wh)*(wh != -3) + "ALL"*(wh == -3)
        txtmb = str(mb)*(mb != -1) + "ALL"*(mb == -1)
        txtse = str(se)*(se != -1) + "ALL"*(se == -1)
        
        print "> Processing for wheel", txtwh + ",", "sector", txtse + ",", "and chamber", txtmb + "."
        
        occupancy   = r.TGraph()  # Geometrical plot: it will be the direct (or point) space plot.
        actualocc   = r.TGraph()  # Geometrical plot: it will be the direct (or point) space plot.
        linespace   = r.TH2D()    # Dual (or line) space plot.
        linespace1D = r.TH1D()
        
        f = r.TFile.Open(self.inputfile)
        
        ## Loop over the events:
        nevents = 0
        filled  = False
        for event in f.Get("DTTree"):
            # We produce an artificial event taking segment info from various chambers and/or sectors and/or wheels until 20 hits are collected.
            if filled: break
            
            nevents += 1
            if n_events_limit and nevents >= n_events_limit: break
            if (nevents + 1)%1000 == 0: print "Event: ", nevents+1
            
            nhits = 0; nsegs = 0; actualhits = 0;
            
            for i in range(len(event.dtsegm4D_wheel)):                               # i := segment index
                if (event.dtsegm4D_wheel[i]   != wh and wh != -3): continue
                if (event.dtsegm4D_station[i] != mb and mb != -1): continue
                if (event.dtsegm4D_sector[i]  != se and se != -1): continue
                for idtsegments in range(event.dtsegm4D_phinhits[i]):
                    if (event.dtsegm4D_phi_hitsSuperLayer[i][idtsegments] == 2): continue # we demand that more than two hits per segment coincide
    ###SF                 print "dtsegm4D_phi_hitsPos %f" %(event.dtsegm4D_phi_hitsPos[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsPosCh %f" %(event.dtsegm4D_phi_hitsPosCh[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsPosErr %f" %(event.dtsegm4D_phi_hitsPosErr[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsSide %f" %(event.dtsegm4D_phi_hitsSide[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsWire %f" %(event.dtsegm4D_phi_hitsWire[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsLayer %f" %(event.dtsegm4D_phi_hitsLayer[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsSuperLayer %f" %(event.dtsegm4D_phi_hitsSuperLayer[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsTime %f" %(event.dtsegm4D_phi_hitsTime[i][idtsegments])
    ###SF                 print "dtsegm4D_phi_hitsTimeCali %f" %(event.dtsegm4D_phi_hitsTimeCali[i][idtsegments])
                    
                    z        = self.getHitPosition(event.dtsegm4D_phi_hitsSuperLayer[i][idtsegments], event.dtsegm4D_phi_hitsLayer[i][idtsegments])

                    (x1, x2) = self.getLRHitPosition(event.dtsegm4D_phi_hitsSuperLayer[i][idtsegments],event.dtsegm4D_phi_hitsLayer[i][idtsegments],
                                                     event.dtsegm4D_phi_hitsPos[i][idtsegments], event.dtsegm4D_phi_hitsSide[i][idtsegments],
                                                     event.dtsegm4D_phi_hitsWire[i][idtsegments])
                    
                    #print "lodelao:", event.dtsegm4D_phi_hitsSide[i][idtsegments]
                    #print "\nx1", x1, "x2:", x2
                    #print "resta:", abs(x1-x2), "\n"
                    
                    occupancy.SetPoint(nhits,     x1, z)
                    actualocc.SetPoint(actualhits,x1, z)
                    occupancy.SetPoint(nhits + 1, x2, z)
                    nhits      += 2
                    actualhits += 1
                    nsegs      += 1
                    
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
            if nhits > self.nhits:
                filled = True
            else:
                for hit in range(nhits):
                    occupancy.RemovePoint(hit)
                for hit in range(actualhits):
                    actualocc.RemovePoint(hit)
        
        f.Close()
        if not filled: wr.warn("> Not filled properly!", UserWarning, stacklevel=2)
        print "> Number of hits plotted:", nhits
        print "> Number of segments selected:", nsegs
        
        ## Now Draw and save TGraph
        c1 = r.TCanvas("c1", "Occupancy", 700, 700)
        c1.SetLeftMargin(0.15)
        
        occupancy.SetTitle("x-z hits")
        occupancy.SetMarkerColor(r.kBlack)
        occupancy.SetMinimum(-2)
        occupancy.SetMaximum(30)
        occupancy.SetMarkerStyle(5)
        occupancy.SetMarkerSize(1)
        occupancy.GetXaxis().SetRange(0, 250)
        occupancy.GetXaxis().SetLimits(-140, 80)
        occupancy.GetXaxis().SetTitle("x coordinate (in cm)")
        occupancy.GetYaxis().SetTitle("z coordinate (in cm)")
        occupancy.Draw("AP")
        self.print_canvas(c1, "Occupancy_MB{c}_Wh{w}_S{s}".format(c = txtmb, w = txtwh, s = txtse))
        del c1
        
        ## now make HT
        linespace = self.makeHoughTransform(occupancy)
        if self.doVarious: maxima, maximaerr = self.findLocalMaxima(linespace)
        c2 = r.TCanvas("c2", "HT", 700, 700)
        c2.SetLeftMargin(0.15)
        
        linespace.SetTitle("Linespace")
        linespace.GetXaxis().SetTitle("#theta (rad)")
        linespace.GetYaxis().SetTitle("#rho (cm)")
        linespace.Draw("COLZ")
        self.print_canvas(c2, "HoughTransform_MB{c}_Wh{w}_S{s}".format(c = txtmb, w = txtwh, s = txtse))
        del c2
        
        c1 = r.TCanvas("c1", "Occupancy", 700, 700)
        c1.SetLeftMargin(0.15)
        occupancy.SetTitle("x-z hits")
        occupancy.SetMarkerColor(r.kBlack)
        occupancy.SetMinimum(-2)
        occupancy.SetMaximum(30)
        occupancy.SetMarkerStyle(5)
        occupancy.SetMarkerSize(1)
        occupancy.GetXaxis().SetRange(0, 250)
        occupancy.GetXaxis().SetLimits(-140, 80)
        occupancy.GetXaxis().SetTitle("x coordinate (in cm)")
        occupancy.GetYaxis().SetTitle("z coordinate (in cm)")
        occupancy.Draw("AP")
        
        actualocc.SetMarkerColor(r.kAzure)
        actualocc.SetMarkerStyle(5)
        actualocc.SetMarkerSize(1)
        actualocc.Draw("P")
        ## Now try to find segments
        if not self.doVarious:
            (m,n) = self.findBestSegment(linespace)
            segm = r.TF1("seg", "[0]*x+[1]", -150, 150)
            segm.SetParameters(m,n)
            segm.SetLineColor(r.kRed)
            segm.SetLineWidth(2)
            segm.Draw("SAME")
        else:
            hmax = []
            for i in range(len(maxima)):
                tmpsegm = r.TF1("segmax_" + str(i), "[0]*x+[1]", -150, 150)
                print maxima[i][0], maxima[i][1]
                tmpsegm.SetParameters(maxima[i][0], maxima[i][1])
                tmpsegm.SetLineColor(r.kRed)
                tmpsegm.SetLineWidth(2)
                hmax.append(copy.deepcopy(tmpsegm.Clone("segmax_" + str(i))))
                del tmpsegm
                hmax[-1].Draw("SAME")
                
        self.print_canvas(c1, "Occupancy_MB{c}_Wh{w}_S{s}_fit".format(c = txtmb, w = txtwh, s = txtse))
        del c1
        return



### ========= MAIN =======================
if __name__ == "__main__":
    print "> Beginning to study DTntuples!"
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
    
    hlpr.binsangle = 500
    #hlpr.binsrho   = 500
    hlpr.binsrho   = hlpr.binsangle
    hlpr.doVarious = True
    hlpr.nhits     = 20
    
    hlpr.run()
    
    print "> Done!\n"

#  LocalWords:  dtsegm
