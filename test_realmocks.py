#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################
# PROGRAM NAME
#
#
##############################
"""
Created on Mon Apr 21 08:24:25 2014

@author: clackner

program details
"""
import sys, os, re
import argparse
import numpy as np
from collections import defaultdict

import matplotlib.figure as figure
from matplotlib.backends.backend_ps import FigureCanvasPS as FigCanvasPS
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvasA

from scipy.stats import spearmanr, percentileofscore
import pyfits
import line_fits

import astropy.table
from astropy.cosmology import FlatLambdaCDM
import configParams as config

def readMergers(inputfolder, origfile, clean=False):
    """
    read in a set of mock mergers, including where the merger was finally detected
    """
    with open(inputfolder+'cleaned_peaks_sources.txt') as pSrc, open(inputfolder+'cleaned_peaks.txt') as pXY:
        linesSrc = pSrc.readlines()
        linesPxy = pXY.readlines()
        peaks = {}
        for i in range(len(linesSrc)):
            fields = linesSrc[i].split()
            galid = int(fields[0])
            npeaks = int(fields[5])
            peaks[galid] = defaultdict(list)
            for p in range(npeaks):
                x0 = linesPxy[i].split()[2+2*p]
                y0 = linesPxy[i].split()[2+2*p+1]
                flux = fields[6+npeaks+p]
                peaks.append((x0,y0,flux))
    
    mergers = list()
    origGals = astropy.table.Table().read(origfile)
    cosmos = FlatLambdaCDM(H0=config.H0, Om0=config.Om0)
    ang_dist_kpc = cosmos.angular_diameter_distance(origGals['Z1']).to('kpc').value
    for indg, g in enumerate(origGals):
        merger_kwargs = dict([(c, g[c]) for c in origGals.colnames])
        merger_kwargs['flux'] = (10.0**(-0.4*(g['mag1']-config.zeropt)) + 
                                10.0**(-0.4*(g['mag2']-config.zeropt)))
        merger_kwargs['sepas'] = merger_kwargs['sep_kpc']/ang_dist_kpc[indg] * 206265.
        merger_kwargs['ang_dist_kpc'] = ang_dist_kpc[indg]
        mergers.append(Merger(**merger_kwargs))
        
        if g['ID'] in peaks:
            for p in peaks[g['ID']]:
                mergers[-1].addPeak(**zip(['x0', 'y0', 'flux'], p))
        if clean:
            #this will remove peaks, but right now, I'm assuming that's already done
            #mergers[-1].cleanPeaks()
            pass
        
        mergers[-1].assignRealPeaks()
                                
    return mergers


def makesubplots(fig, nx=3, ny=3, **kwargs):
    for window in range(nx*ny):
        yield fig.add_subplot(nx,ny,window+1, **kwargs)
 
    
class Peak:
    
    def __init__(self, **kwargs):
        for key, val in kwargs.iteritems():
            setattr(self, key, val)
        
    def distTo(self, x, y, scaling=1.0):
        return np.sqrt((x-self.x0)**2 + (y-self.y0)**2) *scaling
        
    def __cmp__(self, other):
        return int(np.sign(self.distTo(config.cutoutX*0.5, config.cutoutY*0.5)-
                       other.distTo(config.cutoutX*0.5, config.cutoutY*0.5)))
        
    def __repr__(self):
        return '(%.2f, %.2f)-%.3e'%(self.x0, self.y0, self.flux)


class Merger:
    
    def __init__(self, **kwargs):
        for key, val in kwargs.iteritems():
            setattr(self, key, val)
        self.peaklist = list()
        self.mag = -2.5*np.log10(self.flux)+25.959
        
    def setReal(self):
        self.isdbl = True
    
    def toKpc(self, arcsec):
        return arcsec / 206265.0 * self.angdist_kpc
        
    def toArcsec(self, kpc):
        return kpc / (self.angdist_kpc) * 206265.
        
    def addPeak(self, **kwargs):
        self.peaklist.append(Peak(**kwargs))
        self.peaklist[-1].isGal1 = False
        self.peaklist[-1].isGal2 = False
        self.peaklist.sort()
        
    def cleanPeaks(self, ecut=0.0, centcut=100.0, totfluxcut=0.0, fluxRcut=0.0,
                   distcut=(0.0,100.0)):
                       
        if len(self.peaklist) > 0:
            self.peaklist = [p for p in self.peaklist if 
            ((self.toKpc(p.distTo(config.cutoutX/2., config.cutoutY/2., 
                                  scaling=config.pixelscale)) < centcut) &
            (p.ellip > ecut) & (p.flux/self.flux >= totfluxcut))]
           
        if len(self.peaklist) > 0:
            maxI = np.argmax([p.flux for p in self.peaklist])
            maxFlux = self.peaklist[maxI].flux
            xM, yM = self.peaklist[maxI].x0, self.peaklist[maxI].y0

            self.peaklist = [p for p in self.peaklist if 
            (p.flux==maxFlux)|((p.flux/maxFlux > fluxRcut) &
            (self.toKpc(p.distTo(xM, yM, scaling=config.pixelscale)) < distcut[1]) & 
            (self.toKpc(p.distTo(xM, yM, scaling=config.pixelscale)) >= distcut[0]))]
        
    def assignRealPeaks(self, thresh=8):
        
        self.detect1 = False
        self.detect2 = False
        self.MeasP1 = None
        self.MeasP2 = None
        self.isdbl = False
        
        if len(self.peaklist) > 0:
            self.peaklist.sort()#key = lambda p1: p1.distTo(self.x01, self.y01))
            if ((np.abs(self.peaklist[0].x0 - self.x01) < thresh) and
                (np.abs(self.peaklist[0].y0 - self.y01) < thresh)):
                    self.peaklist[0].isGal1 = True
                    self.detect1 = True     
                    self.MeasP1 = self.peaklist[0]
                    
            self.peaklist.sort()#key = lambda p1: p1.distTo(self.x02, self.y02))
            cc = 0
            while (cc < len(self.peaklist)) and (self.peaklist[cc].isGal1):
                cc += 1
            
            if ((cc < len(self.peaklist)) and
                (np.abs(self.peaklist[cc].x0 - self.x02) < thresh) and
                (np.abs(self.peaklist[cc].y0 - self.y02) < thresh)):
                self.peaklist[cc].isGal2 = True
                self.detect2 = True
                self.MeasP2 = self.peaklist[cc]
        
        self.extrapeaks = [p for p in self.peaklist 
                                       if not (p.isGal1 or p.isGal2)]
                
        #assert len(self.extrapeaks)+(1 if self.MeasP1 is not None else 0)+\
        #(1 if self.MeasP2 is not None else 0) ==len(self.peaklist)
            
        self.isdbl = (self.detect1) & (self.detect2)
        

    def getMeasPeak1(self):
        return self.MeasP1
            
    def getMeasPeak2(self):
        return self.MeasP2
            
    def getSeps(self, x0, y0):
        return [self.toKpc(pp.distTo(x0, y0, scaling=0.03))
                for pp in self.peaklist]
                
    def getSepsMaxPeak(self, safe=False):
        self.peaklist.sort(cmp=lambda x, y: int(np.sign(x.flux-y.flux)))
        if len(self.peaklist) > 1:
            maxP = np.argmax([p.flux for p in self.peaklist])
            seps = self.getSeps(self.peaklist[maxP].x0,
                                self.peaklist[maxP].y0)
            return [s for i, s in enumerate(seps) if i != maxP]
        else:
            if not safe:
                return []
            else:
                return [-1]
            
    def getFluxRatioMaxPeak(self, withMax=False, extraonly=False, safe=False):
        if extraonly:
            myplist = self.extrapeaks
        else:
            myplist = self.peaklist
        myplist.sort(cmp=lambda x, y: int(np.sign(x.flux-y.flux)))
        if len(self.peaklist) > 1:
            maxP = np.argmax([p.flux for p in self.peaklist])
            fr = [p.flux/self.peaklist[maxP].flux for p in myplist]
            if not withMax:
                return [s for i, s in enumerate(fr) if i != maxP]
        else:
            if safe:
                return [-1]
            else:
                return []    
    
    def plot(self, ax, path, peaks=True, vmax=None, title=False, realpeaks=False, 
             bar=False, bw=False, **kwargs):
        
        imgfile = path + '%06d_%06d_%04.1f.fits'%(self.id1, 
                                                  self.id2, self.sepkpc)
        data = pyfits.open(imgfile)[0].data
        xoff=70
        yoff=70
        data = data[xoff:268-xoff,yoff:268-yoff]
        if vmax is None:
            vmax = np.percentile(data.ravel(), 99.95)
        im = ax.imshow(np.arcsinh(data), 
                       vmax=vmax, **kwargs)
    
        if peaks:
            ax.plot([p.y0-xoff for p in self.peaklist],
                    [p.x0-yoff for p in self.peaklist], 'cx' if not bw else 'wx',
                    scalex=False, scaley=False, mew=1.2)
        if realpeaks:
            ax.plot([self.y01-xoff, self.y02-xoff], [self.x01-yoff, self.x02-yoff], 'mo', mfc='none',
                    mec='m' if not bw else 'w', mew=1, ms=10,
                    scalex=False, scaley=False)
        
        if title:
            mergertype={1:'early', 2:'late'}
            ax.text(0.05,0.75,"%s-%s merger\nz=%.2f"%(mergertype[self.zest1], mergertype[self.zest2],
                                                             self.z),
                     size=9, transform=ax.transAxes, ha='left', va='bottom')
                     
        if bar:
            ax.set_autoscalex_on(False)
            ax.set_autoscaley_on(False)
            ax.text(0.5,0.04, '8 kpc', horizontalalignment='center',
                    transform=ax.transAxes,
                    color='k', fontsize=9)
            size=ax.get_xlim()[-1]
            scale = self.toArcsec(8.0) / 0.03
            height = ax.transData.inverted().transform(ax.transAxes.transform((0,0.14)))[1]
            ax.plot([size/2 - scale/2.0, size/2+scale/2.0], [height, height],
                    'k-', lw=2)
                    
        ax.tick_params(labelleft='off', labelbottom='off')
        print 'plotted ', self.id1, self.id2, self.sepkpc, self.z, self.fluxratio,
        print self.mag, self.zest1, self.zest2, self.isdbl, self.isdetdbl, self.mag,
        print self.measFlux12(), self.measSep12()
        return im
        
    def measSep12(self):
        p1 = self.getMeasPeak1()
        p2 = self.getMeasPeak2()
        
        if (p1 is not None) and (p2 is not None):
            return self.toKpc(p1.distTo(p2.x0, p2.y0, scaling=0.03))
        else:
            return -1
        
    def measFlux12(self):
        p1 = self.getMeasPeak1()
        p2 = self.getMeasPeak2()
        
        if (p1 is not None) and (p2 is not None):
            return np.amin([p1.flux/p2.flux, p2.flux/p1.flux])
        else:
            return -1
            
    def flux12(self):
        try:
            f1 = self.getMeasPeak1().flux
        except AttributeError:
            f1 = -1
        try:
            f2 = self.getMeasPeak2().flux
        except AttributeError:
            f2 = -1
        return [f1, f2]



def plotPeaks(mergers, path, impath, ending, FigCanvas, ngal=None):

    grpmerge = {}
    for m in mergers:
        key = (m.id1, m.id2)
        if key not in grpmerge.keys():
            grpmerge[key] = []
        grpmerge[key].append(m)

    for key in grpmerge.keys()[:len(grpmerge) if ngal is None else ngal]:
        fig = figure.Figure(frameon=False)
        canv = FigCanvas(fig)
        subs = makesubplots(fig)
        for mm in grpmerge[key]:
            ax = subs.next()
            myim = mm.plot(ax, path, origin='lower',
                           cmap='gray', interpolation='none', aspect='equal')
    
        fig.tight_layout(pad=0.1)
        fig.savefig(impath+'peaksall_%06d_%06d'%key+ending)


def cumComplete(value, selected, bins=10, cum=False):
    """
    return the cumulative completeness (assuming lower bounds)
    """
    
    num, bins1 = np.histogram(value[selected], bins)
    denom, bins1 = np.histogram(value, bins1)
    
    frac = num[::-1]*1.0/denom[::-1]        
    cumfrac = np.cumsum(num[::-1])*1.0/np.cumsum(denom[::-1])
    
    if cum:
        return bins1[:-1]+np.diff(bins1)*0.5, cumfrac[::-1]
    else:
        return bins1[:-1]+0.5*np.diff(bins1), frac[::-1]
        
        
def cumContaim(realvalue, measvalue, bins=10, cum=True):
    """
    returns cumulative contaminations
    """
    #realval, bins1 = np.histogram(realvalue, bins=bins)
    #ratio = np.cumsum(realval)*1.0/(np.cumsum(realval[::-1])[::-1])
    #return bins1[1:,], ratio
    measval, bins1 = np.histogram(measvalue, bins=bins)
    contam = np.array([len(realvalue[(measvalue>bb)&(realvalue<bb)]) for bb in bins1[:-1]])
    frac = contam*1.0/(np.cumsum(measval[::-1])[::-1])
    return bins1[:-1], frac


def plotcomplete(ax, prop, measprop, restrict, isdbl, risdbl, nbins=15,
                 name='property', split=True, bw=False, log=False, points=False):
    """
    plots the completeness and the error in the measured parameter (if split==True)
    """
    denom, bins = np.histogram(prop, bins=nbins)
    d2, bins2 = np.histogram(prop[restrict], bins=nbins)
    d2 *= 1.0
    denom *= 1.0
    d4 = np.histogram(prop[isdbl], bins=bins)[0]*1.0
    dblrest = np.histogram(prop[risdbl&restrict], bins=bins2)[0]*1.0
    err1=np.sqrt(d4/denom**2 + d4**2/(denom**3))
    err2 = np.sqrt(dblrest/d2**2 + dblrest**2/d2**3)
    color1='b' if not bw else 'k'
    color2='c' if not bw else 'darkgray'
    ax.plot(bins[:-1], d4/denom, color=color1, linestyle='solid',
                #yerr=err1, color='b', ls='solid', marker=None,
                label=r'$\mathrm{total sample}$', marker='o' if points else 'None')
    ax.plot(bins2[:-1], dblrest/d2, color=color2, linestyle='solid', 
            lw=2., label=r'$\mathrm{restricted\ sample}$', marker='o' if points else 'None')
    if points:
        ax.errorbar(bins[:-1], d4/denom, yerr=np.sqrt(d4+d4**2/denom)/denom,
                    color=color1, ls='none', marker='None')
        ax.errorbar(bins2[:-1], dblrest/d2, yerr=np.sqrt(dblrest+dblrest**2/d2)/d2,
                    color=color2, ls='none', marker='None', )
    ax.set_ylabel(r'$\mathrm{completeness}$', size=10)
    ax.tick_params(labelsize=9)
    if log:
        ax.set_xscale('log')
    
    if split:
        ax.tick_params(labelbottom='off')
        div = make_axes_locatable(ax)
        axlower = div.append_axes("bottom",1.2, sharex=ax)
        axlower.set_xlabel(name, size=10)
        axlower.plot(prop[isdbl], (measprop[isdbl]/prop[isdbl]-1), marker='.',
                     ls='none', color=color1,
                     ms=3)
        axlower.plot(prop[risdbl&restrict], 
                     (measprop[risdbl&restrict]/prop[risdbl&restrict]-1), 
                    marker='.', color=color2, ls='none',
                     ms=3)
        axlower.axhline(0, color='k')
        axlower.tick_params(labelsize=9)
        axlower.yaxis.set_major_locator(MaxNLocator(prune='upper', nbins=5))
        axlower.set_ylabel(r'$\mathrm{fractional\ error}$', size=10)
        return axlower
    ax.set_xlabel(name, size=10)


def main():

    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("infolder", help='folder with peak finding output')
    parser.add_argument('origfile', help='file with real merger info')
    parser.add_argument('-i', '--impath', dest='impath', default='./',
                        help='path to image output')
    parser.add_argument("-e", "--eps",
                 action='store_true', default=False, help='make eps plots',
                 dest='epsPlot')
    args = parser.parse_args()
    
    
    imCent = (268./2.0, 268.0/2.0)
    mergers = list()
    tfc = 0.01
    desired_fr = 0.25
    cut_fr = 0.00


    mergers = np.asarray(readMergers(args.infolder, args.origfile))
    for m in mergers:
        m.cleanPeaks(centcut=10., totfluxcut=tfc)
        m.assignRealPeaks()
        #m.isdbl = m.isdetdbl & m.isdbl
        
            
    FigCanvas = FigCanvasPS if args.epsPlot else FigCanvasA
    ending='.eps' if args.epsPlot else '.png'

    mergers = np.asarray(mergers)
    if False:
        plotPeaks(mergers, os.path.dirname(args.origfile)+'/imgs/', 
                  args.impath, ending, FigCanvas, ngal=20)
        ee = np.asarray([(m.zest1==1)&(m.zest2==1) for m in mergers])
        for m in mergers[5:200:8]:
            print len(m.peaklist), m.getMeasPeak1(), 
            print m.getMeasPeak2(), m.x01, m.y01, m.x02, m.y02
        plotPeaks(mergers[ee],
                  os.path.dirname(args.origfile)+'/imgs/', 
                args.impath+'zest11_', ending, FigCanvas, ngal=20)


    redshift = np.array([m.z for m in mergers])
    sepkpc = np.array([m.sepkpc for m in mergers])
    fluxratio = np.array([m.fluxratio for m in mergers])
    zest1 = np.array([m.zest1 for m in mergers])
    zest2 = np.array([m.zest2 for m in mergers])
    mass = np.array([m.mass for m in mergers])
    mag = np.array([-2.5*np.log10(m.flux)+25.959 for m in mergers])
    measfluxratio = np.array([m.measFlux12() for m in mergers])
    isdbl = np.array([m.isdbl for m in mergers])
    zest11 = (zest1==1)&(zest2==1)
    zest12 = ((zest1==1)|(zest2==1))&(zest1<3)&(zest2<3)
    zest22 = (zest1==2)&(zest2==2)
    zestcut = np.in1d(zest1, [1,2]) & np.in1d(zest2, [1,2])#(zest1==zest1)&(zest2==zest2)#(zest1!=3)&(zest2!=3)
    masscut = mass>10.6
    magcut = mag < 20
    restrict = ((sepkpc > 2.2)&zestcut&(masscut)&(sepkpc<8.0)&(fluxratio>desired_fr))
    restrictM = zestcut&(masscut)
    print len(mergers[restrict])
    for m in mergers[restrict & (zest1==1) & (zest2==1) & (sepkpc==5.0)]:
        print m.id1, m.id2, m.isdbl, m.sepkpc, m.fluxratio, m.detect1, m.detect2

    fig = figure.Figure((6,6))
    canv = FigCanvas(fig)
    subs = makesubplots(fig, nx=1, ny=1)
    ax = subs.next()
    color = np.array(['r' if (m.zest1==1)&(m.zest2==1) else 
                    ('b' if (m.zest1==2)&(m.zest2==2) else 'g')
                    for m in mergers])
    ax.scatter(np.array([min(m.mass1/m.mass2, m.mass2/m.mass1) for m in mergers])[restrictM],
               fluxratio[restrictM], c=color[restrict], alpha=0.1, lw=0)
    print spearmanr(np.array([min(m.mass1/m.mass2, m.mass2/m.mass1) for m in mergers])[restrictM],
               fluxratio[restrictM])
    print spearmanr(np.array([min(m.mass1/m.mass2, m.mass2/m.mass1) for m in mergers])[restrictM&zest11],
               fluxratio[restrictM&zest11])
    print spearmanr(np.array([min(m.mass1/m.mass2, m.mass2/m.mass1) for m in mergers])[restrictM&zest12],
               fluxratio[restrictM&zest12])
    print spearmanr(np.array([min(m.mass1/m.mass2, m.mass2/m.mass1) for m in mergers])[restrictM&zest22],
               fluxratio[restrictM&zest22])
    print line_fits.robustPolyFit(np.array([min(m.mass1/m.mass2, m.mass2/m.mass1) for m in mergers])[restrictM],
               fluxratio[restrictM], 1)
    ax.set_xlabel('mass ratio')
    ax.set_ylabel('flux ratio')
    ax2 = ax.twinx()
    ax2.hist(np.array([min(m.mass1/m.mass2, m.mass2/m.mass1) for m in mergers])[restrictM],
             bins=7, histtype='step', color='k')
    fig.savefig(args.impath+'flux_mass_ratio'+ending)

    fig = figure.Figure((12,12))
    canv = FigCanvas(fig)
    subs = makesubplots(fig, nx=2)
    ax = subs.next()
    plotcomplete(ax, redshift, redshift, 
                 restrict, isdbl, isdbl, name='redshift', 
                 split=False, nbins=10)
    ax2 = ax.twinx()
    ax2.hist(redshift[restrict], bins=10, 
            histtype='step', color='k', ls='dotted')
    ax.set_ylim((0,1.1))
    
    ax = subs.next()
    plotcomplete(ax, sepkpc, np.array([m.measSep12() for m in mergers]), 
                 (fluxratio>0.4)&zestcut&(masscut), 
                 np.array([m.isdbl for m in mergers]),
                np.array([m.isdbl for m in mergers]), name='separation',
                split=True, nbins=[0,1.0,1.75,2.25,2.75,4.0,6.0,9.0,11.0])
                
    ax.set_ylim((0,1.1))
           
    ax = subs.next()
    plotcomplete(ax, fluxratio, np.array([m.measFlux12() for m in mergers]), 
                 (sepkpc>2.2)&(sepkpc<8.0)&zestcut&(masscut), 
                 isdbl, isdbl, name='flux ratio',
                split=True, nbins=10)
    ax2 = ax.twinx()
    ax2.hist(fluxratio, histtype='step', bins=10, color='k', ls='dotted')
#    print np.percentile(np.array([m.measFlux12() 
#        for m in mergers[(fluxratio<0.3)&restrict&isdbl]])/
#    fluxratio[isdbl&restrict&(fluxratio<0.3)], (10,25,50,75,90))
    #sys.exit()
    ax.set_ylim((0,1.1))
                
    ax = subs.next()
    ax.scatter(fluxratio[isdbl&masscut&zestcut], 
               np.array([m.measFlux12() for m in mergers[isdbl&masscut&zestcut]]),
               c=redshift[isdbl&masscut&zestcut], lw=0, alpha=0.5)
    ax.tick_params(labelsize=8)
    ax.set_xlabel('input flux ratio', size=10)
    ax.set_ylabel('meas. flux ratio', size=10)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim((0.005,1.1))
    ax.set_ylim((0.005,1.1))
    ax.plot([1.e-3,1],[1.e-3,1], 'k-')
    print 'fit to real flux-measured flux ratio',
    print np.polyfit((fluxratio[isdbl&masscut&zestcut]),np.array([m.measFlux12() 
               for m in mergers[isdbl&masscut&zestcut]]),1)
    print 'spearman-r real flux-measured flux ratio',
    print spearmanr(fluxratio[isdbl&masscut&zestcut], 
               np.array([m.measFlux12() 
               for m in mergers[isdbl&masscut&zestcut]]))
    
    ax = subs.next()
    plotcomplete(ax, mag, mag, 
                 restrict, isdbl, isdbl, name='magnitude',
                nbins=10, split=False)
    ax2=ax.twinx()
    ax2.hist(mag[restrict], bins=10,
            histtype='step', color='k', ls='dotted')
    ax.set_ylim((0,1.1))
    
    ax = subs.next()
    ax.scatter(sepkpc[isdbl&zestcut&masscut], 
               np.array([m.measSep12() for m in mergers[isdbl&zestcut&masscut]]),
               c=mag[isdbl&masscut&zestcut], lw=0, alpha=0.5)
    ax.tick_params(labelsize=8)
    ax.set_xlabel('input sep kpc', size=10)
    ax.set_ylabel('meas. sep kpc', size=10)
    ax.plot([0,11],[0,11], 'k-')
    

    fig.savefig(args.impath+'z_complete_2lt_off_lt5'+ending)
    
    
    #contamination plots
    fig = figure.Figure((16,16))
    canv = FigCanvas(fig)
    subs = makesubplots(fig, 4, 4)
    
    origisdbl = isdbl
    for m in mergers:
        m.cleanPeaks(centcut=10.0, totfluxcut=tfc)#, fluxRcut=0.2)
        m.assignRealPeaks()
    isdbl = np.array([m.isdbl for m in mergers])
    fakedbl = np.array([(len(m.peaklist)>1)&(not m.isdbl) for m in mergers])
    alldbl = fakedbl | isdbl
    restrict = zestcut & masscut & (redshift < 1.0)

    ax = subs.next()
    plotcomplete(ax, redshift[restrict], redshift[restrict], 
                 alldbl[restrict], isdbl[restrict], 
                fakedbl[restrict], nbins=10, split=False,
                name='redshift')
    plotcomplete(ax, redshift[restrict], redshift[restrict], alldbl[restrict],
                 isdbl[restrict], isdbl[restrict], nbins=10, split=False,
                name='redshift', bw=True)

    ax = subs.next()
    plotcomplete(ax, mag[restrict], mag[restrict], 
                 alldbl[restrict], isdbl[restrict], 
                fakedbl[restrict], nbins=10, split=False,
                name='magnitude')
    plotcomplete(ax, mag[restrict], mag[restrict], alldbl[restrict],
                 isdbl[restrict], isdbl[restrict], nbins=10, split=False,
                name='magnitude', bw=True)
                
    ax = subs.next()
    plotcomplete(ax, mass[restrict], mass[restrict], 
                 alldbl[restrict], isdbl[restrict], 
                fakedbl[restrict], nbins=10, split=False,
                name='log mass')
    plotcomplete(ax, mass[restrict], mass[restrict], alldbl[restrict],
                 isdbl[restrict], isdbl[restrict], nbins=10, split=False,
                name='log mass', bw=True)

    ax = subs.next()
    allseps, isdblSep, fakedblSep = map(list,
        zip(*[[ss, m.isdbl, (not m.isdbl)&(len(m.peaklist)>1)] 
        for m in mergers[restrict] for ss in m.getSepsMaxPeak()]))
    alldblSep = np.asarray(fakedblSep) | np.asarray(isdblSep)
    allseps = np.asarray(allseps)
    plotcomplete(ax, allseps, allseps, alldblSep, np.asarray(isdblSep), 
                 np.asarray(fakedblSep), nbins=np.linspace(0,12,15),
                 split=False)
    plotcomplete(ax, allseps, allseps, alldblSep, np.asarray(isdblSep), 
                 np.asarray(isdblSep), nbins=np.linspace(0,12,15),
                 split=False, bw=True, name='separation [kpc]')
    
    ax = subs.next()
    ax.hist(np.array([np.asarray(m.flux12())/m.flux
                 for m in mergers[restrict]]).ravel(),
            bins=np.logspace(-4,1,30), histtype='step', label='real peaks')
    ax.hist(np.array([np.asarray(m.flux12())/m.flux
                 for m in mergers[restrict&isdbl]]).ravel(),
            bins=np.logspace(-4,1,30), histtype='step', label='peaks in pairs')
    ax.hist([p.flux/m.flux for m in mergers[restrict] for p in m.extrapeaks],
            bins=np.logspace(-4,1,30), histtype='step', label='all peaks')
    ax.set_xlabel('peak flux/total flux')
    ax.legend(loc=2, prop={'size':8})
    ax.set_xscale('log')
    
    zall, magall, frall = map(np.asarray, zip(*[[m.z, m.mag, 
                                           p.flux/m.flux] 
                    for m in mergers[restrict] for p in m.extrapeaks]))
    zmerge = [zz for m in mergers[restrict] for zz in [m.z]*2]
    magmerge = [mm for m in mergers[restrict] for mm in [m.mag]*2]
    fmerge = np.array([np.array(m.flux12())/m.flux 
                        for m in mergers[restrict]]).ravel()
    zmerge = np.asarray(zmerge)[fmerge>0]
    magmerge = np.asarray(magmerge)[fmerge>0]
    fmerge = fmerge[fmerge>0]

    zalldbl, magalldbl, fralldbl = map(np.asarray, zip(*[[m.z, m.mag, 
                                           p.flux/m.flux] 
                    for m in mergers[restrict&(alldbl)] for p in m.extrapeaks]))
    zmergedbl = [zz for m in mergers[restrict] for zz in [m.z]*2]
    magmergedbl = [mm for m in mergers[restrict] for mm in [m.mag]*2]
    fmergedbl = np.array([np.array(m.flux12())/m.flux 
                        for m in mergers[restrict&(alldbl)]]).ravel()
    zmergedbl = np.asarray(zmerge)[fmerge>0]
    magmergedbl = np.asarray(magmerge)[fmerge>0]
    fmergedbl = fmerge[fmerge>0]
    
    ax = subs.next()
    ax.plot(zall, frall, 'k.')
    ax.plot(zmerge, fmerge, 'c,')
#    for mcurr in [18, 19.5, 21, 22.5]:
#        ax.plot(np.arange(0.2,1.,0.1), 0.04 - 0.02*(mcurr-20.5) + 
#            0.04*(np.arange(0.2,1,0.1)-0.5), 'm-')
    ax.set_ylabel('peak flux/total flux')
    ax.set_yscale('log')
    ax.set_ylim((1.e-3,1.))
    #print '80 percentiles with redshift of non-detections ',
    zbs = np.linspace(0.2,0.9,7)
    for izbs  in range(len(zbs)-1):
        bad80 = np.percentile(np.asarray(frall)[(zall>zbs[izbs]) & 
        (zall<=zbs[izbs+1])], 75)
#        print '%.1f -- %.2e, %.2f'%(zbs[izbs], bad80,
#                                    percentileofscore(fmerge[(zmerge>zbs[izbs])&(
#                                    zmerge<=zbs[izbs+1])], bad80)),
        ax.plot(np.mean(zbs[izbs:izbs+2]), bad80, 'ro')
        good25= np.percentile(np.asarray(fmerge)[(zmerge>zbs[izbs]) & 
        (zmerge<=zbs[izbs+1])], 25)
        #print '%.2f -- %.2e, %.2f'%(zbs[izbs], good25,
#                                    percentileofscore(np.asarray(frall)[(zall>zbs[izbs])&(
#                                    zall<=zbs[izbs+1])], good25)),
        ax.plot(np.mean([zbs[izbs:izbs+2]]), good25, 'y^')
    print ''
    
    ax = subs.next()
    cutoff=[]
    for izbs in range(len(zbs)-1):
        inz_a = (zall > zbs[izbs])&(zall <= zbs[izbs+1])
        inz_m = (zmerge > zbs[izbs])&(zmerge <= zbs[izbs+1])
        bins = np.logspace(-3,1,40)
        bad = np.cumsum(np.histogram(frall[inz_a], bins=bins)[0][::-1])
        allp = bad+np.cumsum(np.histogram(fmerge[inz_m], bins=bins)[0][::-1])
        ax.semilogx(bins[:-1][::-1], bad*1.0/allp)
        xval = bins[:-1][::-1]
        keep = (xval > 0.002) & (xval < 0.3)
        cutoff.append(np.interp(0.25, (1.0*bad/allp)[keep], xval[keep]))
        bad2 = np.cumsum(np.histogram(frall, bins=bins)[0][::-1])
        allp2 = 0.001+bad2 + np.cumsum(np.histogram(fmerge, bins=bins)[0][::-1])
        ax.semilogx(xval, bad2*1.0/allp2, 'k-')
        #print np.interp(0.25, (1.0*bad2/allp2)[keep], xval[keep])
    #print zbs, cutoff
    totcutpoly = np.polyfit((zbs[:-1]-0.5), cutoff, 1)
    #print totcutpoly
    ax.set_xlim((0.002,0.3))
    ax.tick_params(labelsize=9)
    ax.set_xlabel('peak flux ratio to galaxy')
    ax.set_ylabel('cum. frac. (p/t > x) of non-real peaks')
    print np.mean(cutoff), np.median(cutoff)
    print 'total fraction of non-real peaks with cutoff at %.3f'%tfc,
    print '%.3e'%(len(frall[frall>tfc])*1.0/(len(frall[frall>tfc])+len(fmerge[fmerge>tfc])))
    print 'total fraction of non-real peaks (in pairs) with cutoff at %.3f'%tfc,
    print '%.3e'%(len(fralldbl[fralldbl>tfc])*1.0/(len(fralldbl[fralldbl>tfc])+len(fmergedbl[fmergedbl>tfc])))
    print 'fraction of real peaks above cutoff=%.3f'%tfc,
    print '%.3e'%(len(fmerge[fmerge>tfc])/(2.0*len(mergers[restrict])))

    ax = subs.next()
    ax.plot(magall, frall, 'k.')
    ax.plot(magmerge, fmerge, 'c,')
#    for zcurr in [0.2,0.45,0.75,1.0]:
#        ax.plot(np.arange(17,23,0.2), 0.04 - 0.02*(np.arange(17,23,0.2)-20.5) + 
#            0.04*(zcurr-0.5), 'm-')
    ax.set_ylabel('peak flux/total flux')
    ax.set_yscale('log')
    ax.set_ylim((1.0e-3,1.0))
    #print '80 percentiles with magnitude of non-detections ',
    zbs = np.linspace(18.5,20,5)
    for izbs  in range(len(zbs)-1):
        bad80 = np.percentile(np.asarray(frall)[(magall>zbs[izbs]) & 
        (magall<=zbs[izbs+1])], 75)
        #print '%.2f -- %.2e, %.2f'%(zbs[izbs], bad80,
#                                    percentileofscore(fmerge[(magmerge>zbs[izbs])&
#                                    (magmerge<=zbs[izbs+1])], bad80)),
        ax.plot(np.mean(zbs[izbs:izbs+2]), bad80, 'ro')
        good25= np.percentile(np.asarray(fmerge)[(magmerge>zbs[izbs]) & 
        (magmerge<=zbs[izbs+1])], 25)
        #print '%.2f -- %.2e, %.2f'%(zbs[izbs], good25,
#                                    percentileofscore(np.asarray(frall)[(magall>zbs[izbs])&(
#                                    magall<=zbs[izbs+1])], good25)),
        ax.plot(np.mean(zbs[izbs:izbs+2]), bad80, 'ro')
        ax.plot(np.mean(zbs[izbs:izbs+2]), good25, 'y^')
    #print ''
    np.savetxt('foo_all.txt', np.asarray([zall, magall, frall]).T, fmt='%.2f %.2f %.3e')
    np.savetxt('foo_mergeisdbl.txt', np.asarray([zmerge, magmerge, fmerge]).T,
               fmt='%.2f %.2f %.3e')
    
    ax = subs.next()
    allfr, isdblfr, fakedblfr = map(list,
        zip(*[[ss, m.isdbl, (not m.isdbl)&(len(m.peaklist)>1)] 
        for m in mergers[restrict] for ss in m.getFluxRatioMaxPeak()]))
    alldblfr = np.asarray(fakedblfr) | np.asarray(isdblfr)
    allfr = np.asarray(allfr)
    plotcomplete(ax, allfr, allfr, alldblfr, np.asarray(isdblfr), 
                 np.asarray(fakedblfr), nbins=15,
                 split=False)
    plotcomplete(ax, allfr, allfr, alldblfr, np.asarray(isdblfr), 
                 np.asarray(isdblfr), nbins=15,
                 split=False, bw=True, name='flux ratio to max peak')
    
#    ax = subs.next()
#    magall, zall, frmaxall = map(list, zip(*[[m.mag, m.z, ss]
#        for m in mergers[restrict] for ss in m.getFluxRatioMaxPeak()]))
#    ax.plot(magall, frmaxall, 'k.')
#    ax.plot([m.mag for m in mergers[restrict]],
#            [m.measFlux12() for m in mergers[restrict]], 'c.')
#    ax.set_ylim((1.e-4,1.1))
#    ax.set_yscale('log')
#    ax.set_ylabel('peak flux ratios', size=9)
    
    ax = subs.next()
    ax.hist([m.measFlux12() for m in mergers[restrict]],
            bins=np.logspace(-4,0.1,30), histtype='step', label='real peaks')
    ax.hist([ff for m in mergers[restrict]
            for ff in m.getFluxRatioMaxPeak(extraonly=True)],
            bins=np.logspace(-4,0.1,30), histtype='step', label='all peaks')
    ax.set_xlabel('peak flux/max peak flux')
    ax.legend(loc=2, prop={'size':8})
    ax.set_xscale('log')
    
    for m in mergers:
        m.cleanPeaks(centcut=10.0, totfluxcut= tfc, fluxRcut=cut_fr)
        m.assignRealPeaks()
    isdbl = np.array([m.isdbl for m in mergers])
    fakedbl = np.array([(len(m.peaklist)>1)&(not m.isdbl) for m in mergers])
    alldbl = fakedbl | isdbl
    minordbl = np.array([(m.fluxratio < desired_fr) for m in mergers])
    restrict = zestcut & masscut
    seps = np.array([m.measSep12() if isdbl[restrict&alldbl][im] 
            else m.getSepsMaxPeak()[0] for im,m in enumerate(mergers[restrict&alldbl])])
    frs = np.array([m.measFlux12() if isdbl[restrict&alldbl][im] else m.getFluxRatioMaxPeak()[0]
                    for im, m in enumerate(mergers[restrict&alldbl])])
    print 'contaim from fake sources with MP cut',
    print len(mergers[restrict&fakedbl])*1.0/len(mergers[restrict&alldbl])
    print 'completeness of real, major mergers with MP cut',
    print len(mergers[restrict&isdbl&(~minordbl)])*1.0/len(mergers[restrict&(~minordbl)])
    print 'contaim from minor mergers with MP cut',
    print len(mergers[restrict&isdbl&minordbl])*1.0/len(mergers[restrict&alldbl])
    print 'completeness for all mergers with MP cut',
    print len(mergers[restrict&isdbl])*1.0/len(mergers[restrict])          
    
    
    ax = subs.next()
    intlpr = np.array([m.isdbl & ((m.fluxratio < desired_fr)) 
                for m in mergers])
    plotcomplete(ax, redshift[restrict&alldbl], redshift[restrict&alldbl],
                 alldbl[restrict&alldbl], isdbl[restrict&alldbl], intlpr[restrict&alldbl],
                split=False, nbins=10, name='redshift')
    
    ax = subs.next()
    plotcomplete(ax, seps,
                np.array([m.sepkpc for m in mergers])[restrict&alldbl],
                alldbl[restrict&alldbl], isdbl[restrict&alldbl], intlpr[restrict&alldbl],
                split=True, nbins=np.linspace(2.0,10.,10), name='separation [kpc]')
    
    
    ax = subs.next()
    plotcomplete(ax, frs,
                np.array([m.fluxratio for m in mergers])[restrict&alldbl],
                alldbl[restrict&alldbl], isdbl[restrict&alldbl], intlpr[restrict&alldbl],
                split=True, nbins=10, name='flux ratio')
    #ax.set_xscale('log') 
    
#    ax = subs.next()
#    ax.hist(measfluxratio[restrict&(fluxratio>desired_fr)],
#             bins=np.linspace(0,1,25), histtype='step', color='c',
#            label='pairs w/ 2/1 > %.2f'%desired_fr)
#    print 'percentiles measured flux ratios with real ones larger than %.2f'%desired_fr,
#    print np.percentile(measfluxratio[restrict&(fluxratio>desired_fr)&isdbl&(origisdbl)], (10,25,50,75,90))
#    print 'percentiles measured flux ratios with real ones less than %.2f'%desired_fr,
#    print np.percentile(measfluxratio[restrict&(fluxratio<=desired_fr)&isdbl&(origisdbl)], (10,25,50,75,90))
#    ax.hist(measfluxratio[restrict&(fluxratio<desired_fr)],
#             bins=np.linspace(0,1,25), histtype='step', color='b',
#            label='pairs w/ 2/1 < %.2f'%desired_fr)
#    ax.legend(loc=2, prop={'size':8})
#    ax.set_xlabel('flux ratio')
    
    rlow = measfluxratio[restrict&(fluxratio<=desired_fr)&origisdbl]
    rReallow = measfluxratio[restrict&(fluxratio<=0.2)&origisdbl]
    rhigh = measfluxratio[restrict&(fluxratio>desired_fr)&origisdbl]
    frbins = np.linspace(0,1,40)
    fraclow = np.cumsum(np.histogram(rlow, bins=frbins)[0][::-1]) * 1.0 / \
    np.cumsum(np.histogram(np.r_[rhigh, rlow], bins=frbins)[0][::-1])
    xval = frbins[::-1][1:]
    keep = (xval < 0.5)
    low30 = np.interp(0.3, fraclow[keep], xval[keep])
    ax = subs.next()
    ax.plot(frbins[::-1][1:], fraclow)
    ax.plot(frbins[::-1][1:], [percentileofscore(rhigh, bbb)*1e-2 for bbb in frbins[::-1][1:]])
    print low30, percentileofscore(rhigh, low30)
    #print zip(xval[keep], fraclow[keep])
    frcut=0.25
    print frcut, percentileofscore(rhigh, frcut), np.interp(frcut, xval[keep][::-1],
                                                fraclow[keep][::-1]),
    print len(measfluxratio[(measfluxratio>frcut)&restrict&(fluxratio<=frcut)&origisdbl])* \
        1.0/(len(measfluxratio[(measfluxratio>frcut)&restrict&(fluxratio>frcut)&origisdbl]))
    ax.set_ylim((0,0.4))
    
    if False:
        #plot the completeness and contamination for different cuts
        ax1 = subs.next()
        ax2 = subs.next()
        
        for cc in [0.0, 0.1, 0.2, 0.25, 0.35, 0.4, 0.5]:
            for m in mergers:
                m.cleanPeaks(centcut=10.0, distcut=(2.2,8.0), totfluxcut= tfc, fluxRcut=cc)
                m.assignRealPeaks()
            isdbl = np.array([m.isdbl for m in mergers])
            measfluxratio = np.array([m.measFlux12() for m in mergers])
            ax1.plot(*(cumComplete(fluxratio[restrict], isdbl[restrict], bins=12, cum=True)), label='%.2f'%cc)
            ax2.plot(*(cumContaim(fluxratio[isdbl], measfluxratio[isdbl], bins=12, cum=True)))
            print np.median(measfluxratio[restrict&isdbl]), np.median(fluxratio[restrict&isdbl])
        ax1.legend(loc=2, prop={'size':6}, ncol=2, borderpad=0.0, columnspacing=0.5)
        ax1.set_ylabel('completeness', size=9)
        ax2.set_ylabel('contaim frac', size=9)
        ax1.set_xlabel('real flux ratio')
        ax2.set_xlabel('meas flux ratio')
        ax1.tick_params(labelsize=6)
        ax2.tick_params(labelsize=6)
    
#    print low30, percentileofscore(rhigh, low30)
#    print zip(xval[keep], fraclow[keep])
#    print 0.2, percentileofscore(rhigh, 0.2), np.interp(0.2, xval[keep][::-1], 
#                                                fraclow[keep][::-1])
#    per = 25
#    medfr = np.percentile(fluxratio[restrict&(isdbl)],per)
#    print '(1) %d percentile real fr in observed sample'%per, medfr,
#    print 'median ',np.median(fluxratio[restrict&isdbl])
#    intl2 = np.array([m.isdbl & (fluxratio < medfr) for m in mergers])
#    rlow2 = measfluxratio[restrict&(fluxratio<=medfr)&isdbl]
#    rhigh2 = measfluxratio[restrict&(fluxratio>medfr)&isdbl]
#    print 'contamination of mergers above median', len(rlow2)*1.0/(len(rlow2)+len(rhigh2))
#    print 'completeness above median', len(measfluxratio[restrict&(fluxratio>medfr)&isdbl])*1.0 / \
#    len(measfluxratio[restrict&(fluxratio>medfr)&origisdbl])
    fig.savefig(args.impath+'contaim_2lt_off'+ending)



    return 0


if __name__=='__main__':
    main()
    
