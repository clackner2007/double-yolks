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
import os
import numpy as np
from collections import defaultdict
import matplotlib.figure as figure
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator

import astropy.io.fits
import astropy.table
from astropy.cosmology import FlatLambdaCDM

import configParams as config

def readMergers(inputfolder, origfile, clean=False):
    """
    read in a set of mock mergers, including where the merger was finally detected
    """
    with open(os.path.join(inputfolder,'cleaned_peaks_sources.txt')) as pSrc, open(os.path.join(inputfolder+'cleaned_peaks.txt')) as pXY:
        linesSrc = pSrc.readlines()
        linesPxy = pXY.readlines()
        peaks = defaultdict(list)
        for i in range(1,len(linesSrc)):
            fields = linesSrc[i].split()
            galid = int(fields[0])
            npeaks = int(fields[5])
            for p in range(npeaks):
                x0 = linesPxy[i].split()[2+2*p]
                y0 = linesPxy[i].split()[2+2*p+1]
                flux = fields[6+npeaks+p]
                peaks[galid].append(map(float, (x0,y0,flux)))
    
    mergers = list()
    origGals = astropy.table.Table().read(origfile, format='ascii')
    cosmos = FlatLambdaCDM(H0=config.H0, Om0=config.Om0)
    ang_dist_kpc = cosmos.angular_diameter_distance(origGals['z1']).to('kpc').value
    for indg, g in enumerate(origGals):
        merger_kwargs = dict([(c, g[c]) for c in origGals.colnames])
        merger_kwargs['flux'] = (10.0**(-0.4*(g['mag1']-config.zeropt)) + 
                                10.0**(-0.4*(g['mag2']-config.zeropt)))
        merger_kwargs['sepas'] = merger_kwargs['sep_kpc']/ang_dist_kpc[indg] * 206265.
        merger_kwargs['angdist_kpc'] = ang_dist_kpc[indg]
        merger_kwargs['z'] = g['z1']
        mergers.append(Merger(**merger_kwargs))
        
        if g['ID'] in peaks:
            for p in peaks[g['ID']]:
                mergers[-1].addPeak(**dict(zip(['x0', 'y0', 'flux'], p)))
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
        self.mag = -2.5*np.log10(self.flux)+config.zeropt
        
    def setReal(self):
        self.isdbl = True
    
    def toKpc(self, arcsec):
        return arcsec / 206265.0 * self.angdist_kpc
        
    def toArcsec(self, kpc):
        return kpc / (self.angdist_kpc) * 206265.
        
    def getNPeaks(self):
        return len(self.peaklist)
        
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
        
        imgfile = os.path.join(path, '%06d_%06d_%04.1f.fits'%(self.id1, 
                                                              self.id2, self.sepkpc))
        data = astropy.io.fits.open(imgfile)[0].data
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
        print 'plotted ', self.id1, self.id2, self.sepkpc, self.z, self.flux_ratio,
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
    """
    makes images showing the peaks (real and fake) for mock mergers
    """
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
        fig.savefig(os.path.join(impath,'peaksall_%06d_%06d'%key+ending))


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

