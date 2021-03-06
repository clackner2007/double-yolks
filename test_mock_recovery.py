#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################
# PROGRAM NAME
#
#
##############################
"""
Created on Wed Mar  4 17:07:07 2015

@author: clackner

program details
"""
import sys, os, re
import argparse
import numpy as np

import matplotlib.figure as figure
from matplotlib.backends.backend_ps import FigureCanvasPS as FigCanvasPS
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvasA
from mock_recovery import *

from scipy.stats import spearmanr, percentileofscore


def main():
    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("infolder", help='folder with peak finding output')
    parser.add_argument('origfile', help='file with real merger info')
    parser.add_argument('-o', '--outpath', dest='impath', default='./',
                        help='path to image output')
    parser.add_argument("-e", "--eps",
                 action='store_true', default=False, help='make eps plots',
                 dest='epsPlot')
    parser.add_argument("-m", "--images", nargs=1, dest='nimages', 
                        help='plot images of mocks with detected peaks\n(default=0)',
                        default=0)
    args = parser.parse_args()

    try:
        os.makedirs(args.impath)
    except OSError:
        pass
    
    FigCanvas = FigCanvasPS if args.epsPlot else FigCanvasA
    ending='.eps' if args.epsPlot else '.png'

    desired_fr = 0.25 #desired flux ratio we want to be sensitive down to (exclude 'minor' mergers)
    sep_kpc_min = 3.5 #minimum REAL separation we expect to measure
    sep_kpc_max = 8.0 #maximum REAL separation we expect to measure 
                        #(note the sims aren't separated by more than 10 kpc)
    
    mergers = np.asarray(readMergers(args.infolder, args.origfile))

    #plot images for the first nimages mocks you could imagine only doing certain kinds here, etc.
    if args.nimages:
        plotPeaks(mergers, os.path.dirname(args.origfile)+'/imgs/', 
                  args.impath, ending, FigCanvas, ngal=args.nimages)
    


    #what follows is a list of possible tests looking at the galaxies
    #i've included tests with different zest (morphology parameters)
    #but you can imagine doing lots of things if you've stored those parameters in
    #the simulatedSample file, they'll be read in in the mergers
    redshift = np.array([m.z for m in mergers])
    sepkpc = np.array([m.sep_kpc for m in mergers])
    fluxratio = np.array([m.flux_ratio for m in mergers])
    zest1 = np.array([m.zest1 for m in mergers])
    zest2 = np.array([m.zest2 for m in mergers])
    mag = np.array([-2.5*np.log10(m.flux)+config.zeropt for m in mergers])
    measfluxratio = np.array([m.measFlux12() for m in mergers])
    isdbl = np.array([m.isdbl for m in mergers])
    zestcut = np.in1d(zest1, [1,2]) & np.in1d(zest2, [1,2])
    
    #only look at sources with a REAL separation btw. 2.2 and 8 kpc
    #and reasonable morphologies and flux ratios > desired_fr
    #NOTE: for the completenesses to be reasonable these shouldn't be more restrictive than
    #the cuts used to select the peaks
    restrict = ((sepkpc > sep_kpc_min) & (sepkpc<sep_kpc_max) & (fluxratio>desired_fr))
    print len(mergers[restrict])
    
    #show plots of the completeness as a function of lots of things
    fig = figure.Figure((12,12))
    fig.set_size_inches((14,10))
    fig.subplots_adjust(left=0.05, right=0.95, wspace=0.25, bottom=0.07,
                        top=0.95)
    canv = FigCanvas(fig)
    subs = makesubplots(fig, nx=2)
    ax = subs.next()
    #redshift completeness
    plotcomplete(ax, redshift, redshift, 
                 restrict, isdbl, isdbl, name='redshift', 
                 split=False, nbins=10)
    ax.legend()
    ax2 = ax.twinx()
    ax2.hist(redshift[restrict], bins=10, 
            histtype='step', color='k', ls='dotted')
    ax2.tick_params(labelsize=8)
    ax.set_ylim((0,1.1))
    
    ax = subs.next()
    #completeness as a fcn of separation
    plotcomplete(ax, sepkpc, np.array([m.measSep12() for m in mergers]), 
                 zestcut&(fluxratio>desired_fr), isdbl, isdbl, name='separation',
                 split=True, nbins=[0,1.0,1.75,2.25,2.75,4.0,6.0,9.0,11.0])

    ax.set_ylim((0,1.1))
           
    ax = subs.next()
    #completeness as a function of flux ratio
    plotcomplete(ax, fluxratio, np.array([m.measFlux12() for m in mergers]), 
                 restrict, 
                 isdbl, isdbl, name='flux ratio',
                split=True, nbins=10)
    ax2 = ax.twinx()
    ax2.hist(fluxratio, histtype='step', bins=10, color='k', ls='dotted')
    ax2.tick_params(labelsize=8)
    ax.set_ylim((0,1.1))
    
    #comparison of measured and real flux ratios
    ax = subs.next()
    ax.scatter(fluxratio[isdbl&restrict], 
               measfluxratio[isdbl&restrict],
               c=redshift[isdbl&restrict], lw=0, alpha=0.5)
    ax.tick_params(labelsize=8)
    ax.set_xlabel('input flux ratio', size=10)
    ax.set_ylabel('meas. flux ratio', size=10)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim((0.005,1.1))
    ax.set_ylim((0.005,1.1))
    ax.plot([1.e-3,1],[1.e-3,1], 'k-')
    if sum(isdbl&restrict) > 0:
        print 'fit to real flux-measured flux ratio',
        print np.polyfit((fluxratio[isdbl&restrict]),
                         measfluxratio[isdbl&restrict], 1)
        print 'spearman-r real flux-measured flux ratio',
        print spearmanr(fluxratio[isdbl&restrict], 
                        np.array([m.measFlux12() 
                                  for m in mergers[isdbl&restrict]]))
    
    #completeness as a function of magnitude
    ax = subs.next()
    plotcomplete(ax, mag, mag, 
                 restrict, isdbl, isdbl, name='magnitude',
                nbins=10, split=False)
    ax2=ax.twinx()
    ax2.hist(mag[restrict], bins=10,
            histtype='step', color='k', ls='dotted')
    ax.set_ylim((0,1.1))
    
    #meausred vs. real separation in kpc
    ax = subs.next()
    ax.scatter(sepkpc[isdbl&zestcut], 
               np.array([m.measSep12() for m in mergers[isdbl&zestcut]]),
               c=mag[isdbl&zestcut], lw=0, alpha=0.5)
    ax.tick_params(labelsize=8)
    ax2.tick_params(labelsize=8)
    ax.set_xlabel('input sep kpc', size=10)
    ax.set_ylabel('meas. sep kpc', size=10)
    ax.plot([0,11],[0,11], 'k-')
    
    
    fig.savefig(os.path.join(args.impath,'completeness_scatter'+ending))
    
 
    #do contamination tests here --in progress

    #images with 2 measured peaks
    measdbl = np.array([m.getNPeaks()>1 for m in mergers])
    #images where the 2 measured peaks ARE the real merging galaxies
    realdbl = np.array([m.isdbl for m in mergers])
    #images where the 2 measured peaks aren't the real galaxies (one could be)
    fakedbl = measdbl & (~realdbl)
    #images where the 2 measured peaks are the real galaxies, but the 
    #real flux ratios and real separations are outside the desired range
    interloperdbl = realdbl & (fluxratio <= desired_fr) & ((sepkpc < sep_kpc_min) | (sepkpc >= sep_kpc_max))
    print "number of mergers", len(mergers)
    print 'number measured pairs', sum(measdbl)
    print 'number real pairs', sum(realdbl)
    print 'number false pairs', sum(fakedbl)
    
    print "number of mergers (restricted)", sum(restrict)
    print 'number measured pairs', sum(measdbl&restrict)
    print 'number real pairs', sum(realdbl&restrict)
    print 'number false pairs', sum(fakedbl&restrict)
    print "number of real pairs with wrong separation or flux ratio", sum(interloperdbl)
    
    #plot showing the peak-flux/total-flux ratios for 'real' galaxy peaks and
    #other peaks found in the image. In this case, the 'real' peaks can come
    #from images in which both galaxies or only one was detected
    ratio_false = np.array([p.flux/m.flux for m in mergers 
                              for p in m.extrapeaks])
    ratio_real = np.array([ip/m.flux for m in mergers for ip in m.flux12()])

    #total flux cut we want to use:
    totfluxcut = 0.03

    fig = figure.Figure((3,3))
    canv = FigCanvas(fig)
    left=0.2
    bot=0.2
    top=0.96
    right=0.96
    width = (right-left)
    ax = fig.add_axes((left,bot,width,top-bot))
    hm = ax.hist(ratio_real,
                 bins=np.logspace(-4,1.3,30), histtype='step', label='real galaxies',
                 color='k', lw=2)[0]
    ha = ax.hist(ratio_false,
                 bins=np.logspace(-4,1.3,30), histtype='step', label='other peaks', 
                 color='r')[0]
    ax.axvline(totfluxcut, color='#A0A0A0', lw=2)
    ax.set_xlabel(r'$\mathrm{peak\ flux/total\ galaxy\ flux}$', size=11)
    ax.set_ylabel(r'$\mathrm{number}$', size=11)
    ax.legend(loc=2, prop={'size':9})
    ax.set_xscale('log')
    ax.set_xlim((1.e-4, 1.0))
    ax.tick_params(labelsize=9)
    print 'fraction of other peaks to all peaks above ratio=%.3f:   '%totfluxcut,
    print sum(ratio_false > totfluxcut) * 1.0/(sum(ratio_false > totfluxcut) + 
                                                 sum(ratio_real > totfluxcut))
    print 'fraction of real peaks above ratio=%.3f:   '%totfluxcut,
    print sum(ratio_real>totfluxcut) * 1.0 / (len(ratio_real))
    fig.savefig(os.path.join(args.impath, 'peak_tot_flux_ratio'+ending))

    return 0


if __name__=='__main__':
    main()
    
