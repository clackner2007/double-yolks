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
    
