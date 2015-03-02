#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################
# PROGRAM NAME
#
#
##############################
"""
Created on Sun Apr 20 13:48:47 2014

@author: clackner

program details
"""
import sys, os, re, copy
import argparse
import numpy as np
from cosm_distTable import CosmoDist
import pyfits
import sigma_clip

import matplotlib.figure as figure
from matplotlib.backends.backend_ps import FigureCanvasPS as FigCanvasPS
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvasA

#######################################
def moveImage(image, delta_r, phi, noiseadd):
    """
    offset an image by a given number of pixels in a certain direction
    don't worry about sinc-interpolation, do whole pixels only
    add noise to the blank pixels  by the amount noiseadd
    """
    
    delta_x = int(delta_r * np.cos(phi))
    delta_y = int(delta_r * np.sin(phi))
    newimage = copy.copy(image)
    newimage = np.roll(newimage, delta_x, axis=0)
    newimage = np.roll(newimage, delta_y, axis=1)
    
    resamp_range = lambda delta, dim: np.arange(delta) if delta > 0 \
                                        else np.arange(dim+delta,dim)
    newimage[resamp_range(delta_x, image.shape[0]),:] = \
        np.random.normal(scale=noiseadd, size=(np.abs(delta_x),image.shape[1]))
    newimage[:, resamp_range(delta_y, image.shape[1])] = \
        np.random.normal(scale=noiseadd, size=(image.shape[0], np.abs(delta_y)))
                                                    
    return (newimage, delta_x, delta_y)

########################
def getMatch(selected, gal_list, delz=0.03, zests=[1,2,3]):
    """
    returns a list of galaxy pairs, which match in 
    redshift to delz=0.03
    """
    pair = np.zeros_like(np.atleast_1d(selected))
    for inds, s in enumerate(np.atleast_1d(selected)):
        
        parents = np.where((np.abs(gal_list['z'] - gal_list['z'][s]) < delz) &
                            (gal_list['ident'] != gal_list['ident'][s]) & 
                            (np.in1d(gal_list['zest'], zests)))[0]

        if len(parents) < 1:
            print "Not enouch close z pairs for ", gal_list[s]
            pair[inds] = -1
        else:
            pair[inds] = np.random.permutation(parents)[0]
    return pair




##############################
def listStamps(infile, impath, outpath='.', mag_limit=20.5, ngal=2000):
    """
    makes a list of postage stamp galaxies
    """
    pixel_scale = 0.03
    gal_list = np.loadtxt(infile,
                          #unpack=True,
                          dtype=[('ident', int),
                                 ('z', float), ('mag', float), 
                                 ('filename', '|S80'), ('zest', int),
                                    ('logmass', float)],
                          usecols=(0,2,3,1,96,101))
 

    gal_list = gal_list[(gal_list['z'] > 0.2) &
                        (gal_list['z'] < 1.1)] #& ((gal_list['zest']==1) | (gal_list['zest']==2))]
    #only apply magnitude limit to selected function
    selected = np.random.permutation(np.where((gal_list['mag']<mag_limit))[0])[:ngal]
    #selected = np.random.randint(0, len(gal_list), ngal)
    paired = getMatch(selected, gal_list, zests=[1,2,3])
    
    
    #offsets in kpc
    offsets = np.array([0.5, 1.5, 2.0, 2.5, 3.0, 5.0, 7.0, 10.0])
    to_pix = 1.0e-3/CosmoDist().ang_dist(gal_list['z']) * 206265. / pixel_scale
    
    outfile = open(outpath+'simulatedSample_%.2f.dat'%mag_limit, 'w')
    outfile.write('#gal1_id gal2_id mag1 mag2 flux_ratio z1 z2 zest1 zest2 sep_kpc x01 y01 x02 y02 mass1 mass2\n')
    
    peakoutfile = open(outpath+'input_peakfilter_%.2f.dat'%mag_limit, 'w')

    for pair in range(len(selected)):
        gal1 = gal_list[selected[pair]]
        gal2 = gal_list[paired[pair]]
        
        fluxratio = 10.0**(-0.4*np.abs(gal1['mag']-gal2['mag']))
        magtot = -2.5*np.log10(10.0**(-0.4*gal1['mag']) +
                                10.0**(-0.4*gal2['mag']))
        
        try:
            img1 = pyfits.open(impath + gal1['filename'])[0].data
            img2 = pyfits.open(impath + gal2['filename'])[0].data
        except IOError:
            print 'either ', impath + gal1['filename'],' or ',
            print impath + gal2['filename'],'is missing'
    
    
        img2noise = sigma_clip.sigclip(img2)[1]

        for off in offsets:  
            img2_off, delx, dely = moveImage(img2, off*to_pix[selected[pair]],
                                       np.random.random()*2.0*np.pi, img2noise)
            compimg = img1 + img2_off
            hdu = pyfits.PrimaryHDU(compimg)
            hdu.writeto(outpath+"imgs/%06d_%06d_%04.1f.fits"%(gal1['ident'],
                                                   gal2['ident'], off),
                        clobber=True)
    
            outfile.write('%d %d %.2f %.2f %.2e %.2f %.2f %d %d %.1f %.1f %.1f %.1f %.1f %.2e %.2e\n' % (gal1['ident'], 
             gal2['ident'], gal1['mag'], gal2['mag'],
             fluxratio, 
             gal1['z'], gal2['z'], gal1['zest'], gal2['zest'], off,
            img1.shape[0]/2.0, img1.shape[1]/2.0, img1.shape[0]/2.0+delx,
            img1.shape[1]/2.0+dely, gal1['logmass'], gal2['logmass']))
            
            peakoutfile.write('%06d%06d %06d_%06d_%04.1f.fits %.2f %.3e\n'%(
            gal1['ident'],gal2['ident'], gal1['ident'], gal2['ident'],
            off, gal1['z'], magtot))
            
    outfile.close()
    peakoutfile.close()
                                   
                                            
if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputlist", help='input list of galaxies')
    parser.add_argument('-m', '--maglimit', help='magnitude limit of selected galaxies',
                        default=23.0, dest='maglimit', type=float)
    parser.add_argument('-i', '--impath', dest='impath',
                        default='/data/cosmos_data/photodata/imgs/',
                        help='path to input images')
               
    parser.add_argument("-e", "--eps",
                 action='store_true', default=False, help='make eps plots',
                 dest='epsPlot')
    args = parser.parse_args()

    FigCanvas = FigCanvasPS if args.epsPlot else FigCanvasA
    ending='.eps' if args.epsPlot else '.png'
    
    listStamps(args.inputlist, args.impath, mag_limit=args.maglimit,
               outpath=os.path.dirname(args.inputlist)+'/')
    
