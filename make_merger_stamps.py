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
import pyfits
import sigma_clip
import astropy.io

from astropy import table
from astropy.cosmology import FlatLambdaCDM
import astropy.wcs

from ConfigParser import SafeConfigParser
import configParams as config

#######################################
def makeWCS(image):
    """
    make a WCS centered on the central pixel in the image
    """
    sx, sy = image.shape
    wcs = astropy.wcs.WCS(naxis=2)
    wcs.wcs.crpix = np.array([sx/2.0, sy/2.0])
    wcs.wcs.crval = np.array([0.0,0.0])
    wcs.wcs.cunit = np.array(["pix", "pix"])
    #wcs.wcs.cd = np.array([[1.0,0.0],[0.0,1.0]])
    #wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    
    return wcs

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
        
        parents = np.where((np.abs(gal_list['Z'] - gal_list['Z'][s]) < delz) &
                            (gal_list['ID'] != gal_list['ID'][s]) & 
                            (np.in1d(gal_list['ZEST'], zests)))[0]

        if len(parents) < 1:
            print "Not enouch close z pairs for ", gal_list[s]
            pair[inds] = -1
        else:
            pair[inds] = np.random.permutation(parents)[0]
    return pair




##############################
def listStamps(infile, impath, outpath='.', H0=70, Om0=0.3, 
               mag_limit=20.5, ngal=2000, 
               zlims=[0.2,1.1], offsets = np.array([0.5, 1.5, 2.0, 2.5, 3.0, 5.0, 7.0, 10.0]), 
                morph_class=[1,2,3]):
    """
    makes a list of postage stamp galaxies
    """
    #input needs columns IDENT, FILENAME, Z, and MAG, optionally ZEST (MORPHOLOGY CLASS)
    gal_list = table.Table().read(infile,
                            format='fits' if os.path.splitext(infile)[1]=='.fits' else 'ascii')
    if 'ZEST' not in gal_list.colnames:
        gal_list.add_column(table.Column(name='ZEST', 
                                         data=np.ones(gal_list['ID'].shape)*morph_class[0]))
    gal_list = gal_list[(gal_list['Z'] > zlims[0]) &
                        (gal_list['Z'] < zlims[1])]
    #only apply magnitude limit to selected function
    selected = np.random.permutation(np.where((gal_list['MAG']<mag_limit))[0])[:ngal]
    paired = getMatch(selected, gal_list, zests=morph_class)
    print selected
    print paired
    
    #offsets in kpc
    cosmos = FlatLambdaCDM(H0=H0, Om0=Om0)  

    to_pix = 1.0e-3/cosmos.angular_diameter_distance(gal_list['Z']).value * 206265. / config.pixelscale
    
    outfile = open(outpath+'simulatedSample_%.2f.dat'%mag_limit, 'w')
    outfile.write('#ID gal1_id gal2_id mag1 mag2 flux_ratio z1 z2 ')
    outfile.write('zest1 zest2 sep_kpc x01 y01 x02 y02\n')
    
    peakoutfile = open(outpath+'input_peakfilter_%.2f.dat'%mag_limit, 'w')
    peakoutfile.write('ID FILENAME Z MAG\n')

    count = 0
    for pair in range(len(selected)):
        gal1 = gal_list[selected[pair]]
        if paired[pair] < 0:
            continue
        gal2 = gal_list[paired[pair]]
        
        fluxratio = 10.0**(-0.4*np.abs(gal1['MAG']-gal2['MAG']))
        magtot = -2.5*np.log10(10.0**(-0.4*gal1['MAG']) +
                                10.0**(-0.4*gal2['MAG']))
        
        try:
            img1 = astropy.io.fits.open(impath + gal1['FILENAME'])[0].data
            img2 = astropy.io.fits.open(impath + gal2['FILENAME'])[0].data
        except IOError:
            print 'either ', impath + gal1['FILENAME'],' or ',
            print impath + gal2['FILENAME'],'is missing'
    
        img2noise = sigma_clip.sigclip(img2)[1]

        for off in offsets:  
            img2_off, delx, dely = moveImage(img2, off*to_pix[selected[pair]],
                                       np.random.random()*2.0*np.pi, img2noise)
            compimg = img1 + img2_off
            wcs = makeWCS(compimg)
            hdu = astropy.io.fits.PrimaryHDU(compimg, header=wcs.to_header())
            hdu.writeto(outpath+"imgs/%06d_%06d_%04.1f.fits"%(gal1['ID'],
                                                   gal2['ID'], off),
                        clobber=True)
    
            outfile.write('%6d %d %d %.2f %.2f %.2e %.2f %.2f %d %d %.1f %.1f %.1f %.1f %.1f\n' % (
                count, gal1['ID'], gal2['ID'], gal1['MAG'], gal2['MAG'],
                 fluxratio, 
                 gal1['Z'], gal2['Z'], gal1['ZEST'], gal2['ZEST'], off,
                img1.shape[0]/2.0, img1.shape[1]/2.0, img1.shape[0]/2.0+delx,
                img1.shape[1]/2.0+dely))
            
            peakoutfile.write('%6d %06d_%06d_%04.1f.fits %.2f %.3e\n'%(
            count, gal1['ID'], gal2['ID'],
            off, gal1['Z'], magtot))
            count += 1
            
    outfile.close()
    peakoutfile.close()
                                   
                                            
if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputlist", help='input list of galaxies')
    parser.add_argument("param_file", help="config parameter file (*.ini) for making mocks")
    parser.add_argument('-o' '--outputpath', dest='outputpath', 
                        help='output path for images and list file',
                        default='')
    parser.add_argument('-i', '--impath', dest='impath',
                        help='path to input images')
    args = parser.parse_args()

    cParser = SafeConfigParser()
    cParser.read(args.param_file)
    
    config_mocks = dict([(name, cParser.getfloat('list_stamps', name))
                        for name in cParser.options('list_stamps')])
    for category in ('offsets', 'morph_class', 'zlims'):
        mylist = np.asarray([cParser.getfloat(category, n) for n in cParser.options(category)])
        config_mocks[category] = mylist                    
    
    listStamps(args.inputlist, args.impath,
               outpath=args.outputpath, H0=config.H0, Om0=config.Om0, **config_mocks)
    
