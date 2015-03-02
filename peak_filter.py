# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 19:20:52 2012

@author: clackner
"""

import argparse
import sys,os,re
import pyfits
#import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import matplotlib.cm
import copy
import scipy.ndimage
import glob
from cosm_distTable import CosmoDist
import flood_fill as ff
import fluxSinc
from multiprocessing import Pool

def sidePlot(ax1,img):
    mymap = copy.copy(matplotlib.cm.jet)
    mymap.set_bad('w')
    ax1.imshow(img, origin='lower',cmap=mymap, interpolation='none')
    div=make_axes_locatable(ax1)
    axRight=div.append_axes("right",0.5, sharey=ax1)
    axRight.plot(np.sum(img,axis=1), np.arange(img.shape[1]))
    axTop = div.append_axes("top",0.5,sharex=ax1)
    axTop.plot(np.sum(img,axis=0))
    axTop.tick_params(labelleft='off', labelbottom='off')
    axRight.tick_params(labelleft='off', labelbottom='off')
    #plt.show()

def makeRing(nx,ny,r_in,r_out):
    """
    makes a ring filter footprint
    nx,ny=x,y dimensions of footprint
    r_in, r_out=inner/outer radii of ring
    the footprint is +1 in the ring and 0 outside
    """

    x,y = np.meshgrid(np.arange(ny),
                      np.arange(nx))
    x -= ny/2
    y -= nx/2
    rad = np.sqrt(x*x+y*y)
    footprint = (rad >= r_in)&(rad <= r_out)
    return footprint


def getSigma(img, clip=3.0, iterate=3):

    attempt = 0
    keep = np.ones_like(img, dtype='bool')
    while(attempt < iterate):
        sigma = np.std(img[keep].flatten())
        keep = np.abs(img) < clip*sigma
        attempt += 1
    return sigma


class Peak:
    """
    class containing peak and doing calculation on peak
    area is the x,y coordinate of pixels
    vals is the image values at those points
    can be switched out easily for another image
    """
    def __init__(self, pixels, values):
        self.area = pixels
        self.size = len(self.area)
        self.vals = values


    def getSize(self):
        return self.size

    def getFlux(self,img=None):
        if img==None:
            return np.sum(self.vals)
        else:
            return np.sum(img[ff.make_pix_array(self.area)])
    

    def getSincFlux(self, img, size=16.7, maskpeaks=[]):
        x, y = self.getCenter()
        mask = np.ones_like(img, dtype=float)
        for pmask in maskpeaks:
            mask[ff.make_pix_array(pmask.area)] *= 0.0
        return fluxSinc.getFlux( img*mask, x, y, 0.0, size, 0.0, 0.0 )


    def getCenter(self):
        x,y = ff.make_pix_array(self.area)
        Itot=self.getFlux()
        x0 = np.sum(self.vals*x)/Itot
        y0 = np.sum(self.vals*y)/Itot
        return (x0,y0)

    def getMoments(self):
        x,y = ff.make_pix_array(self.area)
        x0,y0 = self.getCenter()
        invItot = 1.0/self.getFlux()
        qxx = np.sum((x-x0)**2 * self.vals) * invItot
        qxy = np.sum((x-x0)*(y-y0)*self.vals) * invItot
        qyy = np.sum((y-y0)**2 * self.vals) * invItot

        temp = np.sqrt( (qxx-qyy)**2 + 4.0*qxy**2 )
        m1 = np.sqrt( 0.5 * (qxx + qyy + temp) )
        m2 = np.sqrt( 0.5 * (qxx + qyy - temp) )

        #for completeness, but unused
        phi = np.arctan( 2.0 * qxy / (qxx -qyy) ) * 0.5
        #print phi

        return m1, m2



def getPeaks(image, thresh):
    total_list = np.where(image > thresh)

    test_pix = set(ff.make_pix_list(total_list))
    areas = []

    while len(test_pix) > 0:
        tp = list(test_pix)[0]
        area = ff.floodfill(image,tp,
                            comp=lambda x: x > thresh,
                            diagonal=False)
        areas.append(Peak(area,image[ff.make_pix_array(area)]))
        test_pix = test_pix.difference(set(areas[-1].area))

    return areas
    
    
def mapMeasPeak(galobj):
    """
    mapping function that processes images for multithreading
    returns a list of peaks
    f is the filename
    """
    f = galobj['file']
    fwhm = 0.09/0.03
    fileid = f.split('/')[-1].split('_')[-3]#[-8] #real
    #fileid = f.split('/')[-1].split('_')[-1].split('.')[0] #fake mocks
    raS, decS = os.path.splitext(os.path.basename(f))[0].split("_")[-2:]#[-7:-5]
    #raS, decS = (fileid,fileid)
    
    if galobj['inkpc']:
        #kpc in pixels
        kpc_1 = 1.0/(galobj['z'])*206265.0/0.03
        if kpc_1 > 100.:
            kpc_1= 5*fwhm

    img=pyfits.open(f)[0].data
    #r_in = fwhm*1.5 if not args.inkpc else max(kpc_1,fwhm*1.3)
    r_in = fwhm*1.3 if not galobj['inkpc'] else max(kpc_1,fwhm*1.3)
    #r_in = fwhm*2.3 if not args.inkpc else max(kpc_1,fwhm*1.3)
    r_in = fwhm/2.35*3.0
    r_out = r_in + 0.7
    nx=np.round((r_out + 1.5)*2)
    lowpass=scipy.ndimage.median_filter(img,
                                        footprint=
                                        makeRing(nx,nx,r_in,r_out),
                                        mode='constant')

    new_img = img - lowpass
    smooth_img = new_img
    #get some measure of the noise
    sigma = getSigma(smooth_img, iterate=6)
    
    #find peaks (5 sigma peaks)
    peak_list = np.asarray(getPeaks(smooth_img, 5*sigma))
    #print len(peak_list)
    peak_list = peak_list[np.where(np.array([p.size>8
                            for p in peak_list]))[0]]
    
    galine = '%d %.2f %.2e %s %s '%(galobj['id'], galobj['z'], galobj['delz'], raS, decS)
    galine += '%.2f %s %s %s %d\n'%(galobj['mag'], fileid, raS, decS, len(peak_list))
                            

    pstring = ''
    for ip, p in enumerate(peak_list):
            pstring += '%d %s %s %d '%(galobj['id'], raS, decS, ip)
            ms=p.getMoments()
            pstring += '%d %.2e '%(p.size, (ms[1]/ms[0]))
            pstring +='%.7f '%p.getCenter()[0]+'%.7f '%p.getCenter()[1]
            pstring +='%.3e '%p.getFlux()
            pstring += '%.3e\n'%p.getFlux(img)
            
    return (galine, pstring)


###########################################################
############################################################
def main():

    parser=argparse.ArgumentParser()
    parser.add_argument('infile',  help='input file list')
    parser.add_argument('-p', '--path', help='path to output files',
                        default='./', dest='path')
    parser.add_argument('-k', '--kpc', help='do rings in kpc',
                        default=False, action='store_true',
                        dest='inkpc')
    parser.add_argument('-i', '--impath', dest='impath',
                        default='/data/cosmos_data/photodata/imgs/',
                        help='path to input images')
    parser.add_argument('-d', '--deltaz', default=0.0, type=float,
                        dest='delz',
                        help='redshift offset from original image')

    args=parser.parse_args()

    #get the files
#    id, ra, dec, zphot = np.loadtxt(args.infile,
#                                    unpack=True,
#                                    dtype=[('a', int),
#                                           ('b', float),
#                                           ('c', float),
#                                           ('e', float)],
#                                    usecols=(0,1,2,3))
    ident, zphot, mag, filename = np.loadtxt(args.infile,
                                     unpack=True,
                                     dtype=[('a', int),
                                            ('b', float), ('c', float),
                                            ('f', '|S80')],
                                            usecols=(0,2,3,1))
    try:
        delz = np.loadtxt(args.infile, dtype=float, usecols=(4,))
    except:
        delz = np.zeros_like(ident) + args.delz

    impath = args.impath #'/data/cosmos_data/photodata/imgs/'
    if args.delz != 0:
        files = [impath+'z%0.2f_'%(zphot[i])+filename[i]
                for i in range(len(ident))]
        for nf in range(len(files)):
            openOK = False
            try:
                with open(files[nf]): pass
                openOK =True
            except IOError:
                files[nf] = impath+'z%0.2f_'%(zphot[nf]+0.0001)+filename[nf]
            if not openOK:
                try:
                    with open(files[nf]): pass
                    openOK = True
                except IOError:
                    files[nf] = impath+'z%0.2f_'%(zphot[nf]-0.0001)+filename[nf]

    else:
        files = [impath+filename[i] for i in range(len(ident))]

#    files = [impath+'%06d_%.5f_%.5f.fits'%(id[i], ra[i], dec[i])
 #            for i in range(len(id))]
    if args.inkpc:
        z_list = np.loadtxt('/data/cosmos_data/id_ra_dec_z',
                            usecols=[1,2,3])
        cc = CosmoDist()
        ad = cc.ang_dist(z_list[:,2])*1000.0
        z_dict = dict(zip(zip(z_list[:,0], z_list[:,1]), ad))


    gals = np.array(zip(ident, files, zphot, mag, 
                        delz, np.repeat(args.inkpc, len(ident))), 
                        dtype=np.dtype([('id', int), ('file', '|S128'),
                        ('z', float), ('mag', float), ('delz', float), ('inkpc', bool)]))

        
    gal_file = open(args.path+"gal_list", 'w')
    peak_file = open(args.path+'peak_list', 'w')    

    #multi-thread
    pool = Pool(processes=12)
    results = pool.map(mapMeasPeak, gals[:], len(gals[:])/(12*2)+1)
    pool.close()
    pool.join()
    
    for r in results:
        gal_file.write(r[0])
        peak_file.write(r[1])
    

    #single thread    
#    for numf, f in enumerate(files[:1000]):     
#        galine, peaklines = mapMeasPeak(gals[numf])         
#        print >>gal_file, galine
#        print >>peak_file, peaklines

    peak_file.close()
    gal_file.close()

    return 0

if __name__=='__main__':
    main()
