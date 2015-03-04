#!/usr/bin/env python

##############################
# detections.py
#
# Claire Lackner
# 2012-11-29
#
# Class file for nuclei detections
##############################
"""program details"""

import os
import argparse
import numpy as np

from astropy import table
from collections import defaultdict

import matplotlib.figure as figure
from matplotlib.backends.backend_ps import FigureCanvasPS as FigCanvasPS
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvasA
from matplotlib.patches import Ellipse

from astropy.cosmology import FlatLambdaCDM 
import astropy.wcs
import pyfits
from ConfigParser import SafeConfigParser

import configParams as config



class Nuclei():

    def __init__(self,x0,y0,npix,flux,ellip,allflux):
        self.x0=x0
        self.y0=y0
        self.npix=npix
        self.flux=flux
        self.ellip=ellip
        self.allflux=allflux

    def getDist(self,x,y,scale=1.0):
        return np.sqrt((x-self.x0)**2 +
                       (y-self.y0)**2)*scale

class Galaxy():

    def __init__(self, ident, z, mag, ang_dist, ang_dist_orig, ra=None, dec=None, imgfile=None,
                 delz=0.0):
        self.ra=ra
        self.dec=dec
        self.imgfile = imgfile
        self.z=z
        self.ident = ident
        self.nuclei=[]
        self.ad = ang_dist
        self.kpc_p_pix = (self.ad*1000.0)/(206265.0)*config.pixelscale
        self.flux = 10**(-0.4*(mag - config.zeropt))
        self.mag = mag
        self.delz = delz
        self.ad_orig = ang_dist_orig
        self.goodMask = None
        try:
            self.wcs = astropy.wcs.WCS(self.imgfile)
        except:
            print "No WCS available in image"
            self.wcs = None

    def add_nuclei(self, n):
        self.nuclei.append(n)

    def num_nuclei(self, mask=True, **kwargs):
        if len(self.nuclei)==0:
            return 0
        elif mask:
            return len(np.asarray(self.nuclei)[self.good_nuclei(**kwargs)])
        else:
            return len(self.nuclei)

    def good_nuclei(self, ecut=0.2, dist_cut=(0.0,8.0),
                    flux_cut=0.25,
                    totflux_cut=0.03, #0.03, #0.2, #0.05,
                    cent_dist_cut=10.0,
                    center=(8/0.03*0.5,8/0.03*0.5), 
                    pearsonr_cut=0.5,
                    redo=False):
        if (self.goodMask is not None) & (~redo):
            return self.goodMask

            
        center = (8.0/0.03*0.5*self.ad_orig/self.ad,#ang_dist(self.z-self.delz)/self.ad,
                  8.0/0.03*0.5*self.ad_orig/self.ad)#ang_dist(self.z-self.delz)/self.ad)
        #distcent = np.array([nn.getDist(center[0], center[1], scale=0.03)
        #                    for nn in self.nuclei])
        mask = np.array([(nn.ellip >= ecut) &
                         (nn.allflux/self.flux >= totflux_cut) &
                         #max(totflux_cut - 
                             #0.02*(self.mag-20.5) + 
                             #0.04*(self.z-0.5), 0.01)) &
                         #fitting formulat from mock n=4 galaxies
                         #(totflux_cut )) & ##* (self.ad**2/ang_dist(1.0)**2))) &
                         (nn.getDist(center[0],center[1],
                                     scale=self.kpc_p_pix) < cent_dist_cut)
                         for nn in self.nuclei])

        if mask.any():
            nucs=np.asarray(self.nuclei)[mask]
            maxI = np.argmax([nn.allflux for nn in nucs])
            maxFlux = nucs[maxI].allflux
            xM, yM = nucs[maxI].x0, nucs[maxI].y0
            mask2 = np.zeros_like(mask, dtype=bool)
            for numNn, nn in enumerate(self.nuclei):
                if (np.abs(nn.allflux-maxFlux) < 1.0e-8):
                    mask2[numNn] = True
                else:
                    d=nn.getDist(xM,yM,scale=self.kpc_p_pix)
                    mask2[numNn] = (d < dist_cut[1]) & \
                    (d >= dist_cut[0]) & (nn.allflux/maxFlux>flux_cut)
            self.goodMask = mask & mask2
        else:
            self.goodMask = mask

        #if (distcent[self.goodMask] >= 1.0).all():
        #    self.goodMask = np.zeros(len(self.nuclei), dtype=bool)
        if np.abs(self.getPearsonR(masked=True)) > pearsonr_cut:
            self.goodMask = np.zeros(len(self.nuclei), dtype=bool)

        return self.goodMask

    ###############
    def getPearsonR(self,masked=True):
        """
        get's the Pearson linear correlation
        coefficient of the peak coordinates

        Galaxies with a line of peaks are
        probably edge-on spirals. This measure
        is only sensible for at least 3 peaks
        """

        if masked:
            gn = self.good_nuclei()
            if len(gn) < 1:
                nucs=[]
            else:
                nucs = np.asarray(self.nuclei)[gn]
        else:
            nucs = np.asarray(self.nuclei)
        if len(nucs) < 3:
            return 0

        x = np.array([nn.x0 for nn in nucs])
        y = np.array([nn.y0 for nn in nucs])
        return np.sum( (x - np.mean(x))*(y - np.mean(y)) ) / \
            np.sqrt( np.sum((x - np.mean(x))**2) *
                     np.sum((y - np.mean(y))**2) )



    def getMaxNuc(self,masked=True, **kwargs):
        if masked:
            nucs=np.asarray(self.nuclei)[self.good_nuclei(**kwargs)]
        else:
            nucs=np.asarray(self.nuclei)
        mymax=np.argmax([nn.allflux for nn in nucs])
        return mymax
        
        
    def plot(self, FigCanvas, ending, outdir=''):
        fig=figure.Figure((8,4))
        canv=FigCanvas(fig)

        ax1=fig.add_subplot(121, frameon=False)
        img = pyfits.open(self.imgfile)[0].data
        vmax = 14*np.std(img.flatten())
        ax1.imshow(np.arcsinh(img), origin='lower', cmap='gray',
                   vmax=np.arcsinh(vmax), aspect='equal',
                   interpolation='none')
        ax2=fig.add_subplot(122, sharex=ax1, sharey=ax1, frameon=False)
        ax2.imshow(np.arcsinh(img),origin='lower', cmap='gray',
                   vmax=np.arcsinh(vmax),aspect='equal',
                   interpolation='none')
        #show the nuclei
        if self.num_nuclei() > 0:
            nucs = np.asarray(self.nuclei)[self.good_nuclei()]
            imax=self.getMaxNuc(masked=True)
            ax2.scatter([p.y0 for p in nucs],
                        [p.x0 for p in nucs],
                        c='r', marker='x',
                        s=[p.flux/nucs[imax].flux*25 for p in nucs])

            ax2.add_artist(Ellipse((nucs[imax].y0, nucs[imax].x0),
                                   10.0/self.kpc_p_pix, 10./self.kpc_p_pix,
                                   0.0, ec='m', fc='None', lw=1))

        #scale bar
        ax1.plot([5,5+5.0/self.kpc_p_pix], [5,5], 'b-', lw=2)
        ax1.text(5+5.0/self.kpc_p_pix*0.2, 10, '5 kpc', color='b')
        for a in fig.get_axes():
            a.tick_params(bottom='off', top='off', left='off', right='off',
                          labelbottom='off', labelleft='off')
        fig.subplots_adjust(top=0.98,left=0.0,right=0.98,bottom=0.00,
                            hspace=0.01, wspace=-0.1)
        fig.savefig(outdir+'peaks_{0}'.format(self.ident)+ending)


###############################################################
def makehtml(gals, mylist, FigCanvas, ending, outpath):

    s  = "<html>\n"
    s += "<body>\n"
    s += "<table>\n"
    for i in mylist:
        gals[i].plot(FigCanvas,ending, outdir=outpath)
        peaks = sorted(np.asarray(gals[i].nuclei)[gals[i].good_nuclei()], 
                       key=lambda p: p.allflux,
                       reverse=True)
                        
        s += "<tr>\n"
        s2 = "<table>\n"
        s2 += "<tr><td>%s:</td><td>%s</td></tr>\n" % \
            ('COSMOS_ID',repr(gals[i].ident))
        s2 += "<tr><td>%s:</td><td>%.3f</td></tr>\n" % ('Z',gals[i].z)
        if gals[i].ra:
            s2 += "<tr><td>%s:</td><td>%.6f</td></tr>\n" % ('RA',gals[i].ra)
        if gals[i].dec:
            s2 += "<tr><td>%s:</td><td>%.6f</td></tr>\n" % ('DEC',gals[i].dec)
        s2 += "<tr><td>%s:</td><td>%.3e</td></tr>\n" % ('P1/P2', 
                                                        peaks[1].allflux/peaks[0].allflux)
        s2 += "<tr><td>%s:</td><td>%.3e</td></tr>\n" % ('P2/tot', 
                                                        peaks[1].allflux/gals[i].flux)
        s2 += "<tr><td>%s:</td><td>%d</td></tr>\n" % \
            ('NPEAK',gals[i].num_nuclei())
        s2 += "</table>\n"
        s += "<td>%s</td>\n" % (s2)
        s += "<td><img src={0}></td>\n".\
            format('peaks_{0}{1}'.format(gals[i].ident,ending))
        s += "</tr>\n"


    s += "</table>\n"
    s += "</body>\n"
    s += "</html>"

    return s



################################################################

def main():

    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputfile", help='input file table')
    parser.add_argument("file_gal", help="galaxy listing from peak_filter (gal_list)")
    parser.add_argument("file_peaks", help="peak listing from peak_filter (peak_list)")
    parser.add_argument("param_file", help="config parameter file (*.ini) for cleaning peaks")
    parser.add_argument("-p", "--path", default='./', help='path to output', dest='path')
    parser.add_argument('-x', '--pixels', dest='pixels', default=False, action='store_true',
                        help='output pair coordinates in clean_pairs.txt in pixels instead of ra/dec')
    parser.add_argument('-l', '--plot', dest='makePlots', default=False,
                        action='store_true', help='make images of peaks')
    parser.add_argument("-i","--imgpath", help='path to images, needed for plotting')
    parser.add_argument("-e", "--eps",
                        action='store_true', default=False, help='make eps plots', dest='epsPlot')
    args = parser.parse_args()
    
    cParser = SafeConfigParser()
    cParser.read(args.param_file)
    gal_kwargs = dict([(name, cParser.getfloat('galaxy_init', name)) 
                        for name in cParser.options('galaxy_init')])


    gals=[]

    g1 = table.Table().read(args.file_gal, format='ascii')
    p1 = table.Table().read(args.file_peaks, format='ascii')

    cosmos = FlatLambdaCDM(H0=config.H0, Om0=config.Om0)  
    
    g1.add_column(table.Column(name='ang_dist', data=cosmos.angular_diameter_distance(g1['Z']).value))
    g1.add_column(table.Column(name='ang_dist_orig', 
                               data=cosmos.angular_diameter_distance(g1['Z']-g1['DELTA_Z']).value))
                              
    extra = table.Table().read(args.inputfile, 
                                format='fits' if os.path.splitext(args.inputfile)[1]=='.fits' else 'ascii')
    if ("RA" in extra.colnames) and ("DEC" in extra.colnames):
        g1.add_column(extra["RA"])
        g1.add_column(extra["DEC"])
        
    g1.add_column(table.Column(name="FILENAME", data=map(lambda x: args.imgpath+"/"+x, extra["FILENAME"])))
 
    dict_g = dict(zip(zip(g1['ID'], g1['DELTA_Z']),g1))
    peak_dict = defaultdict(list)

    for e in p1:
        peak_dict[(e['ID'], e['DELTA_Z'])].append(e)

    for k in dict_g:
        gg=dict_g[k]
        gals.append(Galaxy(gg['ID'], gg['Z'], gg['MAG'], gg['ang_dist'], 
                           gg['ang_dist_orig'],
                           ra=gg['RA'] if 'RA' in gg.colnames else None,
                            dec=gg['DEC'] if 'DEC' in gg.colnames else None, 
                            delz=gg['DELTA_Z'], imgfile=gg['FILENAME'], **gal_kwargs))
        if k in peak_dict:
            for p in peak_dict[k]:
                gals[-1].add_nuclei(Nuclei(p['PEAK_X0_PIX'], p['PEAK_Y0_PIX'], p['PEAK_NPIX'],
                                           p['PEAK_FLUX_FILTERED'], p['PEAK_B-A'], p['PEAK_FLUX']))
                gals[-1].file_id = int(p[0])

    gals = np.asarray(gals)

    clean_params = dict((name, cParser.getfloat("good_nuclei", name)) for name in cParser.options("good_nuclei"))
    clean_params['dist_cut'] = (clean_params['dist_cut_1'], clean_params['dist_cut_2'])
    del clean_params['dist_cut_1']
    del clean_params['dist_cut_2']
    clean_params['center'] = (clean_params['imsize_x']*0.5/config.pixelscale,
                             clean_params['imsize_y']*0.5/config.pixelscale)
    del clean_params['imsize_x']
    del clean_params['imsize_y']
    
    n_nucs = np.array([g.num_nuclei(**clean_params) for g in gals])

    in_pairs = np.where(n_nucs >= 2)[0]

    f = open(args.path+'cleaned_peaks_sources.txt', 'w')
    f1 = open(args.path+'all_peaks_sources.txt', 'w')
    fpeaks = open(args.path+'cleaned_peaks.txt', 'w')

    fpairs = open(args.path+"clean_pairs.txt", 'w')    
    
    print >>f1, "ID Z RA DEC PearonRho N_PEAK DIST_FROM_BRIGHTEST(xN_PEAK) FLUX(xN_PEAK)"
    print >>f, "ID Z RA DEC PearonRho N_PEAK DIST_FROM_BRIGHTEST(xN_PEAK) FLUX(xN_PEAK)"
    print >>fpeaks, "ID Z (PEAK_X0_PIX PEAK_Y0_PIX)xN_PEAK"
    print >>fpairs, "ID Z RA DEC PearsonRho N_PEAK SEP_KPC FLUX_RATIO",
    
    doRaDec = (~args.pixels) & (sum([g.wcs is None for g in gals])==0)
    if doRaDec:
        print >>fpairs, "RA_P0 DEC_P0 RA_P1 DEC_P1 FLUX_P0 FLUX_P1"
    else:
        print >>fpairs, "X0_P0 Y0_P0 X0_P1 Y0_P1 FLUX_P0 FLUX_P1"
    
    for g in gals:
        maxN = g.getMaxNuc(masked=False)
        print >>f1,"%6d %.3f %.6f %.6f %.3f %d"%(g.ident, g.z, g.ra if g.ra else np.NaN, 
                                                 g.dec if g.dec else np.NaN,
                                                 g.getPearsonR(masked=False),     
                                                g.num_nuclei(mask=False)),
        for p in np.asarray(g.nuclei):
            print >>f1, "%.5e"%(p.getDist(g.nuclei[maxN].x0,
                                         g.nuclei[maxN].y0,
                                         scale=config.pixelscale)),
        for p in g.nuclei:
            print >>f1, "%.5e"%(p.allflux),
        print >>f1, ""
    
    
    for g in gals[n_nucs>0]:
        print >>f,"%6d %.3f %.6f %.6f %.3f %d"%(g.ident, g.z, g.ra if g.ra else np.NaN, 
                                                g.dec if g.dec else np.NaN,
                                                g.getPearsonR(),
                                               g.num_nuclei(**clean_params)),
        if g in gals[in_pairs]:
            print >>fpairs,"%6d %.3f %.6f %.6f %.3f %d"%(g.ident, g.z, g.ra if g.ra else np.NaN, 
                                                         g.dec if g.dec else np.NaN,
                                                         g.getPearsonR(),
                                                       g.num_nuclei(**clean_params)),
        print >>fpeaks, "%6d %.3f"%(g.ident, g.z),
        
        maxN = g.getMaxNuc(masked=True, **clean_params)
        
        goodp = np.asarray(g.nuclei)[g.good_nuclei(**clean_params)]
        for p in goodp:
            print >>f, "%.5e"%(p.getDist(goodp[maxN].x0,
                                         goodp[maxN].y0, scale=config.pixelscale)),
            print >>fpeaks, "%.3f %.3f"%(p.x0, p.y0),

        for p in goodp:
            print >>f, "%.5e"%(p.allflux),

        print >>f,''
        print >>fpeaks,""
      
        if len(goodp) >= 2:
            order = sorted(range(len(goodp)), key=lambda k: goodp[k].flux, reverse=True)
            p0 = goodp[order[0]]
            p1 = goodp[order[1]]
            print >>fpairs, "%.5f %.5e "%(p1.getDist(p0.x0, p0.y0, scale=g.kpc_p_pix), p1.allflux/p0.allflux),
            if doRaDec:
                coords = g.wcs.all_pix2world([p0.x0, p1.x0], [p0.y0, p1.y0], 0)
                print >>fpairs, "%.6f %.6f %.6f %.6f %.3e %.3e"%(coords[0][0], coords[1][0],
                                                                 coords[0][1], coords[1][1],
                                                                p0.allflux, p1.allflux)
            else:
                print >>fpairs, "%.3f %.3f %.3f %.3f %.3e %.3e"%(p0.x0, p0.y0, p1.x0, p1.y0, 
                                                                 p0.allflux, p1.allflux)


    f.close()
    fpeaks.close()
    f1.close()
    fpairs.close()

    if args.makePlots:
        FigCanvas = FigCanvasPS if args.epsPlot else FigCanvasA
        ending='.eps' if args.epsPlot else '.png'
        my_list = np.arange(len(gals))[(n_nucs >=2)]
        #rand = np.random.permutation(my_list)[:100]
        if not os.path.exists(args.path+'/imgs'):
            os.makedirs(args.path+'/imgs')
        s = makehtml(gals, my_list, FigCanvas, ending, args.path)
        htmlFile = open(args.path+'/imgs/imgs.html', 'w')
        htmlFile.write(s)
        htmlFile.close()
        
    return 0


if __name__=='__main__':
    main()
