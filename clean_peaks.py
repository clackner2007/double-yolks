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

import sys, os, re
import argparse
import numpy as np

import matplotlib.figure as figure
from matplotlib.backends.backend_ps import FigureCanvasPS as FigCanvasPS
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvasA
from matplotlib.patches import Ellipse

from cosm_distTable import CosmoDist
from cosm_dist import ang_dist
import pyfits
from scipy.stats import scoreatpercentile

########################
def addCol(data, newarray, newname):
    """
    add a column to a record array
    """
    newarray=np.asarray(newarray)
    newdtype = np.dtype(data.dtype.descr+[(newname, newarray.dtype)])
    newdata = np.empty(data.shape, dtype=newdtype)
    for field in data.dtype.fields:
        newdata[field] = data[field]
    newdata[newname] = newarray

    return newdata

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

    def __init__(self,ra,dec,z,cosmos_id,gal_mag,
                 ang_dist, ang_dist_delz, file_id=0, delz=0.0):
        self.ra=ra
        self.dec=dec
        self.z=z
        self.cosmos_id = cosmos_id
        self.file_id=file_id
        self.nuclei=[]
        self.ad = ang_dist
        self.ad_delz = ang_dist_delz
        self.kpc_p_pix = (self.ad*1000.0)/(206265.0)*0.03
        self.flux = 10**(-0.4*(gal_mag - 25.959))
        self.mag = gal_mag
        self.delz = delz
        self.goodMask = None


    def add_nuclei(self, n):
        self.nuclei.append(n)

    def num_nuclei(self, mask=True):
        if len(self.nuclei)==0:
            return 0
        elif mask:
            return len(np.asarray(self.nuclei)[self.good_nuclei()])
        else:
            return len(self.nuclei)

    def good_nuclei(self, ecut=0.2, dist_cut=(0.0, 8.0), 
                    flux_cut=0.25,
                    totflux_cut=0.03, #0.03, #0.2, #0.05,
                    cent_dist_cut=10.0,
                    center=(8/0.03*0.5,8/0.03*0.5), redo=False):
        if (self.goodMask is not None) & (~redo):
            return self.goodMask

        #if self.cosmos_id==272688:
        #    import pdb
        #    pdb.set_trace()
            
        center = (8.0/0.03*0.5*self.ad_delz/self.ad,#ang_dist(self.z-self.delz)/self.ad,
                  8.0/0.03*0.5*self.ad_delz/self.ad)#ang_dist(self.z-self.delz)/self.ad)
        distcent = np.array([nn.getDist(center[0], center[1], scale=0.03)
                            for nn in self.nuclei])
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

        if (distcent[self.goodMask] >= 1.0).all():
            self.goodMask = np.zeros(len(self.nuclei), dtype=bool)
        if np.abs(self.getPearsonR(masked=True)) > 0.5:
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



    def getMaxNuc(self,masked=True):
        if masked:
            nucs=np.asarray(self.nuclei)[self.good_nuclei()]
        else:
            nucs=np.asarray(self.nuclei)
        mymax=np.argmax([nn.allflux for nn in nucs])
        return mymax

    def plot(self, FigCanvas, ending, path='../imgs/', outdir='wwwnosinc/'):
        fig=figure.Figure((8,4))
        canv=FigCanvas(fig)

        ax1=fig.add_subplot(121, frameon=False)
        img = pyfits.open(path+'%04d_%.6f_%.6f_acs_I_mosaic_30mas_sci.fits'%\
                              (self.file_id, self.ra, self.dec))[0].data
        size=img.shape
        vmax = 14*np.std(img.flatten())
#scoreatpercentile(img[size[0]/4:3*size[0]/4,#
#                                    size[1]/4:3*size[1]/4].flatten(),
#                              99.9)
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
        fig.savefig(outdir+'peaks_{0}'.format(self.cosmos_id)+ending)


###############################################################
def makehtml(gals, mylist, FigCanvas, ending):

    s  = "<html>\n"
    s += "<body>\n"
    s += "<table>\n"
    for i in mylist:
        gals[i].plot(FigCanvas,ending)
        peaks = sorted(np.asarray(gals[i].nuclei)[gals[i].good_nuclei()], 
                       key=lambda p: p.allflux,
                       reverse=True)
                        
        s += "<tr>\n"
        s2 = "<table>\n"
        s2 += "<tr><td>%s:</td><td>%s</td></tr>\n" % \
            ('COSMOS_ID',repr(gals[i].cosmos_id))
        s2 += "<tr><td>%s:</td><td>%.3f</td></tr>\n" % ('Z',gals[i].z)
        s2 += "<tr><td>%s:</td><td>%.6f</td></tr>\n" % ('RA',gals[i].ra)
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
            format('peaks_{0}{1}'.format(gals[i].cosmos_id,ending))
        s += "</tr>\n"


    s += "</table>\n"
    s += "</body>\n"
    s += "</html>"

    print s



################################################################

def main():

    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#    parser.add_argument("-t", "--temp",
#               action='store_true', default=False, help="temp opt",
#               dest='temp')
    parser.add_argument("file_gal", help="galaxy listing")
    parser.add_argument("file_peaks", help="peak listing")
    parser.add_argument("-e", "--eps",
                 action='store_true', default=False, help='make eps plots',
                 dest='epsPlot')
    parser.add_argument("-d", "--deltaz",
                        help='change in redshift from original',
                        default=0.0, type=float, dest='delz')
    parser.add_argument('-l', '--plot', dest='plotid',
                        type=int)
    parser.add_argument("-p", "--path", default='./',
                        help='path to output',
                        dest='path')
    args = parser.parse_args()

    FigCanvas = FigCanvasPS if args.epsPlot else FigCanvasA
    ending='.eps' if args.epsPlot else '.png'


    gals=[]

    g1=np.loadtxt(args.file_gal)
    p1=np.loadtxt(args.file_peaks)
    g1ra=3#2#3
    g1dec=4#3#4
    g1z=1
    g1id=0
    g1file=6#0#6
    g1imag=5#4#5
    g1dz=2
    g1[:,g1dz] *= args.delz
    g1 = np.append(g1, ang_dist(g1[:,g1z]).reshape(len(g1),1), axis=1)#addCol(g1, CosmoDist().ang_dist(g1[g1z]), 'ang_dist')
    g1 = np.append(g1, ang_dist(g1[:,g1z]-g1[:,g1dz]).reshape(len(g1),1), axis=1)#addCol(g1, CosmoDist().ang_dist(g1[g1z]+g1[g1dz]))
    keyp = zip(g1[:,0], g1[:,4]) #for mocks
#    keyp = zip(g1[:,g1ra], g1[:,g1dec]) #for real
    dict_g = dict(zip(keyp,g1))
    dictp = {}
    

    for e in p1:
#        kp = (e[1], e[2])  #for real
        kp = (e[0], e[2]) #for mocks
        if kp not in dictp:
            dictp[kp] = []
        dictp[kp].append(e)

    for k in keyp:
        gg=dict_g[k]
        gals.append(Galaxy(gg[g1ra], gg[g1dec], gg[g1z],
                           int(gg[g1id]), gg[g1imag], 
                            gg[-2], gg[-1], 0, delz=gg[g1dz]))
        if k in dictp:
            for p in dictp[k]:
                gals[-1].add_nuclei(Nuclei(p[6], p[7], int(p[4]),
                                           p[8], p[5], p[9]))
                gals[-1].file_id = int(p[0])

    gals = np.asarray(gals)
    n_nucs = np.array([g.num_nuclei() for g in gals])

    not_pairs = np.where(n_nucs < 2)[0]
    in_pairs = np.where(n_nucs >= 2)[0]

    if not args.plotid:
        f = open(args.path+'foo2', 'w')
        f1 = open(args.path+'nocutfoo2', 'w')
        fpeaks = open(args.path+'foo2_peaks', 'w')
        
        for g in gals[np.where(np.array([g.num_nuclei(mask=False) for g in gals])>=2)[0]]:
            maxN = g.getMaxNuc(masked=False)
            print >>f1,"%6d %.6f %.6f %.3f %d"%(g.cosmos_id, g.ra, g.dec,
                                               g.getPearsonR(),
                                               g.num_nuclei()),
            for p in np.asarray(g.nuclei):
                print >>f1, "%.5e"%(p.getDist(g.nuclei[maxN].x0,
                                             g.nuclei[maxN].y0,
                                             scale=0.03)),
            for p in g.nuclei:
                print >>f1, "%.5e"%(p.allflux),
            print >>f1, ""
        
        
        for g in gals[in_pairs]:
            print >>f,"%6d %.6f %.6f %.3f %d"%(g.cosmos_id, g.ra, g.dec,
                                               g.getPearsonR(),
                                               g.num_nuclei()),
           
            maxN = g.getMaxNuc(masked=True)
            
            goodp = np.asarray(g.nuclei)[g.good_nuclei()]
            #goodp = sorted(goodp, key=lambda pp:pp.allflux)
            #print >>fpeaks,"%6d "%g.cosmos_id,
            for p in goodp:
                print >>f, "%.5e"%(p.getDist(goodp[maxN].x0,
                                             goodp[maxN].y0,
                                             scale=0.03)),
                print >>fpeaks, "%.3f %.3f"%(p.x0, p.y0),

            for p in goodp:
                print >>f, "%.5e"%(p.allflux),
            print >>f,''
            print >>fpeaks,""
            
        f.close()
        fpeaks.close()
        f1.close()

    if args.plotid:

#        for g in gals[in_pairs]:
#            goodp = np.asarray(g.nuclei)[g.good_nuclei()]
#            flux = np.asarray(sorted([n.allflux for n in goodp]))
#            rflux = flux/flux[-1]
#            tflux = flux/g.flux
#
#            for i in range(len(flux)):
#                print "%d %d "%(g.cosmos_id, len(goodp)),
#                print "%.2f %.2f "%(-2.5*np.log10(g.flux)+25.959, g.z),
#                print "%.3e %.3e "%(tflux[-1], np.sum(tflux)),
#                print "%.3e %.3e"%(rflux[i], tflux[i])

        if False:
            pp=np.where(np.array([g.cosmos_id for g in gals])
                        == args.plotid)[0][0]
            print pp
            print gals[pp].file_id
            gals[pp].plot(FigCanvas,ending, outdir=args.path)

        my_list = np.arange(len(gals))[(n_nucs >=2)]
        print len(my_list), len(gals), len(in_pairs)
        rand = np.random.permutation(my_list)[:100]
        makehtml(gals,
                 rand,
                 FigCanvas,ending)

    #if 800821 in [g.cosmos_id for g in gals[in_pairs]]:
    #    print '800821 is in pairs'
    #print len(in_pairs), len(not_pairs)




    return 0


if __name__=='__main__':
    main()
