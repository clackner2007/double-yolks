#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################
# PROGRAM NAME
#
# AUTHOR: clackner
#
# DATE: Wed Jul 17 16:48:09 2013
##############################
"""
Created on Wed Jul 17 16:48:09 2013

@author: clackner

program details
"""
import argparse
import numpy as np

import cosm_distTable as cdt

def main():

    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("folder", help='output folder')
    parser.add_argument("inputfile", help='input file')
    parser.add_argument('-l', '--latex', default=False, action='store_true', dest='latex')

    args = parser.parse_args()

    f1 = args.folder + "foo2"
    f2 = args.inputfile #"/data/cosmos_data/photodata/input_photocatxc_mass"

    try:
        gid, z, mag, ra, dec, xid, xmmid = np.loadtxt(f2, unpack=True,
                                                      usecols=(0,2,3,4,5,7,10),
                                                      dtype=[('a',int), ('2',float),
                                                             ('3',float), ('4',float),
                                                             ('5', float), ('6',int), ('7',int)])
    except IndexError:
        gid, z, mag, ra, dec = np.loadtxt(f2, unpack=True,
                                          usecols=(0,2,3,4,5),
                                            dtype=[('a',int), ('2',float),
                                                   ('3',float), ('4',float), ('5', float)])
        xid = np.zeros_like(gid)
        xmmid = np.zeros_like(gid)
    
    with open(f1) as fileH:
        peaks = fileH.readlines()

    backwards = dict(zip(gid, range(len(gid))))
    isx = np.array((xid > 0), dtype=int)
    isxmm = np.array((xmmid>0), dtype=int)

    isdbl = np.zeros(len(gid), dtype=int)
    fluxratio = np.zeros_like(isdbl, dtype=float)
    sep = np.zeros_like(fluxratio)
    tokpc = cdt.CosmoDist(Omega_M=0.25, Omega_L=0.75).ang_dist(z)*1000.0/206265.0
    for p in peaks:
        entries = p.split()
        index = backwards[int(entries[0])]
        isdbl[index] = 1
        numpeaks = int(entries[4])
        #get the separations
        seps = map(float, entries[5:5+numpeaks])
        #fluxes
        flux = map(float, entries[5+numpeaks:])
        #sort by flux, and take the second brightest peak
        #and it's distance to the first
        order = np.argsort(flux)[::-1]
        fluxratio[index] = flux[order[1]]/flux[order[0]]
        sep[index] = tokpc[index]  * seps[order[1]]

    if not args.latex:
    #print "#1-ID 2-RA 3-DEC 4-Z 5-MAGAUTO 6-ISDBL 7-FLUXRATIO 8-SEPKPC 9-ISCHANDRA"
        np.savetxt(args.folder+'pair_info', np.array([gid, ra, dec, z, mag,
                                           isdbl, fluxratio, sep, isx]).T,
                    fmt="%7d %.5f %.5f %.3f %.2f %d %.2f %.2f %d",
                   header="#1-ID 2-RA 3-DEC 4-Z 5-MAGAUTO 6-ISDBL 7-FLUXRATIO 8-SEPKPC 9-ISCHANDRA")
        np.savetxt(args.folder+'all_pairs', np.array([gid[isdbl>0], ra[isdbl>0], dec[isdbl>0], 
                                                      z[isdbl>0], mag[isdbl>0],
                                                      isdbl[isdbl>0], fluxratio[isdbl>0], 
                                                      sep[isdbl>0], isx[isdbl>0]]).T,
                    fmt="%7d %.5f %.5f %.3f %.2f %d %.2f %.2f %d",
                   header="#1-ID 2-RA 3-DEC 4-Z 5-MAGAUTO 6-ISDBL 7-FLUXRATIO 8-SEPKPC 9-ISCHANDRA")

    if args.latex:
        for i in np.where(isdbl>0)[0]:
            xray='--'
            if isx[i]:
                xray = 'Chandra'
                if isxmm[i]:
                    xray += ", XMM"
            elif isxmm[i]:
                xray = "XMM"
            print r"$%.5f$ & $%.5f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%s$ \\"%(ra[i], dec[i], z[i],
                                                                                    mag[i], sep[i]/tokpc[i],
                                                                                    fluxratio[i], xray)


    return 0


if __name__=='__main__':
    main()
