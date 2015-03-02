#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Flood fill algorithm for detecting pixels above a certain
threshold in contiguous regions

Created on Fri Oct 26 16:28:45 2012

@author: clackner
"""

import numpy as np

def is_in(image, xy):
    """
    check that a point is inside the image
    """
    x,y=xy
    return (0 <= x < image.shape[0]) and ( 0 <= y < image.shape[1])

def floodfill(image, start_xy, comp=lambda x: x > 0,
              diagonal=True):
    """
    return lists of pixels containing point start_xy
    above the limit thresh
    includes pixels at a diagonal
    """
    startx, starty = start_xy
    queue = []
    hot_pix = []

    if not is_in(image, (startx, starty)) or not comp(image[startx,starty]):
        return []

    queue.append((startx,starty))

    for pixel in queue:
        if comp(image[pixel]):
            west= pixel
            east = pixel

            #get row of pixels above threshold
            while (west not in hot_pix) and (is_in(image,west)) \
                and comp(image[west]):
                west = (west[0]-1,west[1])
            while (east not in hot_pix) and (is_in(image,east)) \
                and comp(image[east]):
                east = (east[0]+1,east[1])

            if east==west:
                continue

            for longitude in range(west[0]+1,east[0]):
                hot_pix.append((longitude,pixel[1]))

                north = (longitude, pixel[1]+1)
                if (north not in hot_pix) and (is_in(image,north)) \
                    and comp(image[north]):
                    queue.append(north)

                south = (longitude,pixel[1]-1)
                if (south not in hot_pix) and (is_in(image,south)) \
                    and comp(image[south]):
                    queue.append(south)

            if diagonal:
                # does the corners, to get diagonals:
                nw = (west[0],pixel[1]+1)
                if (nw not in hot_pix) and (is_in(image,nw)) \
                    and comp(image[nw]):
                    queue.append(nw)

                sw = (west[0],pixel[1]-1)
                if (sw not in hot_pix) and (is_in(image,sw)) \
                    and comp(image[sw]):
                    queue.append(sw)

                ne = (east[0],pixel[1]+1)
                if (ne not in hot_pix) and (is_in(image,ne)) \
                    and comp(image[ne]):
                    queue.append(ne)

                se = (east[0],pixel[1]-1)
                if (se not in hot_pix) and (is_in(image,se)) \
                    and comp(image[se]):
                    queue.append(se)

    return hot_pix

def make_pix_array(pix_list):
    """
    return array of pixels, a from a where statement
    """
    if len(pix_list) == 0:
        return (np.array([]),np.array([]))
    unzip = zip(*(pix_list))
    return (np.array(unzip[0]), np.array(unzip[1]))

def make_pix_list(pix_array):
    """
    return list of (x,y) pixels from arrays from where statement
    """
    return zip(pix_array[0],pix_array[1])
    
def get_convexhull(pix_list):
    """
    returns a list of x,y which will draw outline around region
    """
    #convex hull
    from scipy.spatial import ConvexHull
    hull = ConvexHull(np.asarray(pix_list))
    line_segs = []
    for s in hull.simplices:
        line_segs.append([hull.points[s,0], hull.points[s,1]])
    return line_segs
    

def main():

    ##test flood fill

    x,y=np.meshgrid(range(15),range(15))
    x -= 7
    y -= 7
    r = np.sqrt(x*x+y*y)*(1-np.random.rand(15,15)**3)

    ff=floodfill(r,(5,6),lambda x: ((x > 2)and(x < 6)),
                diagonal=False)
    print ff
    plt.ion()
    plt.subplot(211)
    plt.imshow(r, interpolation='none')
    plt.colorbar()
    jj=np.ones_like(r)
    jj[make_pix_array(ff)] *= 2
    plt.subplot(212)
    plt.imshow(jj, interpolation='none')
    plt.savefig('foo')
    return 0

if __name__=='__main__':
    import matplotlib.pyplot as plt
    import numpy as np


    main()