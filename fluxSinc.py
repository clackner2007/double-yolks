import numpy
import scipy.special as special

#########################################################
#
#
#########################################################    

def wijCoefficients(n, x, y, rad1, rad2, theta, ellipticity):
    wid, xcen, ycen = n, (n+1)/2, (n+1)/2
    dx, dy = x - xcen, y - ycen
    if n%2:
        dx, dy = dx + 1, dy + 1
        
    ftWij = numpy.zeros([wid, wid], dtype=complex)
    scale     = 1.0 - ellipticity
    for iy in range(wid):
        ky = float(iy - ycen)/wid
        for ix in range(wid):
            kx = float(ix - xcen)/wid
            
            # rotate and rescale
            cosT, sinT = numpy.cos(theta), numpy.sin(theta)
            kxr, kyr =  kx*cosT + ky*sinT, scale*(-kx*sinT + ky*cosT)
            k = numpy.sqrt(kxr**2 + kyr**2)
            
            # compute the airy terms, and apply shift theorem
            if k != 0.0:
                airy1 = rad1*special.j1(2.0*numpy.pi*rad1*k)/k
                airy2 = rad2*special.j1(2.0*numpy.pi*rad2*k)/k
            else:
                airy1, airy2 = numpy.pi*rad1**2, numpy.pi*rad2**2
            airy = airy2 - airy1
            phase = numpy.exp(-1.0j*2.0*numpy.pi*(dx*kxr + dy*kyr))
            
            ftWij[iy,ix] = phase*scale*airy

    ftWijShift = numpy.fft.fftshift(ftWij)
    wijShift   = numpy.fft.ifft2(ftWijShift)
    wij        = numpy.fft.fftshift(wijShift)
    return wij.real



def getFlux(img, x, y, r_in, r_out, theta, ellip):

    # the width the postage stamp
    nx, ny = img.shape
    n = 2*int(r_out + 1.0) + 4
    nhalf = n/2
    ix, iy = int(x), int(y)
    ixlo, ixhi = numpy.clip([ix - nhalf, ix + nhalf], 0, nx - 1)
    iylo, iyhi = numpy.clip([iy - nhalf, iy + nhalf], 0, ny - 1)

    stamp = img[ixlo:ixhi,iylo:iyhi]
    
    xcen = x - float(ixlo)
    ycen = y - float(iylo)
    wij = wijCoefficients(n, xcen, ycen, r_in, r_out, theta, ellip)

    # if we clipped at the edge
    if ixlo == 0:
        xshift = nhalf - ix
        wij = wij[xshift:,:]
    elif ixhi < ix + nhalf:
        xshift = nx - 1 - ix - nhalf
        wij = wij[:xshift,:]

    if iylo == 0:
        yshift = nhalf - iy
        wij = wij[:,yshift:]
    elif iyhi < iy + nhalf:
        yshift = ny - 1 - iy - nhalf
        wij = wij[:,:yshift]
        
    return (stamp*wij).sum()
