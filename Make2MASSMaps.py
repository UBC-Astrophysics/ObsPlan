#!/usr/bin/env python
#
# Make2MASSMaps.py
#
# Elisa Antolini
# Jeremy Heyl
# UBC Southern Observatory
#
# This script converts the 2MASS Extended Source Catalogue to a
# several skymaps of the galaxy density for different redshift ranges.
#
# The redshifts of the sources are calculated assuming that all of the sources
# have a K-short absolute magnitude of -24 in the spirit of
#
# Jarrett, T 2004, PASA, 21, 396
# 
# http://adsabs.harvard.edu/abs/2004PASA...21..396J
#
# This is a work in progress.  We have not masked the Galactic plane but
# we have smoothed the galaxy counts with a Gaussian with sigma of 0.01
# radians, yielding an oversampled map of the same format as the LIGO-Virgo
# healpix skymap.
#
# Questions: heyl@phas.ubc.ca
#
#     Copyright 2015, Elisa Antolini and Jeremy Heyl
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
import math as mt
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt


def IndexToDeclRa(NSIDE,index):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return np.degrees(mt.pi/2.0-theta),np.degrees(phi)

def DeclRaToIndex(decl,RA,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(90.-decl),np.radians(RA))


#------------------------------------------------------------------------------
# main
#

def _main():

    """
    This is the main routine.
    """

    h_nside=512       # The same value a bayesstar.fits
    smoothing=0.01    # Smooth the resulting map to 0.01 radians with a Gaussian
    hubbleconstant=72
    speedoflight=3E5
    

    filenameCat = 'XSC_Completed.tbl.gz'

    RA,DEC,JMAG,HMAG,KMAG = np.loadtxt(filenameCat,skiprows= 50,
                                       usecols = (0,1,3,4,5),
                                       dtype=[('f0',float),
                                              ('f1',float),
                                              ('f2',float),
                                              ('f3',float),
                                              ('f4',float)], unpack = True)
    
    MK_Star = -24.0
    DIST = np.power(10,(KMAG-MK_Star)/5)*1E-05

    REDSHIFT = (DIST*hubbleconstant)/speedoflight
    
    print("Interval of K MAG = ",min(KMAG)," - ",max(KMAG))
    print("Interval of distances [Mpc] = ",min(DIST)," - ",max(DIST))
    print("Interval of redshift = ",min(REDSHIFT)," - ",max(REDSHIFT))
    
    # Select Galaxies by redshift under the assumption of MK=MK_Star

    for rmin in np.arange(4)*0.01:
        rmax=rmin+0.01
        if rmin == 0:
            kmin=0
        else:
            kmin=5*np.log10(rmin/hubbleconstant*speedoflight*1E5)+MK_Star
            
        kmax=5*np.log10(rmax/hubbleconstant*speedoflight*1E5)+MK_Star
        galpixels_Range= np.zeros(hp.nside2npix(h_nside))
        include_me = np.logical_and((REDSHIFT > rmin),
                                    np.logical_and((REDSHIFT<rmax),
                                                   (KMAG != 0.0)))
        ra_Range         = RA[include_me]
        dec_Range        = DEC[include_me]
        pix_num_Range    = (DeclRaToIndex(dec_Range,ra_Range,h_nside))
        galpixels_Range[pix_num_Range]+=1
    
        print("Number of objects with %g < z < %g : %d" % (rmin,rmax,len(ra_Range)))
        map = hp.sphtfunc.smoothing(galpixels_Range,sigma = smoothing)
        hp.write_map("%g-%g_raw.fits" % (rmin,rmax),galpixels_Range)
        hp.write_map("%g-%g.fits" % (rmin,rmax),map)
        hp.mollview(map,coord='C',rot = [0,0.3], title='Relative Surface Density of Galaxies: %0.1f < K < %0.1f (%g < z < %g)' % (kmin,kmax,rmin,rmax), unit='prob',xsize = 2048)
        hp.graticule()
        plt.savefig("%g-%g.png" % (rmin,rmax))
        plt.show()

#------------------------------------------------------------------------------
# Start program execution.
#

if __name__ == '__main__':

    _main()


