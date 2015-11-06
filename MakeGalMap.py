#!/usr/bin/env python
#
# MakeObsPlan.py
#
# Elisa Antolini
# Jeremy Heyl
# UBC Southern Observatory
#
# This script takes the LIGO-Virgo Skymap (P(d|m)) and optionally a
# galaxy-density map (P(m)) and finds the most likely fields to
# observe (P(m|d)).  The fields are assumed to be healpix regions from a
# tesselation with a given value of nside (the value of nside
# depends of the field of view of the telescope).
#
#   P(position|data) = P(position) P(data|position) / P(data)
#
#   P(position) is the galaxy density map ( P(m) )
#   P(data|position) is the skymap from LIGO-Virgo ( P(d|m) )
#   P(data) is constant with position so we neglect it.
#
#
# Usage:
#
#    python3 MakeObsPlan.py _nside_ _SkyMap_ [_GalMap_]
#
#
#    nside = ceil ( sqrt (3/Pi) 60 / s )
#
#    where s is the length of one side of the square field of view in degrees.
#
#
# Questions: heyl@phas.ubc.ca
#


import math as mt
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from sys import argv
import pyfits


def IndexToDeclRa(NSIDE,index):
    
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return -np.degrees(theta-mt.pi/2.),np.degrees(mt.pi*2.-phi)

def DeclRaToIndex(decl,RA,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(360.-RA))



#------------------------------------------------------------------------------
# main
#
def main():
    """
    This is the main routine.
    """
    
    #### Input Parameters #####
    
    FitsGalCat_name  = argv[1]
    nvalues          = int(argv[2])
    
    radians_to_deg = 57.2957795
    Lsun = 3.846e26 # Watt
    MagSun = -26.832 # Bolometric


    #Load the Galaxy Catalog
    
    hdulist = pyfits.open(FitsGalCat_name)
    tabledata = hdulist[1].data
    
    #Get RA and DEC in degrees
    
    Gal_RA  = tabledata.field('RA') * radians_to_deg
    Gal_DEC = tabledata.field('DEC') * radians_to_deg
    K_mag   = tabledata.field('KCORR')
    z = tabledata.field('ZPHOTO')
    dist    = (z*3E05)/72
    
    #Create the pixel map
    
    galpixels_GalMap = np.zeros(hp.nside2npix(nvalues))
    
    pixels = DeclRaToIndex(Gal_DEC,Gal_RA,nvalues)

    #LumK     = Lsun*np.power(10,(-MagSun-K_mag)/2.5)
    #LumK    = np.power(10,(-0.4*(-K_mag-5*np.log10(dist)-6.35)))
    #galpixels_GalMap[pixels] += (np.log10(LumK))/np.log10(np.amax(LumK))

    galpixels_GalMap[pixels[(z>5e-2) & (z<6e-2)]] += 1
    
    GalMap_smoothed = hp.sphtfunc.smoothing(galpixels_GalMap,sigma = 0.01)
    hp.mollview(GalMap_smoothed,coord='C',rot = [0,0.3], title='Histogram equalized Ecliptic', unit='prob', xsize=nvalues)
    hp.graticule()
    plt.savefig('./2MPZ_DensityMap.png')
    plt.show()
    
    
    hp.write_map("./2MPZ_DensityMap.fits", galpixels_GalMap)


#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

