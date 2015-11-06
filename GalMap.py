#!/usr/bin/env python
#
# MakeObsPlan.py
#
# Elisa Antolini
# Jeremy Heyl
# UBC Southern Observatory
#
# This script 
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

def isPower(num, base):
    if base == 1 and num != 1: return False
    if base == 1 and num == 1: return True
    if base == 0 and num != 1: return False
    power = int (mt.log (num, base) + 0.5)
    return base ** power == num

def MakeGalMap(FitsGalCat_name,nvalues,z_min,z_max,showMap):
    #Check if the nside is a power of two
    val = isPower(nvalues,2)
    print(isPower(nvalues,2))
    
    if val == False:
        print(" **************** WARNING  **************** ")
        print("The inserted NSIDE is not a power of two")
        y = np.log2(nvalues)
        exp = int(y)
        
        if (exp + 0.5) < y :
            exp = exp +1
        
        nvalues = int(np.power(2,exp))
        print("The nearest NSIDE applicable is "+str(nvalues))
        print(" ****************************************** ")
    
    FitsMapCat_name  = FitsGalCat_name+"_"+str(z_min)+"_"+str(z_max)
    
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
    
    galpixels_GalMap[pixels[(z>z_min) & (z<z_max)]] += 1
    
    GalMap_smoothed = hp.sphtfunc.smoothing(galpixels_GalMap,sigma = 0.01)
    
    if showMap == 'yes':
        hp.mollview(galpixels_GalMap,coord='C',rot = [0,0.3], title='Histogram equalized Ecliptic', unit='prob', xsize=nvalues)
        hp.graticule()
        plt.savefig(FitsMapCat_name+".png")
        plt.show()
        
        hp.mollview(GalMap_smoothed,coord='C',rot = [0,0.3], title='Histogram equalized Ecliptic', unit='prob', xsize=nvalues)
        hp.graticule()
        plt.savefig(FitsMapCat_name+"_smoothed.png")
        plt.show()
    
    
    hp.write_map(FitsMapCat_name+".fits", galpixels_GalMap)
    hp.write_map(FitsMapCat_name+"_smoothed.fits", GalMap_smoothed)


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
    z_min            = float(argv[3])
    z_max            = float(argv[4])
    showMap          = argv[5]
    
    MakeGalMap(FitsGalCat_name,nvalues,z_min,z_max,showMap)
    
#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

