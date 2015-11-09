#!/usr/bin/env python
#
# GalMap.py
#
# Elisa Antolini
# Jeremy Heyl
# UBC Southern Observatory
#
# This script generates a healpix map from a galaxy catalogue
#
# usage: GalMap.py [-h] [--zcolumn ZCOLUMN] [--savefigures] [--no-savefigures]
#                   galaxy-catalogue nside zmin zmax
#
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

from argparse import ArgumentParser

import math as mt
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pyfits


def IndexToDeclRa(NSIDE,index):
    
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return np.degrees(mt.pi/2.0-theta),np.degrees(phi)

def DeclRaToIndex(decl,RA,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(90.-decl),np.radians(RA))

def isPower(num, base):
    if base == 1 and num != 1: return False
    if base == 1 and num == 1: return True
    if base == 0 and num != 1: return False
    power = int (mt.log (num, base) + 0.5)
    return base ** power == num

def MakeGalMap(FitsGalCat_name,nvalues,z_min,z_max,showMap,zcolumn=None,sigma=0.01):
    #Check if the nside is a power of two
    val = isPower(nvalues,2)
    
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
    # Lsun = 3.846e26 # Watt
    # MagSun = -26.832 # Bolometric
    
    
    #Load the Galaxy Catalog
    
    hdulist = pyfits.open(FitsGalCat_name)
    tabledata = hdulist[1].data
    
    #Get RA and DEC in degrees
    
    Gal_RA  = tabledata.field('RA') * radians_to_deg
    Gal_DEC = tabledata.field('DEC') * radians_to_deg
    # K_mag   = tabledata.field('KCORR')
    if zcolumn==None:
        z = tabledata.field('ZPHOTO')
    else:
        z = tabledata.field(zcolumn)
    # dist    = (z*3E05)/72
    
    #Create the pixel map
    
    galpixels_GalMap = np.zeros(hp.nside2npix(nvalues))
    
    pixels = DeclRaToIndex(Gal_DEC,Gal_RA,nvalues)
    
    #LumK     = Lsun*np.power(10,(-MagSun-K_mag)/2.5)
    #LumK    = np.power(10,(-0.4*(-K_mag-5*np.log10(dist)-6.35)))
    #galpixels_GalMap[pixels] += (np.log10(LumK))/np.log10(np.amax(LumK))
    
    galpixels_GalMap[pixels[(z>z_min) & (z<z_max)]] += 1
    
    GalMap_smoothed = hp.sphtfunc.smoothing(galpixels_GalMap,sigma = sigma)
    
    if showMap:
        hp.mollview(galpixels_GalMap,coord='C',rot = [0,0.3],
                    title='Relative Surface Density of Galaxies: %g < z < %g' % (z_min,z_max), unit='prob', xsize=nvalues)
        hp.graticule()
        plt.savefig(FitsMapCat_name+".png")
        plt.show()
        
        hp.mollview(GalMap_smoothed,coord='C',rot = [0,0.3],
                    title='Relative Surface Density of Galaxies: %g < z < %g' % (z_min,z_max), unit='prob', xsize=nvalues)

        hp.graticule()
        plt.savefig(FitsMapCat_name+"_smoothed.png")
        plt.show()
    
    
    hp.write_map(FitsMapCat_name+".fits", galpixels_GalMap)
    hp.write_map(FitsMapCat_name+"_smoothed.fits", GalMap_smoothed)


def _parse_command_line_arguments():
    """
    Parse and return command line arguments
    """
    parser = ArgumentParser(
        description=(
            'Command-line tool to generate a galaxy map from a FITS catalogue'
        ),
    )
    parser.add_argument(
        'galaxy-catalogue',
        type=str,
        help=(
            'A FITS file containing the galaxy catalogue'
        ),
    )
    parser.add_argument(
        'nside',
        type=int,
        help=(
            'nside for the output map'
            'nside = ceil(sqrt(3/Pi) 60 / s)'
            'where s is the length of one side of the square field of view in degrees.'
            'It will be rounded to the nearest power of two.'
        ),
    )
    parser.add_argument(
        'zmin',
        type=float,
        help='The minimum redshift for a galaxy to appear in the output map'
    )
    parser.add_argument(
        'zmax',
        type=float,
        help='The maximum redshift for a galaxy to appear in the output map'
    )
    parser.add_argument(
        '--zcolumn',
        required=False,
        type=str,
        help='A name of the column in FITS file that contains the redshift (default ZPHOTO)'
    )
    parser.add_argument('--smooth',
                        type=float,
                        help='smoothing scale in radians (default 0.01)',
                        required=False)
    parser.set_defaults(smooth=0.01)

    parser.add_argument('--savefigures',dest='savefigures',action='store_true',
                        help='output the healpix data in a png file')
    parser.add_argument('--no-savefigures',dest='savefigures',action='store_false',
                            help='do not output the healpix data in a png file (default)')
    parser.set_defaults(savefigures=False)


    arguments = vars(parser.parse_args())
    return arguments

    
#------------------------------------------------------------------------------
# main
#
def _main():
    """
    This is the main routine.
    """
    
    '''    
    #### Input Parameters #####
    
    FitsGalCat_name  = argv[1]
    nvalues          = int(argv[2])
    z_min            = float(argv[3])
    z_max            = float(argv[4])
    showMap          = argv[5]
    
    MakeGalMap(FitsGalCat_name,nvalues,z_min,z_max,showMap)
    '''

    args=_parse_command_line_arguments()
    MakeGalMap(args['galaxy-catalogue'],args['nside'],args['zmin'],args['zmax'],
                    args['savefigures'],zcolumn=args['zcolumn'],sigma=args['smooth'])
    
#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    _main()

