#!/usr/bin/env python
#
# ObsPlan.py
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
#
# Usage: ObsPlan.py [-h] [--gal-map GAL_MAP] [--nvalues NVALUES]
#                   [--cumprob CUMPROB] [--savefigures]
#                   [--no-savefigures]
#                   sky-map nside
#
#
#    nside = ceil ( sqrt (3/Pi) 60 / s )
#
#    where s is the length of one side of the square field of view in degrees.
#
#
# Questions: heyl@phas.ubc.ca
#


from argparse import ArgumentParser

import math as mt
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt


def IndexToDeclRa(NSIDE,index):
    
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return -np.degrees(theta-mt.pi/2.),np.degrees(mt.pi*2.-phi)

def DeclRaToIndex(decl,RA,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(360.-RA))


def PlotMap(Map,NsideMap,MapName):
    hp.mollview(Map,coord='C',rot = [0,0.3], title='Histogram-Equalized Probability Density Map', unit='prob', xsize=NsideMap)
    hp.graticule()
    plt.savefig(MapName)

def isPower(num, base):
    if base == 1 and num != 1: return False
    if base == 1 and num == 1: return True
    if base == 0 and num != 1: return False
    power = int (mt.log (num, base) + 0.5)
    return base ** power == num

def MakeObsPlan(SkyMap_name,nside,SaveFigures,nvalues=None,
                cumprob=None,DensityMap_name=None):
    
    #Check if the nside is a power of two
    val = isPower(nside,2)
    
    if val == False:
        print(" **************** WARNING  **************** ")
        print("The inserted NSIDE is not a power of two")
        y = np.log2(nside)
        exp = int(y)
        
        if (exp + 0.5) < y :
            exp = exp +1
        
        nside = int(np.power(2,exp))
        print("The nearest NSIDE applicable is "+str(nside))
        print(" ****************************************** ")
    
    nside_DensityMap = 0
    
    if DensityMap_name != None :
        #Load the Glaxy Density Map P(m)
        
        Densitymap_Ring       = hp.read_map(DensityMap_name,0)
        nside_DensityMap      = hp.pixelfunc.get_nside(Densitymap_Ring)
        galpixels_DensityMap = np.asarray(Densitymap_Ring)
        
        if SaveFigures :
            PlotMap(galpixels_DensityMap,nside_DensityMap,'./GalaxyDensityMap.png')
    
    
    #Load the Sky Map from LIGO-Virgo ( P(d|m) )
    
    Skymap_Ring  = hp.read_map(SkyMap_name,0)
    nside_SkyMap = hp.pixelfunc.get_nside(Skymap_Ring)
    galpixels_SkyMap = np.asarray(Skymap_Ring)
    
    if SaveFigures:
        PlotMap(galpixels_SkyMap,nside_SkyMap,'./LIGOSkyMap.png')
    
    
    #Resize the Sky Map if necessary
    
    if nside_SkyMap != nside:
        
        galpixels_SkyMap = hp.pixelfunc.ud_grade(galpixels_SkyMap,nside_out = nside, order_in = 'RING', order_out = 'RING')
        
        if SaveFigures:
            PlotMap(galpixels_SkyMap,nside,'./LIGOSkyMapResized.png')
    
    
    #Resize Galaxy Density Map if necessary
    
    if DensityMap_name != None :
        if nside_DensityMap != nside:
            galpixels_DensityMap = hp.pixelfunc.ud_grade(galpixels_DensityMap,nside_out = nside, order_in = 'RING', order_out = 'RING')
            if SaveFigures:
                PlotMap(galpixels_DensityMap,nside,'./GalaxyDensityMapResized.png')
    
    Map_Position_Data = np.zeros(hp.nside2npix(nside))
    
    # Multiply the resulting maps together ->
    # P(position|data) = P(position) P(data|position)
    
    if DensityMap_name != None :
        Map_Position_Data = galpixels_SkyMap * galpixels_DensityMap
    else :
        Map_Position_Data = galpixels_SkyMap
    
    
    # Normalize to 1 the sum of the pixels
    Map_Position_Data/=np.sum(Map_Position_Data)
    
    if SaveFigures:
        PlotMap(Map_Position_Data,nside,'./MapPositionData.png')
    
    
    # Sort the array by the probability
    
    # Sort from the largest to the smallest
    
    healpixno=np.argsort(-Map_Position_Data)
    
    # Output Largest to smallest
    Map_Position_Data=Map_Position_Data[healpixno]
    probsum=np.cumsum(Map_Position_Data)
    dec, ra = IndexToDeclRa(nside,healpixno)
    np.savetxt("SkyMap_OutFile.txt.gz",
               np.transpose([healpixno,ra,dec,
                             Map_Position_Data,probsum]),
               fmt="%16d %10.5f %10.5f %10.5f %10.5f",

               header="Healpix Number         Ra        Dec Probability Cumulative Prob")

    # fLigofile = open('./SkyMap_OutFile.txt', 'w')
    # fLigofile.write("# Healpix Number, Ra, Dec, Probability, Cumulative P "+"\n")

    if nvalues != None: 
        print("\n")
        print("# Most "+str(nvalues)+" probable values :")
        import sys
        ii=range(nvalues)
        np.savetxt(sys.stdout,
                   np.transpose([healpixno[ii],ra[ii],dec[ii],
                                 Map_Position_Data[ii],probsum[ii]]),
                   fmt="%16d %10.5f %10.5f %10.5f %10.5f",
                   header="Healpix Number         Ra        Dec Probability Cumulative Prob")

    if cumprob != None: 
        print("\n")
        print("# Most probable values with cumprob < %g" % cumprob)
        import sys
        ii=(probsum<cumprob)
        np.savetxt(sys.stdout,
                   np.transpose([healpixno[ii],ra[ii],dec[ii],
                                 Map_Position_Data[ii],probsum[ii]]),
                   fmt="%16d %10.5f %10.5f %10.5f %10.5f",
                   header="Healpix Number         Ra        Dec Probability Cumulative Prob")



def _parse_command_line_arguments():
    """
    Parse and return command line arguments
    """
    parser = ArgumentParser(
        description=(
            'Command line generating an observing plan from a LIGO/Virgo probability map (with an optional galaxy map too)'
        ),
    )
    parser.add_argument(
        'sky-map',
        type=str,
        help=(
            'A FITS file containing the LIGO/Virgo probability map in HEALPIX format'
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
        '--gal-map',
        required=False,
        type=str,
        help='A FITS file containing the galaxy density map in HEALPIX format'
    )
    parser.add_argument(
        '--nvalues',
        required=False,
        type=int,
        help='Number of Maximum Probability pixels to be shown'
    )
    parser.add_argument(
        '--cumprob',
        required=False,
        type=float,
        help='Output up to the given cumulative probability'
    )
    parser.add_argument('--savefigures',dest='savefigures',action='store_true')
    parser.add_argument('--no-savefigures',dest='savefigures',action='store_false')
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

    args=_parse_command_line_arguments()
    
    MakeObsPlan(args['sky-map'],args['nside'],args['savefigures'],
                nvalues=args['nvalues'],cumprob=args['cumprob'],
                DensityMap_name=args['gal_map'])

'''    
    #### Input Parameters #####

    DensityMap_name  = argv[1]      # Density Map Name or none
    SkyMap_name      = argv[2]      # Sky Map Name
    nside            = int(argv[3]) # NSIDE of probability Map
    SaveFigures      = argv[4]      # Yes or No
    nvalues          = int(argv[5]) # Number of Maximum Probability pixels to be shown

    MakeObsPlan(SkyMap_name,nside,SaveFigures,nvalues,DensityMap_name)
'''




#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    _main()

