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
# Usage:
#
#    python3 ObsPlan.py [_GalMap_] _SkyMap_ _nside_ _PlotFigure_ _Nprob
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
from sys import argv


def IndexToDeclRa(NSIDE,index):
    
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return -np.degrees(theta-mt.pi/2.),np.degrees(mt.pi*2.-phi)

def DeclRaToIndex(decl,RA,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(360.-RA))


def PlotMap(Map,NsideMap,MapName):
    hp.mollview(Map,coord='C',rot = [0,0.3], title='Histogram equalized Ecliptic Galaxy Density Map', unit='prob', xsize=NsideMap)
    hp.graticule()

def isPower(num, base):
    if base == 1 and num != 1: return False
    if base == 1 and num == 1: return True
    if base == 0 and num != 1: return False
    power = int (mt.log (num, base) + 0.5)
    return base ** power == num

def MakeObsPlan(DensityMap_name,SkyMap_name,nside,SaveFigures,nvalues):
    
    #Check if the nside is a power of two
    val = isPower(nside,2)
    print(isPower(nside,2))
    
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
    
    
    
    
    GalMap = 0
    nside_DensityMap = 0
    
    if DensityMap_name != 'none' :
        #Load the Glaxy Density Map P(m)
        
        GalMap = 1
        Densitymap_Ring       = hp.read_map(DensityMap_name,0)
        nside_DensityMap      = hp.pixelfunc.get_nside(Densitymap_Ring)
        galpixels_DensityMap = np.asarray(Densitymap_Ring)
        
        if SaveFigures == 'yes' :
            PlotMap(galpixels_DensityMap,nside_DensityMap,'./GalaxyDensityMap.png')
    
    
    #Load the Sky Map from LIGO-Virgo ( P(d|m) )
    
    Skymap_Ring  = hp.read_map(SkyMap_name,0)
    nside_SkyMap = hp.pixelfunc.get_nside(Skymap_Ring)
    galpixels_SkyMap = np.asarray(Skymap_Ring)
    
    if SaveFigures == 'yes' :
        PlotMap(galpixels_SkyMap,nside_SkyMap,'./LIGOSkyMap.png')
    
    
    #Resize the Sky Map if necessary
    
    if nside_SkyMap != nside:
        
        galpixels_SkyMap = hp.pixelfunc.ud_grade(galpixels_SkyMap,nside_out = nside, order_in = 'RING', order_out = 'RING')
        
        if SaveFigures == 'yes' :
            PlotMap(galpixels_SkyMap,nside,'./LIGOSkyMapResized.png')
    
    
    #Resize Galaxy Density Map if necessary
    
    if GalMap == 1:
        
        if nside_DensityMap != nside:
            
            galpixels_DensityMap = hp.pixelfunc.ud_grade(galpixels_DensityMap,nside_out = nside, order_in = 'RING', order_out = 'RING')
            if SaveFigures == 'yes' :
                PlotMap(galpixels_DensityMap,nside,'./GalaxyDensityMapResized.png')
    
    
    
    
    Map_Position_Data = np.zeros(hp.nside2npix(nside))
    
    # Multiply the resulting maps together -> P(position|data) = P(position) P(data|position)
    
    if GalMap == 1:
        Map_Position_Data = galpixels_SkyMap * galpixels_DensityMap
    else :
        Map_Position_Data = galpixels_SkyMap
    
    
    #Normalize to 1 the pixels
    Map_Position_Data/=np.sum(Map_Position_Data)
    
    if SaveFigures == 'yes' :
        PlotMap(Map_Position_Data,nside,'./MapPositionData.png')
    
    
    #Sort the array by the probability
    
    #Smallest to the largest
    
    healpixno=np.argsort(Map_Position_Data)
    
    #Largest to smallest
    
    sum = 0
    count = 0
    fLigofile = open('./SkyMap_OutFile.txt', 'w')
    fLigofile.write("# Healpix Number, Ra, Dec, Probability, Cumulative P "+"\n")
    
    print("\n")
    print("Most "+str(nvalues)+" probables values :")
    print("# Healpix Number, Ra, Dec, Probability, Cumulative P "+"\n")
    
    for i in healpixno[::-1]:
        # Convert Pixels to RA and DEC in Normalized MAP
        dec, ra = IndexToDeclRa(nside,i)
        sum+=Map_Position_Data[i]
        fLigofile.write(str(i)+" "+str(ra)+" "+str(dec)+" "+str(Map_Position_Data[i])+" "+str(sum)+"\n")
        count += 1
        
        if count <= nvalues :
            print(str(i)+" "+str(ra)+" "+str(dec)+" "+str(Map_Position_Data[i])+" "+str(sum)+"\n")
    
    fLigofile.close()


#------------------------------------------------------------------------------
# main
#
def main():
    """
    This is the main routine.
    """
    
    #### Input Parameters #####
    
    DensityMap_name  = argv[1]      # Density Map Name or none
    SkyMap_name      = argv[2]      # Sky Map Name
    nside            = int(argv[3]) # NSIDE of probability Map
    SaveFigures      = argv[4]      # Yes or No
    nvalues          = int(argv[5]) # Number of Maximum Probability pixels to be shown


    MakeObsPlan(DensityMap_name,SkyMap_name,nside,SaveFigures,nvalues)



#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

