#!/usr/bin/env python
#
# MakeObsPlan.py
#
# Elisa Antolini
# Jeremy Heyl
# UBC Southern Observatory
#
# usage: LIGOClient.py [-h] [--gal-map GAL_MAP] [--grace-file GRACE_FILE]
#                      [--nvalues NVALUES] [--cumprob CUMPROB] [--savefigures]
#                      [--no-savefigures] [--textoutput] [--no-textoutput]
#                      graceid nside
#
# Downloads a probability map from the Grace Database and
# then passes it along to ObsPlan to generate an observation plan.
#
from argparse import ArgumentParser
import math as mt
import numpy as np
import healpy as hp
from ligo.gracedb.rest import GraceDbBasic, HTTPError

import ObsPlan



def GetLIGOMap(grace_id,filename=None):
    # Grab the file from the server and write it
    # grace_id = 'T125738'            # identifier for the event
    if filename==None:
        filename = 'bayestar.fits.gz'  # filename of desired skymap
    
    out_filename = grace_id + '_' + filename

    # Instantiate the GraceDB client
    service_url = 'https://gracedb.ligo.org/apibasic/'
    client = GraceDbBasic(service_url)

    # Grab the file from the server and write it
    out_file = open(out_filename, "w")
    r = client.files(grace_id, filename)
    out_file.write(r.read())
    out_file.close()
    
    #out_filename = grace_id + '_' + filename

    return out_filename


def _parse_command_line_arguments():
    """
    Parse and return command line arguments
    """
    parser = ArgumentParser(
        description=(
            'Command line generating an observing plan from a LIGO/Virgo GraceID (with an optional galaxy map too)'
        ),
    )
    parser.add_argument(
        'graceid',
        type=str,
        help=(
            'The Grace-ID for the event'
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
        '--grace-file',
        required=False,
        type=str,
        help='The name of the FITS file containing the probability in HEALPIX format (default of bayestar.fits.gz)'
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

    parser.add_argument('--textoutput',dest='textoutput',action='store_true')
    parser.add_argument('--no-textoutput',dest='textoutput',action='store_false')
    parser.set_defaults(textoutput=False)

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
    DensityMap_name  = argv[1]      # Density Map Name or none
    #SkyMap_name      = argv[2]      # Sky Map Name
    nside            = int(argv[2]) # NSIDE of probability Map
    SaveFigures      = argv[3]      # Yes or No
    nvalues          = int(argv[4]) # Number of Maximum Probability pixels to be shown
    
    ### TO DO
    ######Make this like a routine....
    
    #GCN to listen to the LIGO event and try the test and download the map
    #and generate the observing program
    
    # Grab the file from the server and write it
    grace_id = 'T125738'            # identifier for the event
    filename = 'bayestar.fits.gz'  # filename of desired skymap
    Unzipfilename = grace_id + '_'+'bayestar.fits'
    
    # Prepend with grace_id for output filename
    out_filename = grace_id + '_' + filename
    
    # Instantiate the GraceDB client
    service_url = 'https://gracedb.ligo.org/apibasic/'
    client = GraceDbBasic(service_url)
    
    # Grab the file from the server and write it
    out_file = open(out_filename, "w")
    r = client.files(grace_id, filename)
    out_file.write(r.read())
    out_file.close()

    '''

    args=_parse_command_line_arguments()

    SkyMap_name  = GetLIGOMap(args['graceid'],
                              filename=args['grace_file'])

    ObsPlan.MakeObsPlan(SkyMap_name,args['nside'],args['savefigures'],
                        nvalues=args['nvalues'],cumprob=args['cumprob'],
                        DensityMap_name=args['gal_map'],
                        TextOutput=args['textoutput'])



#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    _main()

