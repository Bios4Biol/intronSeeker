#!/usr/bin/env python3

import textwrap
import os

# Metadata
__author__    = "Lasguignes E., Oudin F., Cabanettes F., Klopp C."
__copyright__ = "Copyright (C) 2019 INRA"
__license__   = "GNU General Public License"
__version__   = "1.0"
__email__     = "support.bioinfo.genotoul(at)inra.fr"
__status__    = "dev"


########################################
########### Version Printing ###########
########################################
def program_version() :
    print("simulation2HTML {version}".format(version=__version__),end="\n")
    print(__copyright__) 

########################################
##### Global program help printing #####
########################################

def program_help() :
    text = '\
\n\
Program: intronSeeker simulation2HTML\n\
Version: {version}\n\n\
\
This tool ...\
 ...\
 ...\n\
 \n'.format(version=__version__)

    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="",
    )
    
    cw = textwrap.TextWrapper(
        width=60,
        initial_indent="\t",
        subsequent_indent="\t\t\t",
        break_long_words=False
    )
    # Program Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]))
    
    # Usage
    print('Usage: simulation2HTML <command> [arguments] [--help] [--version]')
    print('(To know the detailed usage of each sub-commands use \'intronSeeker <command> --help\')',end='\n\n')
    
    # Detail of the commands
    print('Commands: ')
    
    print()
    print('Program:   intronSeeker')
    print('Version:   {version}'.format(version=__version__))
    print('License:   {license}'.format(license=__license__))
    print('Copyright: {copyright}'.format(copyright=__copyright__))
    print('Authors:   {author}'.format(author=__author__))
    print('Support:   {email}\n'.format(email=__email__))


########################################
####### Command help dispatching #######
########################################

########################################
##### cmd help printing #####
########################################
