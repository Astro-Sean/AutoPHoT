#!/usr/bin/env python3
# -*- coding: utf-8 -*-




# AutoPHoT Specific packages
from __future__ import absolute_import


print(r"""
    _       _       ___ _  _    _____
   /_\ _  _| |_ ___| _ \ || |__|_   _|
  / _ \ || |  _/ _ \  _/ __ / _ \| |
 /_/ \_\_,_|\__\___/_| |_||_\___/|_|

 ---------------------------------------
         AutoPhOT Beta
         S.J.Brennan et al. 2020
         Please provide feedback/bugs to:
         Email: sean.brennan2@ucdconnect.ie""")


'''
Enter something here

'''


from autophot.packages.recover_output import recover
from autophot.packages.run import run_autophot


# Check user is running Python3 if not exit
import sys
if sys.version_info<(3,0,0):
  from platform import python_version
  sys.stderr.write("\nYou need Python3 or later to run AutoPHoT\n")
  sys.stderr.write("Your Version: %s\n" % python_version())
  exit(1)

# Standard packages
import time
from urllib.request import urlopen
import matplotlib
matplotlib.use('Agg') # Does not display figures to users


# =============================================================================
# Internet required to catalog download, template download
# =============================================================================

def internet_on():
    print('Checking internet connection...')
    try:
        urlopen('https://www.google.com/', timeout=10)
        return True
    except:
        return False

if not internet_on():
    print('Not connected to internet - some packages may not work')
else:

    print('Connected to internet')

# =============================================================================
# Run autophot
# =============================================================================



def run(syntax):


    start = time.time()

    print('> Using directory from input.yml:',syntax['fits_dir'] )

    '''
    Run complete autophot package for complete photometric reduction
    '''

    run_autophot(syntax)

    ''' Go through output filepath and see what has already
    been done and produce a human readbale output
    '''
    recover(syntax)

    print('\nTotal Time Taken: %ss' %  round(time.time() - start))
