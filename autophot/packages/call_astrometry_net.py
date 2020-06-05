#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 09:24:59 2018

@author: seanbrennan
"""


'''
- Inputs:
    File: path to file
    syntax: dictionary use for settings


- Output:
    Filename of new, downloaded fits file with updated WCS

'''
def AstrometryNetLOCAL(file, syntax = None):

    import subprocess
    import pathlib
    import os
    import shutil
    import sys
    import numpy as np
    from autophot.packages.functions import getheader
    import signal
    import time

    try:

        #Open file and get header information
        headinfo = getheader(file)

        if syntax == None:
            sys.exit('Astrommetry needs synatx Yaml file')
        # get filename from filepath used to name new WCS fits file contained WCS header with values
        base = os.path.basename(file)
        base = os.path.splitext(base)[0]


        new_file = base + "_WCS.fits"


        # location of executable [/Users/seanbrennan/AutoPhot_Development/AutoPHoT/astrometry.net/bin/solve-field]
        exe = syntax['solve_field_exe_loc']

        # Create input for subprocess:

        # Guess image scale if f.o.v is not known
        if syntax['scale_type'] == '':
            scale = [
                    (str("--guess-scale"))
                    ]
        elif syntax['guess_scale']:

            scale = [
                    (str("--guess-scale"))
                    ]
        else:
            scale = [
                    ("--scale-units=" ,str(syntax['scale_type'])),
                    ("--scale-low="   , str(syntax['scale_low'])),
                    ("--scale-high="  , str(syntax['scale_high'])),
                    ]

        if syntax['target_name'] != None or syntax['target_ra'] != None and syntax['target_dec'] != None:
            try:

                tar_args = [("--ra="       , str(syntax['target_ra'])), # target location on sky, used to narrow down cone search
                            ("--dec="         , str(syntax['target_dec'])),
                            ("--radius="      , str(syntax['search_radius'])) # radius of search around target for matching sources to index deafult 0.5 deg
                            ]
                scale = scale + tar_args
            except:
                pass

        include_args = [

            ('--no-remove-lines'),
            ('--uniformize=' , str(0)),
            	("--overwrite"),
            ("--downsample="  , str(syntax['downsample']) ),  # Downsample image - good for large images
            ("--new-fits="    , str(None)), # Don't download new fits file with updated wcs
            ("--cpulimit="   ,  str(syntax['solve_field_timeout'])), # set time limit on subprocess
            ("--wcs="         , str(new_file)), # filepath of wcs fits header file
            ("--index-xyls="  , str(None)),# don;t need all these files
            ("--axy="         , str(None)),
            ("--scamp="       , str(None)),
            ("--corr="        , str(None)),
            ("--rdl="        ,  str(None)),
            ("--match="      ,  str(None)),
            ("--solved="     ,  str(None)),
            ("--height="      , str(headinfo['NAXIS1'])), #set image height and width
            ("--width="      ,  str(headinfo['NAXIS1'])),
            ("--no-plots"),
            ("--no-verify")
            ]
        # call subprocess using executable location and option prasers

        include_args = include_args + scale

        args= [str(exe) + ' ' + str(file) + ' ' ]

        for i in include_args:
            if isinstance(i,tuple):
                for j in i:
                    args[0] += j
            else:
                args[0] += i
            args[0] += ' '

        start = time.time()

        with  open(syntax['write_dir'] + base + '_astrometry.txt', 'w')  as FNULL:
            pro = subprocess.Popen(args,shell=True, stdout=FNULL, stderr=FNULL,preexec_fn=os.setsid)


            # Timeout
            pro.wait(syntax['solve_field_timeout'])


            try:
                # Try to kill process to avoid memory errors / hanging process
                os.killpg(os.getpgid(pro.pid), signal.SIGTERM)

            except:
                pass

        print('ASTROMETRY finished: %ss' % round(time.time() - start)  )

        # move path to new directory
        new_dir = 'WCS_OUPUT'

        pathlib.Path(new_dir).mkdir(parents=True, exist_ok=True)


        # check if file is there - if so return filepath if not return nan
        if os.path.isfile(str(os.getcwd())+ '/'+ str(new_file)):
            print('> Updated WCS <')

            shutil.move(os.path.join(os.getcwd(), new_file),
                    os.path.join(os.getcwd()+ '/' + new_dir , new_file))
            os.remove(os.getcwd() +'/'+'None')

            return os.getcwd()+ '/' + new_dir + '/' + new_file


        else:
            print("-> FILE CHECK FAILURE - Return NAN <-" )

        return np.nan
    except Exception as e:


        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno,e)
        return np.nan

