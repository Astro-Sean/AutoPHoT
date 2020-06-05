#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 12:10:34 2019

@author: seanbrennan
"""

def HOTterPANTS(file,template,syntax = None):

    import subprocess
    import os
    import sys
    import numpy as np
    from pathlib import Path
    import signal
    import time
    from call_functions import getimage,getheader

    try:

        # Get file extension and template data
        fname_ext = Path(file).suffix

        # Open image and template
        file_image     = getimage(file)
        template_image = getimage(template)

        header = getheader(file)

        # some values set to nan during image alignment - ignore these
        check_values_image = file_image[file_image != 1e-5]

        check_values_template = template_image[template_image != 1e-5]

        # Get filename for saving
        base = os.path.splitext(os.path.basename(file))[0]

        # Location of executable
        exe = syntax['hotpants_exe_loc']

        try:
            saturation = float(header['SATURATE'])
        except:
            saturation = 6e4

# =============================================================================
# Argurments to send to HOTPANTS process - list of tuples
# =============================================================================

        # Arguments to pass to HOTPANTS
        include_args = [
                # Input image
                        ('-inim',   str(file)),
                # Template Image
                        ('-tmplim', str(template)),
                # Output image name
                        ('-outim',  str(file.replace(fname_ext,'_subtraction'+fname_ext))),
                # Image lower values
                        ('-il',     str(np.nanmin(check_values_image))),
                # Template lower values
                        ('-tl',     str(np.nanmin(check_values_template))),
                # Image upper values
                        ('-iu',     str(saturation)),
                # Template upper
                        ('-tu',     str(np.nanmax(check_values_template))),
                # Image gain
                        ('-ig',     str(syntax['gain'])),
                # Normalise to image[i]
                        ('-n',      'i'),
                # Background order fitting
#                        ('-bgo' ,   str(1))
                        ]
        args= [str(exe)]

        for i in include_args:
            args[0] += ' ' + i[0] + ' ' + i[1]

# =============================================================================
# Call subprocess using executable location and option prasers
# =============================================================================

        start = time.time()

        with  open(syntax['write_dir'] + base + '_HOTterPANTS.txt', 'w')  as FNULL:

            pro = subprocess.Popen(args,shell=True, stdout=FNULL, stderr=FNULL)

            # Timeout
            pro.wait(syntax['hotpants_timeout'])

            try:
                # Try to kill process to avoid memory errors / hanging process
                os.killpg(os.getpgid(pro.pid), signal.SIGTERM)
                print('HOTPANTS PID killed')
            except:
                pass

        print('HOTPANTS finished: %ss' % round(time.time() - start) )

        # check if file is there - if so return filepath if not return original filepath
        if os.path.isfile(str(file.replace(fname_ext,'_subtraction'+fname_ext))):

            file_size = os.path.getsize(str(file.replace(fname_ext,'_subtraction'+fname_ext)))

            if file_size == 0:
                print('File was created but nothing written')
                print("> FILE CHECK FAILURE - Return original filepath <" )
                return file
            else:
                print('> Subtraction saved as %s < ' % os.path.splitext(os.path.basename(file.replace(fname_ext,'_subtraction'+fname_ext)))[0])
                return str(file.replace(fname_ext,'_subtraction'+fname_ext))

        else:

            print('File was not created')
            print('> FILE CHECK FAILURE - Return orifinal filepath  <')

            return file

    except Exception as e:

        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno,e)
        try:
                # Try to kill process to avoid memory errors / hanging process
            os.killpg(os.getpgid(pro.pid), signal.SIGTERM)
            print('HOTPANTS PID killed')
        except:
            pass

        return file
