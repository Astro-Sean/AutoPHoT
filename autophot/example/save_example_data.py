#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
load example data and save it to users desktop by default
'''

import os
import inspect
from pathlib import Path
from shutil import copyfile

def save_fits_to_desktop(home = os.path.join(str(Path.home()),'Desktop')):

    # locaiton of save_data_data.py file (parent directory of this file location)
    current_filpath = '/'.join(inspect.getfile(inspect.currentframe()).split('/')[0:-1])

    # create folder on dekstop called autophot example
    example_directory_path = os.path.join(home,'autophot_example')

    # create directory on dekstop if not already created
    os.makedirs(example_directory_path, exist_ok=True)

    # copy example.fits to desktop
    copyfile(os.path.join(current_filpath,'example.fits'),
             os.path.join(example_directory_path,'example.fits'))