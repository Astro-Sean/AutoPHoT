#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Script to load in default commands to autophot
to allow user to update values for their work
'''

def load():

    import os
    from autophot.packages.call_yaml import yaml_syntax as cs

    #  Get location of script
    filepath ='/'.join(os.path.dirname(os.path.abspath(__file__)).split('/'))

    # Name of default input yaml file - do not change

    default_input = 'default_input.yml'

    #  Load default commands
    default_syntax = cs(os.path.join(filepath+'/databases',default_input),'AutoPhOT_input').load_vars()

    return default_syntax
