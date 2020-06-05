
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 10:05:34 2018

@author: seanbrennan
"""

from __future__ import absolute_import

def teledata2yml(syntax,flst,filepath = None):

    '''
    Check all files in directory for headers to ensure smooth execution of autophot

    '''

#  Autophot Packages

    from autophot.packages.functions import getheader
    import os
    import yaml
    import sys
    import re

    syntax = syntax
    filepath  = filepath

    if filepath == None:

        filepath = syntax['wdir']+'telescope.yml'
        print(' ')
        print('Instrument database: ' + str(filepath))

    if not os.path.exists(filepath):
        with open(filepath, 'w'): pass

    # load existing data
    with open(filepath, 'r') as stream:
            existing_var = yaml.load(stream, Loader=yaml.FullLoader)

    if existing_var == None:
        existing_var = {}

    # Unknown telescoeps i.e. new entry
    unknown_tele = [ ]

    print('\nNumber of files: %s' % len(flst))

    updated_filter_keys = []

# =============================================================================
#     Check for any unknown telescope
# =============================================================================

    for name in flst:
        try:
            headinfo= getheader(name)

            # Filter Info

            if 'TELESCOP' in headinfo:
                inst_key = 'TELESCOP'
            elif 'OBSERVAT' in headinfo:
                inst_key = 'OBSERVAT'

            else:
                inst_key = 'TELESCOP'
                print('Cannot find TELESCOPE - File: %s'  % name)
                if syntax['ignore_no_telescop']:
                    print('Renaming as UNKNOWN')

                    inst_input = 'UNKNOWN'

                    headinfo[inst_key] = (inst_input,'updated by autophot')
                else:
                    print(headinfo.keys())
                    inst_input = str(input(' Name of telescope? :' + '\n' + '-> '))
                    if inst_input== '':
                        inst_input = 'UNKNOWN'
                        headinfo[inst_key] = (inst_input,'updated by autophot')
                    else:

                        headinfo[inst_key] = (inst_input,'updated by autophot')

            if headinfo[inst_key] == '':
                    tele = 'UNKNOWN'
            else:
                tele = headinfo[inst_key]
            if str(tele).strip() not in unknown_tele:

                unknown_tele.append(tele)
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname1, exc_tb.tb_lineno,e)
            print('*** Cannot find telescope name',e)
            pass

        '''
        Filter information:
            for each detected telscope name, check filter syntax

            if filter is similar to catalog_filter_syntax it isn't asked for

            *** need to add generic filter names, sloan, bessel etc

    '''

    catalog_filter_syntax = ['B','V','U','I','R','g','r','i','z','y','u','v','J','H','K'] # should corrospond to catalog syntax

    if len(unknown_tele) != 0:

        print('%s telescopes - checking for unique header keywords' % str(len(unknown_tele)) )

    for j in unknown_tele:

        doc = {str(j):{}}

        # if telescope name no already in database [telescope.yml]
        if j not in existing_var:

            # Name of telescope for labelling
            define = input('Name of Instrument [ Default: '+str(j)+'] : '+'\n'+' -> ' )
            if define == '':
                define = str(j)
            doc[str(j)].update({'Name':define})

            # add default header filter keyword
            doc[str(j)].update({str('filter_key_0'):'FILTER'})

# =============================================================================
# Plate Scale
# =============================================================================
            plate= {'scale_high':'Upper','scale_low':'Lower'}

            scale_type = {'scale_type':'scale_type'}

            for key,val in scale_type.items():
                try:
                    print('Accepted types: arcminwidth, arcsecperpix')
                    define = str(input(val+' Scale unit type: '))
                except:
                    define = None
                doc[str(j)].update({key:define})
                existing_var.update({j:doc[str(j)]})


            for key,val in plate.items():
                try:
                    define = float(input(val+' scale for FOV [arcmins] : '+'\n'+' -> ' ))
                except:
                    define = None
                doc[str(j)].update({key:define})

                existing_var.update({j:doc[str(j)]})

        else:
            doc[str(j)]= existing_var[j]
            existing_var.update({j:doc[str(j)]})

# =============================================================================
#
# =============================================================================

        unknown_filters = [ ]
        known_filters = []
        GAIN = None

        print(' ')
        print('> Checking: %s' % j)

        for name in flst:


                fname = os.path.basename(name)

                try:
                    if fname.endswith(('.fits','.fit','.fts')):

                        headinfo= getheader(name)

                        gain_keywords = [i for i in list(headinfo.keys()) if 'GAIN' in i]


                        if 'GAIN' not in existing_var[j]:
                            for key in gain_keywords:
                                if key == 'GAIN':
                                    GAIN = 'GAIN'

                            while GAIN == None:
                                print('Relevant keywords found:',gain_keywords)
                                define = input('Corrosponding Gain Keyword from list:'+' -> ' )
                                GAIN = str(define)


                            doc[str(j)].update({'gain':GAIN})


                        try:
                            tele_name = headinfo['TELESCOP']
                        except:
                            tele_name = 'UNKNOWN'
                        if tele_name.strip() == '':
                            tele_name = 'UNKNOWN'


                        if tele_name == j:
                            headinfo= getheader(name)
                            pe_filter_keys =  [s for s in existing_var[j].keys() if "filter_key_" in s]

                            while True:

                                # Look for pre-existing keywords
                                for pe in pe_filter_keys:
                                    # if keyword in in header, skip
                                    try:
                                        if headinfo['FILTER'].lower().replace(' ','') != 'clear':
                                            syntax['filter_key'] = 'FILTER'
                                            break
                                    except:

                                        if existing_var[j][pe] not in list(headinfo.keys()):
                                            continue
                                        else:

                                            # if it is check that it's not clear if it isn't select it as current 'filter_key'
                                            if headinfo[existing_var[j][pe]].lower().replace(' ','') != 'clear':
                                                syntax['filter_key'] = existing_var[j][pe]
                                try:
                                    # make sure it has been defined
                                    syntax['filter_key']
                                    break
                                except:



                                    # if no filter keyword is found ask for a new one
                                    print('No pre-existing filter key found')

                                    old_n = int(re.findall(r"[-+]?\d*\.\d+|\d+", pe)[0])

                                    new_filter_header = pe.replace(str(old_n),str(old_n+1))

                                    # try to help -= look for words with 'fil' in it
                                    filter_keys = [i for i in dict(headinfo) if 'FIL' in i]

                                    print('Relevant filter keywords found:',filter_keys)

                                    print([[str(i)+ ' -> '+ str(headinfo[i])]for i in filter_keys])
                                    print('fname - > ' +str(fname))

                                    define = input('Corrosponding FILTER Keyword :'+' -> ' )

                                    syntax['filter_key'] = str(define)

#                                    update_filter_keys.append([str(headinfo[i])]for i in filter_keys])
#                                    update_filter_keys.append(define)1
#                                    print(update_filter_keys)


                                    doc[str(j)].update({str(new_filter_header):define})
                                    existing_var[j].update({str(new_filter_header):define})



                                    try:
                                        syntax['filter_key'] =  define
                                        break
                                    except:
                                        sys.exit("Can't find filter header")

                            # Go through corrospsonding filter header keywords and tranlste them into standard filter notation

                            if str(headinfo[syntax['filter_key']]).strip() not in doc[str(j)]:

                                if str(headinfo[syntax['filter_key']]).strip() not in unknown_filters:

                                    unknown_filters.append(str(headinfo[syntax['filter_key']]).strip())
                                else:
                                    known_filters.append(str(headinfo[syntax['filter_key']]).strip())
                                for i in list(set(unknown_filters)):

                                    if i.strip() in catalog_filter_syntax:
                                        doc[str(j)].update({str(i.strip()):str(i.strip())})
                                    else:
                                        if i not in catalog_filter_syntax:
                                            if i in updated_filter_keys:
                                                break
                                            else:
                                                filter_default = 'no_filter'
                                                print('[Default: '+str(filter_default)+']')
                                                define = input('Corrosponding filter? : '+str(i)+'\n'+' -> ' )
                                                if define == '':
                                                    define = filter_default
                                                doc[str(j)].update({str(i):define})
                                                updated_filter_keys.append(i)

                except Exception as e:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname1, exc_tb.tb_lineno,e)
                    pass

        with open(filepath,'r') as yamlfile:

            cur_yaml = yaml.safe_load(yamlfile)

            if cur_yaml == None:
                cur_yaml = {}

            cur_yaml.update(doc)

        with open(filepath,'w') as yamlfile:
            yaml.safe_dump(cur_yaml, yamlfile,default_flow_style=False)

    return syntax







