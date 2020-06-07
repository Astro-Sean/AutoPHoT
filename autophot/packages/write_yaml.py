
#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
 Instrument data

 Script to go through files in direcrtory and check if
 appropiate keywords are available in telescope.yml which is found
 in 'wdir'

'''

def teledata2yml(syntax,flst,filepath = None):


    import os,sys
    import yaml
    import re
    import numpy as np

    # Get header information - priporeity script
    from autophot.packages.functions import getheader

    #  if specific filepath isn't given, use default
    if filepath == None:
        filepath = os.path.join(syntax['wdir'],'telescope.yml')
        print('User instrument database: %s' % str(filepath))

    # create new telescope.yml file in directory if not already made
    if not os.path.exists(filepath):
        with open(filepath, 'w'): pass

    # load telescope.yml as exsiting var
    with open(filepath, 'r') as stream:
            existing_var = yaml.load(stream, Loader=yaml.FullLoader)

    # if it's empty initialise and empty dictionary
    if existing_var == None:
        existing_var = {}

    # Unknown telescopes
    tele_list = [ ]


    #  keys to upload
    updated_filter_keys = []

    print('\nNumber of files: %s' % len(flst))


# =============================================================================
#     Check for any unknown telescope
# =============================================================================

    for name in flst:
        try:

            # Load header
            headinfo= getheader(name)


            '''
            Check for all telescope keywords (whether they are in existing var or not)

            looking for the keywrod that describes the name of the telescope, usually telescop, and assigning
            this keyword to inst key

            if it's not found print 'Cannot find TELESCOP'
                > can ignore this and just label as unknown
                > manually asign 'TELESCOP' with unser input name of telescope

            '''
            if 'TELESCOP' in headinfo:
                inst_key = 'TELESCOP'

            elif 'OBSERVAT' in headinfo:
                inst_key = 'OBSERVAT'

            else:
                inst_key = 'TELESCOP'

                print('Cannot find TELESCOP (or other) keyword: %s' % name)

                # change syntax to True to allow this
                if syntax['ignore_no_telescop']:
                    print('Renaming as UNKNOWN')
                    inst_input = 'UNKNOWN'
                    headinfo[inst_key] = (inst_input,'updated by autophot')
                else:
                    #  print all the keys
                    print(headinfo.keys())
                    inst_input = str(input(' Name of telescope? :' + '\n' + '-> '))
                    if inst_input== '':
                        inst_input = 'UNKNOWN'
                        headinfo[inst_key] = (inst_input,'updated by autophot')
                    else:
                        headinfo[inst_key] = (inst_input,'updated by autophot')


            '''
            We have header keyword that describes the name of the telescope
            now append to list
            '''

            if headinfo[inst_key] == '':
                    tele = 'UNKNOWN'
            else:
                tele = headinfo[inst_key]

            #  add name of telescope (from our keyword) to list of telescopes
            if str(tele).strip() not in tele_list:
                tele_list.append(tele)

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname1, exc_tb.tb_lineno,e)
            print('*** Cannot find telescope name: %s' % e)
            pass

    '''
    Filter information:

    for each telescope go through filter and check we have the
    correct keywords to allow to automatic photometry


    Autophot uses simple filter name syntax:

    if filter is similar to catalog_filter_syntax it isn't asked for
    '''

    catalog_filter_syntax = ['B','V','U','I','R','g','r','i','z','y','u','v','J','H','K'] # should corrospond to catalog syntax


    print('%f telescopes detected - checking header keywords' % len(tele_list))

    # for each detected catalog
    for j in tele_list:

        # assign entry to update existing var
        tele_entry = {str(j):{}}

        # if telescope name not already in database [telescope.yml]
        if j not in existing_var:

            # Name of telescope for labelling or take Telescop as default
            define = str(input('Name of Instrument [ Default: %s]: ' % str(j) ) or str(j))

            #  update tele_entry with name of telescope
            tele_entry[str(j)].update({'Name':define})

            # update tele_entry with default filter keyword as FILTER
            tele_entry[str(j)].update({str('filter_key_0'):'FILTER'})

            accepted_scale_types = ['arcminwidth', 'arcsecperpix']

            print('Accpeted scale units')

            for i in accepted_scale_types:
                print("- '%s'" % i)
            while True:
                scale_units = (input('Scale units [type skip to ignore]: ') or None)


                if scale_units =='skip':
                    break

                if scale_units ==None:
                    print('Error: no entry given')
                    continue

                scale_units = str(scale_units.replace("'",'').replace('"',''))

                if scale_units not in accepted_scale_types:
                    print('Error: %s not in %s' % (scale_units,accepted_scale_types))
                else:
                    break

            # plate scale of telescope
            plate= {'scale_high':'Upper','scale_low':'Lower'}

            # if unit type is skipped, skip upper and lower scales
            if scale_units == 'skip':

                # update with scale type
                tele_entry[str(j)].update({'scale_type':None})
                for key,val in plate.items():
                    tele_entry[str(j)].update({key:None})
            else:

                # update with scale type
                tele_entry[str(j)].update({'scale_type':scale_units})

                #  else ask user for them and ensure upper limit is greater than lower limit
                while True:

                    for key,val in plate.items():
                        scale_limits = float(input(val+' scale for FOV [units: %s]: ' % scale_units) or np.nan)

                    if np.isnan(scale_limits):
                        print('Error: No entry given')
                        continue

                    tele_entry[str(j)].update({key:scale_limits})

                    if tele_entry[str(j)]['scale_low'] >= tele_entry[str(j)]['scale_high']:
                        print('Error: [scale_low] >= [scale_low]')
                    else:
                        break


            # update existing var dictioanry with new inputs
            existing_var.update({j:tele_entry[str(j)]})

        else:
            tele_entry[str(j)]= existing_var[j]
            existing_var.update({j:tele_entry[str(j)]})




        '''
        Now go through and check filter header keywords for all telescopes, whether they are known to
        telescope.yml prior to running this script or not

        Development found that although images come from the same instruments, some keywords can change
        so it is safer to loop over all files
        '''

        unknown_filters = [ ]
        known_filters = []

        print('\nChecking: %s for unique keywords' % j)

        for name in flst:
                fname = os.path.basename(name)
                try:
                    if fname.endswith(('.fits','.fit','.fts')):

                        headinfo = getheader(name)
                        # look for words with gain in it
                        gain_keywords = [i for i in list(headinfo.keys()) if 'GAIN' in i]

                        # if specific gain keword not already entered, use gain as keyword in in header
                        if 'GAIN' not in existing_var[j]:
                            print('Similar gain keywords found:',gain_keywords)
                            while True:

                                gain_key = (input('Key matching instrument gain [type skip to ignore]: ') or None)

                                if gain_key == None:
                                    print('Error: no entry made')
                                elif gain_key == 'skip ':
                                    gain_key = 'GAIN'
                                else:
                                    break

                            #  add gain keyword and value for gain
                            tele_entry[str(j)].update({'GAIN':gain_key})



                        try:
                            tele_name = headinfo['TELESCOP']
                        except:
                            tele_name = 'UNKNOWN'
                        if tele_name.strip() == '':
                            tele_name = 'UNKNOWN'




                        if tele_name == j:

                            '''
                            Check for filter keys

                            Files can have multiple types of filter keyword

                            filter keywords are saved in telescope.yml as filter_key_[1..2..3..etc]

                            with the default key being filter_key_0 = 'FILTER'

                            '''

                            headinfo= getheader(name)

                            #  Load existing filter_key_[] in telescope.yml
                            pe_filter_keys =  [s for s in existing_var[j].keys() if "filter_key_" in s]

                            while True:

                                # Look for pre-existing keywords
                                for pe in pe_filter_keys:

                                    try:
                                        # check if 'FILTER' is in headinfo if so break as this is the default
                                        if headinfo['FILTER'].lower().replace(' ','') != 'clear':
                                            syntax['filter_key'] = 'FILTER'
                                            break
                                    except:
                                        '''
                                        Check pre-exising keys i.e ilter_key_ followed by some number
                                        if it's found set use this key and no need to update (as it already exists)
                                        '''

                                        if existing_var[j][pe] not in list(headinfo.keys()):
                                            # not found, try the next
                                            continue
                                        else:
                                            # if it is check that it's not clear if it isn't select it as current 'filter_key'
                                            if headinfo[existing_var[j][pe]].lower().replace(' ','') != 'clear':
                                                syntax['filter_key'] = existing_var[j][pe]

                                try:
                                    # make sure it has been defined if it has break and move on
                                    syntax['filter_key']
                                    break
                                except:
                                    # if it hasn't we need to look for the correct keyword

                                    # if no filter keyword is found ask for a new one
                                    print('No pre-existing filter key found')

                                    #find lastest filer_key_[] value and add +1 to that
                                    old_n = int(re.findall(r"[-+]?\d*\.\d+|\d+", pe)[0])
                                    new_filter_header = pe.replace(str(old_n),str(old_n+1))

                                    # try to help and look for words with 'fil' in it
                                    filter_keys = [i for i in dict(headinfo) if 'FIL' in i]

                                    '''
                                    Find key that corrosponds to filter name
                                    '''

                                    print('Relevant filter keywords found:')
                                    print('File: %s' % fname)
                                    print('[Key - Value]')
                                    for i in filter_keys:
                                        print('%s - %s' %  (headinfo[i],i))

                                    filter_key_new = input('Corrosponding FILTER Keyword :'+' -> ' )

                                    syntax['filter_key'] = str(filter_key_new)

                                    tele_entry[str(j)].update({str(new_filter_header):define})

                                    existing_var[j].update({str(new_filter_header):define})

                                    try:
                                        #  Double check filter key has been defined
                                        syntax['filter_key'] =  define
                                        break
                                    except:
                                        sys.exit("Can't find filter header")


                            '''
                            Now that we have the correct filter key word - make sure that the value that
                            this gives is in a standardised notation
                            e.g rp -> r
                                r' -> r

                            '''
                            #if entry not already in pre-existing data
                            if str(headinfo[syntax['filter_key']]).strip() not in tele_entry[str(j)]:
                                '''
                                Add to unknown_filters if not already in unknown_filters
                                this is to only label the key even if it appears in in multiple files
                                '''
                                if str(headinfo[syntax['filter_key']]).strip() not in unknown_filters:
                                    unknown_filters.append(str(headinfo[syntax['filter_key']]).strip())
                                else:
                                    known_filters.append(str(headinfo[syntax['filter_key']]).strip())

                                # for unknown filter names
                                for i in list(set(unknown_filters)):

                                    # if it is already in the standard system
                                    if i.strip() in catalog_filter_syntax:

                                        #  update entry
                                        tele_entry[str(j)].update({str(i.strip()):str(i.strip())})

                                    elif i not in catalog_filter_syntax:
                                            #  if key has already been dealt with - skip
                                            if i in updated_filter_keys:
                                                break
                                            else:

                                                filter_default = 'no_filter'
                                                print('[Default: '+str(filter_default)+']')

                                                define = str(input('Corrosponding filter - %s- [Default: %s]:'  %(i,filter_default) ) or filter_default )

                                                tele_entry[str(j)].update({str(i):define})
                                                updated_filter_keys.append(i)

                except Exception as e:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname1, exc_tb.tb_lineno,e)
                    pass


        '''
        Finally load and update telescope.yml
        '''

        with open(filepath,'r') as yamlfile:

            cur_yaml = yaml.safe_load(yamlfile)

            if cur_yaml == None:
                # if file is blank
                cur_yaml = {}

            cur_yaml.update(tele_entry)

        with open(filepath,'w') as yamlfile:
            yaml.safe_dump(cur_yaml, yamlfile,default_flow_style=False)

    return syntax







