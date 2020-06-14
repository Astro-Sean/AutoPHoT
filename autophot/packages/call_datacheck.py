
#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
 Instrument data

 Script to go through files in direcrtory and check if
 appropiate keywords are available in telescope.yml which is found
 in 'wdir'

'''

def checkteledata(syntax,flst,filepath = None):

    import os,sys
    import yaml
    import re
    import numpy as np
    import logging

    # Get header information - priporeity script
    from autophot.packages.functions import getheader

    try:
        logger = logging.getLogger(__name__)
    except:
        import logging as logger



    #  if specific filepath isn't given, use default
    if filepath == None:

        filepath = os.path.join(syntax['wdir'],'telescope.yml')
        logger.info('User instrument database: %s' % str(filepath))

    # create new telescope.yml file in directory if not already made
    if not os.path.exists(filepath):
        with open(filepath, 'w'):
            logger.info('No telescope.yml found: creating new one')
            pass

    # load telescope.yml as exsiting var
    with open(filepath, 'r') as stream:
        existing_var = yaml.load(stream, Loader=yaml.FullLoader)

    # if it's empty initialise and empty dictionary
    if existing_var == None:
        existing_var = {}

    # telescopes list
    tele_inst_dict = {}

    #  keys to upload
    updated_filter_keys = []

    logger.info('\nNumber of files: %.f' % len(flst))

# =============================================================================
#     Check for any unknown telescope
# =============================================================================

    for name in flst:
        try:
            # Load header for every file
            headinfo= getheader(name)

            '''
            Check for all telescope keywords (whether they are in existing var or not)

            looking for the keywrod that describes the name of the telescope, usually telescop, and assigning
            this keyword to inst key

            if it's not found print 'Cannot find TELESCOP'
                > can ignore this and just label as unknown
                > manually asign 'TELESCOP' with unser input name of telescope



            finaly structure of this first loop with create a entry looking like:
            (nested dictionaries)

            telescope name [e.g ESO-NTT]:

                - instrument key [e.g. INSTRUME]:

                    - instrument name [e.g SOFI]:
                        - [e.g. F.o.V, filter headers etc]

                    - instrument name [e.g EFOSC]:
                        -

                - instrument key [if for some file the instrument name comes under something different]:

                    - instrument name [e.g something under another instruemnt key]:
                        -

                    - instrument name :
                        -
            '''
            if 'TELESCOP' in headinfo:
                tele_key = 'TELESCOP'
            elif 'OBSERVAT' in headinfo:
                tele_key = 'OBSERVAT'
            else:

                #  Can't find name of telescope - ask user for name

                tele_key = 'TELESCOP'

                logger.info('Cannot find TELESCOP (or other) keyword: %s' % name)

                # change syntax to True to allow this
                if syntax['ignore_no_telescop']:
                    logger.info('Renaming as UNKNOWN')
                    tele_input = 'UNKNOWN'
                    headinfo[tele_key] = (tele_input,'updated by autophot')
                else:
                    #  print all the keys
                    logger.info(headinfo.keys())
                    tele_input = str(input(' Name of telescope? :' + '\n' + '-> '))
                    if tele_input== '':
                        tele_input = 'UNKNOWN'
                        headinfo[tele_key] = (tele_input,'updated by autophot')
                    else:
                        headinfo[tele_key] = (tele_input,'updated by autophot')

            # Now look for instrument key .... usually 'INSTRUME' by fits standard

            try:
                if 'INSTRUME' in headinfo:
                    inst_key = 'INSTRUME'

                else:

                    while True:
                        inst_key= str(input('Name of instrument key? [type skip to ignore]: '))

                        if inst_key == 'skip':
                            inst_key = 'INSTRUME'
                            inst_input= str(input('Name of instrument? [default: UNKNOWN]: ') or 'UNKNOWN')
                            headinfo[inst_key] = (inst_input,'added by autophot')
                            break

                        elif inst_key in headinfo:
                            break
                        else:
                            print('%s not found in header - please try again' % inst_key)

            except Exception as e:
                logger.exception(e)


            '''
            We have header keyword that describes the name of the telescope
            now append to list
            '''

            tele = headinfo[tele_key]
            inst = headinfo[inst_key]

            # print('Telescope: %s :: Instrument: %s' % (tele,inst))


            #  add name of telescope (from our keyword) to list of telescopes
            if str(tele).strip() not in list(tele_inst_dict.keys()):
                tele_inst_dict[tele] = {}

            # add instrument key to telescope key in tele_inst_dict
            tele_inst_dict[tele][inst_key] = {}

            # add instrument name to tiinstrument key in tele_inst_dict
            tele_inst_dict[tele][inst_key][inst] = {}

        except Exception as e:
            logger.exception('*** Cannot find telescope name: %s' % e)
            pass



    '''
    Filter information:

    for each telescope go through filter and check we have the
    correct keywords to allow to automatic photometry


    Autophot uses simple filter name syntax:

    if filter is similar to catalog_filter_syntax it isn't asked for
    '''

    catalog_filter_syntax = ['B','V','U','I','R','g','r','i',
                             'z','y','u','v','J','H','K'] # should corrospond to catalog syntax


    logger.info('%d telescope(s) detected - checking header keywords' % len(tele_inst_dict))

    # for each detected catalog
    # Master loop
    for i in list(tele_inst_dict.keys()):
        #  List of insturent keys that this telescope uses
        inst_keys = tele_inst_dict[i]

        if i not in existing_var:
            # if telescope doesn't existing add it
            existing_var[i] = {}

        for name in flst:
            fname = os.path.basename(name)
            if fname.endswith(('.fits','.fit','.fts')):

                headinfo = getheader(name)

                try:
                        tele_name = headinfo['TELESCOP']
                except:
                    logger.debug('No TELESCOP  keyword :: setting to UNKOWN')
                    tele_name = 'UNKNOWN'

                if tele_name.strip() == '':
                    logger.debug('Empty TELESCOP keyword :: setting to UNKOWN')
                    tele_name = 'UNKNOWN'


                if tele_name == i:

                    for inst_key in inst_keys:

                        if inst_key  not in existing_var[i]:
                             #  add new entry for instruement key
                            existing_var[i][inst_key] = {}

                        if headinfo[inst_key] not in existing_var[i][inst_key]:

                            # Name of telescope for labelling or take Telescop as default
                            inst_name = str(input('Header Name of Instrument [default: %s]: ' % headinfo[inst_key] ) or headinfo[inst_key])

                            #  update tele_entry with name of telescope
                            existing_var[i][inst_key][inst_name] = {}

                            label_inst_name = str(input('Shorthand name [default: %s]: ' % headinfo[inst_key] ) or headinfo[inst_key])
                            # Name for labelling
                            existing_var[i][inst_key][inst_name]['Name']  = label_inst_name

                            # update  default filter keyword as FILTER
                            existing_var[i][inst_key][inst_name]['filter_key_0']  = 'FILTER'

                            accepted_scale_types = ['arcminwidth', 'arcsecperpix']

                            logger.info('Accpeted scale units')

                            for a in accepted_scale_types:
                                logger.info("- '%s'" % a)

                            while True:
                                scale_units = (input('Scale units [type skip to ignore]: ') or None)

                                if scale_units =='skip':
                                    break

                                if scale_units == None:
                                    logger.info('Error: no entry given')
                                    continue

                                scale_units = str(scale_units.replace("'",'').replace('"',''))

                                if scale_units not in accepted_scale_types:
                                    logger.info('Error: %s not in %s' % (scale_units,accepted_scale_types))
                                else:
                                    break

                            # plate scale of telescope
                            plate= {'scale_high':'Upper','scale_low':'Lower'}

                            # if unit type is skipped, skip upper and lower scales
                            if scale_units == 'skip':

                                # update with scale type seeting to none/null
                                existing_var[i][inst_key][inst_name]['scale_type'] = None
                                for key,val in plate.items():
                                    existing_var[i][inst_key][inst_name][key] = None
                            else:

                                # update with scale type
                                existing_var[i][inst_key][inst_name]['scale_type'] = scale_units

                                #  else ask user for them and ensure upper limit is greater than lower limit
                                while True:

                                    for key, val in plate.items():
                                        scale_limits = float(input(val+' scale for FOV [units: %s]: ' % scale_units) or np.nan)

                                        if np.isnan(scale_limits):
                                            logger.info('Error: No entry given')
                                            continue

                                        existing_var[i][inst_key][inst_name][key] = scale_limits

                                    if existing_var[i][inst_key][inst_name]['scale_low'] >= existing_var[i][key][inst_name]['scale_high']:
                                        logger.info('Error: [scale_low] >= [scale_low] ... try again ...')
                                    else:
                                        break





        '''
        Now go through and check filter header keywords for all telescopes, whether they are known to
        telescope.yml prior to running this script or not

        Development found that although images come from the same instruments, some keywords can change
        so it is safer to loop over all files
        '''

        unknown_filters = [ ]
        known_filters = []

        logger.info('Checking: %s for unique keywords' % i)


        # for every file

        for name in flst:

            fname = os.path.basename(name)

            try:
                if fname.endswith(('.fits','.fit','.fts')):
                    headinfo = getheader(name)
                    try:
                        tele_name = headinfo['TELESCOP']
                    except:
                        logger.debug('No TELESCOP  keyword :: setting to UNKOWN')
                        tele_name = 'UNKNOWN'

                    if tele_name.strip() == '':
                        logger.debug('Empty TELESCOP keyword :: setting to UNKOWN')
                        tele_name = 'UNKNOWN'

                    if tele_name == i:

                        '''
                        Check for filter keys

                        Files can have multiple types of filter keyword

                        filter keywords are saved in telescope.yml as filter_key_[1..2..3..etc]

                        with the default key being filter_key_0 = 'FILTER'

                        '''

                        # find matching intsrument key:
                        inst_keys = list(existing_var[i].keys())

                        logger.debug('Matching "INSTRUME" keys for %s: %s' % (i,inst_keys))

                        for j in inst_keys:
                            if j in headinfo:
                                k = headinfo[j]
                                break


                        # look for words with gain in it
                        gain_keywords = [i for i in list(headinfo.keys()) if 'GAIN' in i]

                        # if specific gain keword not already entered, use gain as keyword in in header
                        if 'GAIN' not in existing_var[i][j][k]:
                            print('Similar gain keywords found: \n%s' % gain_keywords)
                            while True:
                                gain_key = (input('Instrument Gain key [type skip to ignore]:\n[Telescope: %s :: Inst: %s] ' % (i,k)) or None)

                                if gain_key == None:
                                    logger.info('Error: no entry made')
                                    continue
                                elif gain_key == 'skip':
                                    gain_key = None
                                    break
                                else:
                                    break

                            #  add gain keyword and value for gain
                            existing_var[i][j][k]['GAIN'] = gain_key

                        '''
                        Filter keys
                        '''

                        #  Load existing filter_key_[] in under telescope and instrument
                        pe_filter_keys =  [s for s in existing_var[i][j][k].keys() if "filter_key_" in s]

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

                                    if existing_var[i][j][k][pe] not in list(headinfo.keys()):
                                        #  pre-existing key not not found
                                        continue
                                    else:
                                        # if it is check that it's not clear or empty  - if it isn't select  - it as current 'filter_key'
                                        if headinfo[existing_var[i][j][k][pe]].lower().replace(' ','') != 'clear' and headinfo[existing_var[i][j][k][pe]].lower().replace(' ','') != '':
                                            syntax['filter_key'] = existing_var[i][j][k][pe]
                            try:

                                # check filter key is written if it hasn't  been write a new one
                                if syntax['filter_key'] in headinfo:
                                    break
                                else:
                                    raise Exception
                            except:

                                # if no filter keyword is found ask for a new one
                                logger.info('No pre-existing filter key found')

                                #find lastest filer_key_[] value and add +1 to that
                                old_n = int(re.findall(r"[-+]?\d*\.\d+|\d+", pe)[0])
                                new_filter_header = pe.replace(str(old_n),str(old_n+1))

                                '''
                                Find key that corrosponds to filter name
                                '''

                                # try to help and look for words with 'fil' in it
                                filter_keys = [i for i in dict(headinfo) if 'fil' in i.lower()]

                                print('Relevant filter keywords found:')
                                print('File: %s' % fname)

                                print('*** [Key - Value] ***')
                                for f in filter_keys:
                                    print('%s - %s' %  (f,headinfo[f]))

                                filter_key_new = input('Corrosponding FILTER Keyword\n[Telescope: %s :: Inst: %s]: ' % (i,j))

                                syntax['filter_key'] = str(filter_key_new)

                                # tele_entry[str(i)].update({str(new_filter_header):filter_key_new})

                                existing_var[i][j][k][new_filter_header] = filter_key_new

                                # existing_var[i].update({str(new_filter_header):filter_key_new})

                                try:
                                    #  Double check filter key has been defined
                                    syntax['filter_key'] =  filter_key_new
                                    break
                                except Exception as e:
                                    logger.exception(e)
                                    sys.exit("Can't find filter header")


                        '''
                        Now that we have the correct filter key word - make sure that the value that
                        this gives is in a standardised notation
                        e.g rp -> r
                            r' -> r

                        '''

                        #if entry not already in pre-existing data - avoid redoing
                        if str(headinfo[syntax['filter_key']]).strip() not in existing_var[i][j][k]:


                            '''
                            Add to unknown_filters if not already in unknown_filters
                            this is to only label the key even if it appears in in multiple files
                            '''

                            if str(headinfo[syntax['filter_key']]).strip() not in unknown_filters:
                                unknown_filters.append(str(headinfo[syntax['filter_key']]).strip())
                            else:
                                known_filters.append(str(headinfo[syntax['filter_key']]).strip())

                            # for unknown filter  (uf) names
                            for uf in list(set(unknown_filters)):

                                # if it is already in the standard system - without spaces
                                if uf.strip() in catalog_filter_syntax:
                                    #  update entry with filter name
                                    existing_var[i][j][k][uf.strip()] = uf.strip()

                                elif uf not in catalog_filter_syntax:
                                        #  if key has already been dealt with - skip
                                        if uf in updated_filter_keys:
                                            break
                                        else:

                                            filter_default = 'no_filter'
                                            logger.info('[Default: '+str(filter_default)+']')

                                            filter_name_convert = str(input('Corrosponding filter - %s\n[Telescope: %s :: Inst: %s][default: %s]: '  %(uf,i,k,filter_default)) or filter_default )

                                            existing_var[i][j][k][uf.strip()] = filter_name_convert
                                            updated_filter_keys.append(uf)

            except Exception as e:
                logger.exception(e)
                pass


        '''
        Finally load and update telescope.yml
        '''

        with open(filepath,'r') as yamlfile:

            cur_yaml = yaml.safe_load(yamlfile)

            if cur_yaml == None:
                # if file is blank
                cur_yaml = {}

            cur_yaml.update(existing_var)

        with open(filepath,'w') as yamlfile:
            yaml.safe_dump(cur_yaml, yamlfile,default_flow_style=False)

    return syntax







