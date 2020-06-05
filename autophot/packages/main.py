
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# from __future__ import absolute_import


"""
Created on Mon Feb  4 15:23:26 2019

@author: seanbrennan
"""

'''

- Main operations of AutoPHoT

Inputs:

    - object_info - yaml file:
        TNS repsonse for target [dictionary]
    - syntax_iter - yaml file:
        yaml file input [dictionary]
    - fpath - str:
        file path [str]

Outputs:
    - output:
        dictionary out output photmetric data:
            fname:
                Original file path
            telescope:
                Telescope used - needed for plotting
            mjd:
                Modified julian date
            zp_[]:
                Photometric zeropoint. dictionary key will be followed by filter name
            zp_[]_err:
                error on photometric zeropoint
            []_inst:
                Instrumental magnitude
            []:
                Calibarted magnitude
            []_err:
                Error on Calibarated magnitude
            SNR:
                Signal to noise ratio of target
            lmag:
                calibrated limiting magnitude
    - fpath:
        filepath to new fits image file

'''


def main(object_info,syntax,fpath):

    # ensure to use copy of original inputed synatc instruction files
    syntax_iter = syntax.copy()

    # Basic Packages
    import sys
    import shutil
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import astroscrappy
    import warnings
    from os.path import dirname
    import pandas as pd
    import pathlib
    import collections
    import lmfit
    from pathlib import Path
    import time
    from photutils import CircularAperture
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.pyplot import Circle

    # Astrpy and photutils
    from astropy.io import fits
    from astropy.stats import sigma_clip
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    import astropy.wcs as wcs
    from astropy.stats import sigma_clipped_stats
    from astropy.time import Time

    # Proprietary modules
    from autophot.packages.functions import gauss_2d,weighted_avg_and_std, getheader,getimage,zeropoint, mag
    from autophot.packages.aperture import ap_phot
    from autophot.packages.check_wcs import updatewcs,removewcs
    from autophot.packages.call_astrometry_net import AstrometryNetLOCAL
    from autophot.packages.call_hotpants import HOTterPANTS
    from autophot.packages.call_yaml import yaml_syntax as cs
    from autophot.packages.uncertain import SNR
    import autophot.packages.find as find
    import autophot.packages.psf as psf
    import autophot.packages.call_catalog as call_catalog
    from autophot.packages.get_template import get_pstars
    from autophot.packages.uncertain import sigma_mag_err
    from autophot.packages.limit import limiting_magnitude_prob


    # Fix the issatty error - https://stackoverflow.com/questions/47069239/consolebuffer-object-has-no-attribute-isatty
    sys.stdout.isatty = lambda: False


    warnings.simplefilter(action='ignore', category=FutureWarning)

    plt.rcParams['xtick.labelsize'] = 7
    plt.rcParams['ytick.labelsize'] = 7
    plt.rcParams['font.size'] = 7
    #plt.rcParams['font.style'] = 'italic'
    #plt.rcParams['figure.autolayout'] = True

    plt.rcParams['figure.subplot.hspace']= 0.25
    plt.rcParams['figure.subplot.wspace']= 0.25

    plt.rcParams['figure.dpi']= 100
    plt.rcParams['axes.titlesize'] = 7
    plt.rcParams['axes.labelsize'] = 7

    plt.rcParams['lines.linewidth'] = 1
    plt.rcParams['savefig.format'] = 'pdf'
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['lines.markersize'] = 3

    plt.rcParams['legend.markerscale'] = 2
    plt.rcParams['legend.fontsize'] = 7

    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 0.5
    plt.rcParams['xtick.minor.width'] = 0.35
    plt.rcParams['xtick.direction'] = 'in'

    #plt.rcParams['ytick.right'] =  False
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 0.5
    plt.rcParams['ytick.minor.width'] = 0.35
    plt.rcParams['ytick.direction'] = 'in'



    def set_size(width,
                 fraction=1,aspect = 1):
        """ Set aesthetic figure dimensions to avoid scaling in latex.

        Parameters
        ----------
        width: float
                Width in pts
        fraction: float
                Fraction of the width which you wish the figure to occupy
        aspect: Float
                height multiple of width

        Returns
        -------
        fig_dim: tuple
                Dimensions of figure in inches
        """
        # Width of figure
        fig_width_pt = width * fraction

        # Convert from pt to inches
        inches_per_pt = 1 / 72.27

        # Golden ratio to set aesthetic figure height
        golden_ratio = (5**.5 - 1) / 2

        # Figure width in inches
        fig_width_in = fig_width_pt * inches_per_pt
        # Figure height in inches
        fig_height_in = fig_width_in * golden_ratio

        fig_dim = (fig_width_in, fig_height_in * aspect)

        return fig_dim

    colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628','#f781bf','#999999',
         '#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd']


    try:
        # reseting print settings
        default_stdout = sys.stdout

        # Preparing output dictionary
        output = collections.OrderedDict({})


# =============================================================================
# Prepare new file
# =============================================================================

        # Change to new working directory set by 'outdir_name' in yaml file
        if syntax_iter['fname']:
            wdir  = syntax_iter['fname']
            work_loc = str(pathlib.Path(dirname(wdir)))

        else:
            wdir = syntax_iter['fits_dir']
            new_dir = '_' + syntax_iter['outdir_name']

            base_dir = os.path.basename(wdir).replace(new_dir,'')
            work_loc = base_dir + new_dir
            # create new output directory is it doesn't exist
            pathlib.Path(dirname(wdir)+'/'+work_loc.replace(' ','')).mkdir(parents = True, exist_ok=True)
            os.chdir(dirname(wdir)+'/'+work_loc)


        '''
        Take in original data, copy it to new location
        New file will have _APT append to fname

        Create names of directories with name=filename
        and account for correct sub directories locations
        '''

        base_wext = os.path.basename(fpath)
        base = os.path.splitext(base_wext)[0]

        # Replace Whitespaces with dash
        base = base.replace(' ','_')

        base = base.replace('_APT','') # If file already has AUTOPHOT in it - for debugging

        base = base.replace('_ERROR','')

        if not syntax_iter['fname']:

            root = dirname(fpath)

            sub_dirs = root.replace(wdir,'').split('/')


            sub_dirs = [i.replace('_APT','').replace(' ','_') for i in sub_dirs]

            cur_dir = dirname(wdir)+'/'+ work_loc

        else:
            sub_dirs = ['']

            cur_dir = dirname(wdir)

        '''
        Move through list of subdirs,
        where the final copied file will be moved
        creating new sub directories with folder name
        extension until it reaches folder

        it will then create a new folder with name of file to which
        every output file is moved to
        '''
        fname_ext = Path(fpath).suffix

        for i in range(len(sub_dirs)):
            if i: # If directory is not blank
                 pathlib.Path(cur_dir +'/' + sub_dirs[i]+'_APT').mkdir(parents = True,exist_ok=True)
                 cur_dir = cur_dir +'/' + sub_dirs[i]+'_APT'

        # create filepath of write directory
        cur_dir = cur_dir + '/' + base
        pathlib.Path(cur_dir).mkdir(parents = True,exist_ok=True)

        # copy new file to new directory
        shutil.copyfile(fpath, (cur_dir+'/'+base + '_APT'+fname_ext).replace(' ','_'))

        # new fpath for working fits file
        fpath = cur_dir + '/' + base + '_APT'+fname_ext

        # base is name [without extension]
        base = os.path.basename(fpath).split('.')[0]

        # write dir is where all files will be saved, pre-iteration
        syntax_iter['write_dir'] = (cur_dir + '/').replace(' ','_')

        # Create text file where all print statments are recorded
        flog = open(cur_dir+'/'+str(base)+'.txt', 'w')

        # Get image and header from function library
        image    = getimage(fpath)
        headinfo = getheader(fpath)

        if object_info == None:
            sys.exit('No Target Info')

        if syntax_iter == None:
            sys.exit("No syntax_iter input file"+"/n"+"*** I dont know what I'm doing! ***")
# =============================================================================
# Bits & Bobs
# =============================================================================
        '''
        Class used to implement verbose action
        Will only write print statements to txt file


        Need to input logging module
        '''

        class Tee(object):
            def __init__(self, *files):
                self.files = files
            def write(self, obj):
                for f in self.files:
                    f.write(obj)
                    f.flush() # If you want the output to be visible immediately
            def flush(self) :
                for f in self.files:
                    f.flush()


#==============================================================================
# Main YAML input and syntax_iter files
#==============================================================================

        # catalog syntax_iter contains header information of selected
        # header given by 'catalog' keyword in yaml file
        filepath ='/'.join(os.path.dirname(os.path.abspath(__file__)).split('/')[0:-1])

        catalog_syntax_yml = 'catalog.yml'
        catalog_syntax_iter = cs(os.path.join(filepath+'/databases',catalog_syntax_yml),syntax['catalog']).load_vars()

        # Telescope syntax_iter file container telescope header information:
        # - Keywords [gain, exposure time, e.t.c]
        # - F.O.V
        # - Filter keywords

        tele_syntax_yml = 'telescope.yml'
        tele_syntax_iter = cs(os.path.join(syntax['wdir'],tele_syntax_yml)).load_vars()


        if not syntax_iter['master_warnings']:
            warnings.filterwarnings("ignore")

# =============================================================================
# Get time of observation
# =============================================================================

        '''
        Geting date of observation from header info
        if 'MJD-OBS' in headinfo, use that
        if not:
            look for 'DATE-AVG' or 'DATE-OBS' and convert to mjd
        if not:
            return np.nan and continue onwards
        '''
        try:

            if 'MJD-OBS' in headinfo:
                j_day = headinfo['MJD-OBS']
            else:
                if 'DATE-OBS' in headinfo:
                    time_obs = headinfo['DATE-OBS']
                if 'DATE-AVG' in headinfo:
                    time_obs = headinfo['DATE-AVG']

                time_obs_iso = Time(time_obs,  format='isot')
                j_day = time_obs_iso.mjd
                headinfo['MJD-OBS'] = (j_day,'Modified Julian Date by AutoPhOT')

        except:
            print('-- WARNING: Modied Julian Date NOT FOUND ' + '/n' +' -- '  )
            j_day = np.nan
            pass

# =============================================================================
# Telescope and instrument parameters§
# =============================================================================

        # Get f.o.v from tele_syntax_iter if needed for wcs and upaate syntax_iter dictionary
        # if f.o.v parametrs are not present, guess_scale in astometry.net
        try:
            telescope = headinfo['TELESCOP']
        except:
            telescope  = 'UNKNOWN'
            headinfo['TELESCOP'] = 'UNKNOWN'

        if telescope  == '':
            telescope  = 'UNKNOWN'
            headinfo['TELESCOP'] = 'UNKNOWN'


        # F.o.V using in astrometry

        syntax_iter['scale_type']  = tele_syntax_iter[telescope]['scale_type']
        syntax_iter['scale_high']  = tele_syntax_iter[telescope]['scale_high']
        syntax_iter['scale_low']   = tele_syntax_iter[telescope]['scale_low']

        # Update outputs with filename, telescope and observation time in mjd

        output.update({'fname':fpath})
        output.update({'telescope':tele_syntax_iter[telescope]['Name']})
        output.update({'mjd':j_day})

# =============================================================================
# Find filter
# =============================================================================
        '''
        Get correct filter keyword with the default being 'FILTER'
        In collaboration with 'write_yaml function'
        will serach for 'FILTER' using filter_key_0 key

        if found sets 'filter_key' in syntax_iter file
        if not:
            file search for filter_key_1 key in telescope_syntax_iter and check
            if result value in headinfo.
        will continue until no more filter_key_[] are in telescope_syntax_iter or
        right keyword is found

        if fails:
            returns filter_key = 'no_filter'



        Was implemented to allow for, although the same telescope/instrument,
        the header name for the filter keyoward may be different
        '''

        filter_header = 'filter_key_0'

        import re

        while True:

            if tele_syntax_iter[telescope][filter_header] not in list(headinfo.keys()):
                old_n = int(re.findall(r"[-+]?\d*\.\d+|\d+", filter_header)[0])
                filter_header = filter_header.replace(str(old_n),str(old_n+1))
            elif tele_syntax_iter[telescope][filter_header].lower() == 'clear':
                continue
            else:
                break

        syntax_iter['filter_key'] = tele_syntax_iter[telescope][filter_header]

        try:

            headinfo[syntax_iter['filter_key']]
        except:
            syntax_iter['filter_key'] = 'no_filter'
            print('ERROR - Filter keywoard == no_filter')

        '''
        Select filter to cross check with selected catalog

        Specific telescope nomenclature found in telescope_syntax_iter.yml

        if no filter is used, keyword is [Clear, no filter, etc],
        script will use a preselected filter ( force_filter keyword )

        if no force_filter is selected,
        using filter presented in input.yml as force_filter

        crappy code

        -- Future implementation to default based on instrument --
        '''

        use_filter =  tele_syntax_iter[headinfo['TELESCOP']][str(headinfo[syntax_iter['filter_key'] ])]

        if use_filter.lower() == 'no_filter':
            if syntax_iter['force_filter'] != 'None':
                use_filter = syntax_iter['force_filter']
            if syntax_iter['force_filter'].lower() != 'clear':
                use_filter = syntax_iter['force_filter']
            else:
                # List of filters, will check spread and return most appropiate filter
                use_filter = syntax_iter['filter_key']



# =============================================================================
# Begin photometric
# =============================================================================

        print('\n' +'File: '+str(base+fname_ext) + ' - PID: '+str(os.getpid()))

        '''
        if verbose is true, all prints will be printed inline as well as
        to text file

        if False, it will only print to text file

        crappy idea - will implement logging module

        '''

        if syntax_iter['verbose']:
            sys.stdout = Tee(sys.stdout, flog)
        else:
            sys.stdout = flog
            print('File: '+str(base) + '['+str(os.getpid())+']'+ '\n')

        start = time. time()

        import datetime

        print('Start Time:',str(datetime.datetime.now()) )
        print('Telescope: %s' % telescope)
        print('Filter: %s'% use_filter)

#        print('MJD: %s'% str(round(j_day)))

# =============================================================================
# File deduction - VERY dependant on
# =============================================================================
        try:
            try:
                from reduce import reduce
                fpath = reduce(fpath,use_filter,syntax_iter)
                image    = fits.getdata(fpath)
                headinfo = getheader(fpath)

            except:
                pass

            if fpath == None:
                    raise Exception


    # =============================================================================
    # Cosmic ray removal usign astroscrappy
    # =============================================================================

            # Removing cosmic rays using 'astroscrappy'
            # documnetation:
            #   - https://buildmedia.readthedocs.org/media/pdf/astroscrappy/latest/astroscrappy.pdf


            if syntax_iter['remove_cmrays']:
                try:
                # if cosmic rays have no already been removed
                        if 'CRAY_RM'  not in headinfo:

                            print('Detecting/removing cosmic ray sources')

                            headinfo = getheader(fpath)

                            # image with cosmic rays
                            image_old = fits.PrimaryHDU(image)

                            # Settings are chosen to minic la cosmic - see documentation
                            cray_remove = astroscrappy.detect_cosmics(image_old.data,
                                                                      inmask = None,
                                                                      satlevel = 2**16,
                                                                      sepmed=False,
                                                                      cleantype='medmask',
                                                                      fsmode='median')
                            # image with cosmic rays removed
                            image = cray_remove[1]

                            if np.nanmedian(image) <=0:
                                headinfo['apt_offset'] = ('T', 'Additon of image offset ')
                                print('Median of image: %.f - applying offset' % np.nanmedian(image))
                                image = image + abs( (np.nanmin(image) ))

                            # Update header and write to new file
                            headinfo['CRAY_RM'] = ('T', 'Comsic rays w/astroscrappy ')
                            fits.writeto(fpath,image,headinfo,overwrite = True,output_verify = 'silentfix+ignore')

                            print('Cosmic rays removed - image updated')

                        else:
                            print('Cosmic sources pre-cleaned')
                except Exception as e:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname, exc_tb.tb_lineno,e)
# =============================================================================
# Gain
# =============================================================================

            try:
                '''
                Look for gain ketword, with preset as 'GAIN',
                same procedure as finding filter keyword.

                if fail:
                    gain = 0

                Doesn't work as intended 12-06-19
                '''

                gain = headinfo[tele_syntax_iter[headinfo['TELESCOP']]['GAIN']]
                print('Gain: %s' % str(round(gain,3)))
            except:
                print('GAIN NOT FOUND - SETTING TO 1')
                gain = 1
                pass

            syntax_iter['gain'] = gain

    # ============================================v=================================
    # Expsoure time
    # =============================================================================

            exp_time = headinfo['EXPTIME']
            '''
            Known issue -> header values written as string

            Below just splits them apart
            '''
            if isinstance(exp_time, str):
               exp_time = exp_time.split('/')
               exp_time = float(exp_time[0])

            syntax_iter['exp_time'] = exp_time


            try:
                syntax_iter['sat_lvl'] = headinfo['SATURATE']
            except:
                syntax_iter['sat_lvl'] = 2**16
            print('Exposure: '+str(exp_time))
    # =============================================================================
    # WCS and target
    # =============================================================================


            '''
            Autophot can be executed with or without target coords

            If valid targett name i.e. 2009ip -  target ra/dec will be taken from tns query and converted to pixel coords
            '''

            if syntax_iter['target_name'] != None:
                '''
                Get target info from [pre-saved] TNS_response
                '''

                target_ra  = object_info['ra']
                target_dec = object_info['dec']

                target_coords = SkyCoord(target_ra , target_dec ,unit = (u.hourangle,u.deg))


                syntax_iter['target_ra'] = target_coords.ra.degree
                syntax_iter['target_dec']= target_coords.dec.degree

            elif syntax_iter['target_ra'] != None and syntax_iter['target_dec'] != None:

                target_coords = SkyCoord(syntax_iter['target_ra'] , syntax_iter['target_dec'] ,unit = (u.deg,u.deg))

                syntax_iter['target_ra'] = target_coords.ra.degree
                syntax_iter['target_dec']= target_coords.dec.degree

            else:
                try:
                    if syntax_iter['use_header_radec']:
                        target_coords = SkyCoord(headinfo['CAT-RA'] , headinfo['CAT-DEC'] ,unit = (u.hourangle,u.deg))

                        syntax_iter['target_ra'] =  target_coords.ra.degree
                        syntax_iter['target_dec']=  target_coords.dec.degree

                except:
                    print( '> NO RA:DEC keywords found - attempting astrometry without <')
    #
    #==============================================================================
    # WCS check
    #==============================================================================

            print('\n --- WCS ---')

            '''
            -- Various instances of when/if to query using astrometry.net --

            if any instance of wcs_keywords are not found in the header infomation it

            if succesful, it will add the UPWCS = T header/value to the header
            hdu of the newly created file
            '''
            '''
            Option in input.yml to check for & remove pre-existing wcs values
            Doesn't work - 12-06-19
            '''

            if syntax_iter['remove_wcs']:
                new_header = removewcs(headinfo,delete_keys = True)
                fits.writeto(fpath,image,new_header,overwrite = True,output_verify = 'silentfix+ignore')
                headinfo = getheader(fpath)

            # search keywords for wcs validation
            wcs_keywords = ['CD1_1','CD1_2','CD2_1','CD2_2',
                            'CRVAL1','CRVAL2','CRPIX1','CRPIX2']

            # if no wcs values found in headinfo, ignore file and exit loop (raise exception)
            if syntax_iter['ignore_no_wcs']:
                if any(i not in headinfo for i in wcs_keywords):
                    print('No wcs found - ignoring_wcs setting == True')
                    raise Exception('ignore files w/o WCS')



            if 'UPWCS'  in headinfo:
                print('Astrometry.net already excuted')

            # Check if headers keywords are present
            elif any(i not in headinfo for i in wcs_keywords):

                # No WCS header found

                print('No WCS values found - attempting to solve field')

                '''
                Call local instance of astrometry.net alogorithm
                astrometry.net documentation:
                https://buildmedia.readthedocs.org/media/pdf/astrometrynet/latest/astrometrynet.pdf

                '''

                # Run local instance of Astrometry.net - returns filepath of wcs file
                astro_check = AstrometryNetLOCAL(fpath, syntax= syntax_iter)

                # Open wcs fits file with wcs values
                new_wcs  = fits.open(astro_check,ignore_missing_end = True)[0].header

                fit_open = getheader(fpath)

                # script used to update per-existing header file with new wcs values
                headinfo_updated = updatewcs(fit_open,new_wcs)

                # update header to show wcs has been checked
                headinfo_updated['UPWCS'] = ('T', 'Cross checked wcs with astrometry.net')

                # Write new header
                fits.writeto(fpath,image,headinfo_updated,overwrite = True,output_verify = 'silentfix+ignore')

                headinfo = getheader(fpath)

                w1 = wcs.WCS(headinfo)

                print('WCS saved to new file')
            else:
                print('WCS Values found')

#==============================================================================
# Load target information from TNS query if needed
#==============================================================================
            # WCS values
            w1 = wcs.WCS(headinfo)

            if syntax_iter['target_name'] == None :


                '''
                If not target coords are given, will use center of field as location on sky
                '''
                if syntax_iter['target_ra'] == None and syntax_iter['target_dec'] == None:

                   '''
                   if no object is given i.e name,ra,dec then take the middle
                   of the screen as the target and get region around it
                   '''

                   # center of image
                   x = image.shape[1]/2
                   y = image.shape[0]/2

                   # translate pixel values to ra,dec
                   n = w1.all_pix2world([x],[y],0 )

                   # get ra,dec in deg/SkyCoord format
                   target_coords = SkyCoord(n[0][0] , n[1][0] ,unit = (u.deg,u.deg))

                   # update syntax_iter file
                   syntax_iter['target_ra'] = target_coords.ra.degree
                   syntax_iter['target_dec']= target_coords.dec.degree

                   syntax_iter['target_x_pix'] = x
                   syntax_iter['target_y_pix'] = y

                   target_x_pix = x
                   target_y_pix = y

                else:

                   '''
                   if no name is given but ra and dec are, use those instead
                   '''

                   target_ra  = syntax_iter['target_ra']
                   target_dec = syntax_iter['target_dec']

                   target_coords = SkyCoord(target_ra , target_dec ,unit = (u.deg,u.deg))

                   target_x_pix, target_y_pix = w1.all_world2pix(target_coords.ra.degree, target_coords.dec.degree, 1)

                   syntax_iter['target_ra'] = target_coords.ra.degree
                   syntax_iter['target_dec']= target_coords.dec.degree

                   syntax_iter['target_x_pix'] = target_x_pix
                   syntax_iter['target_y_pix'] = target_y_pix

            else:
                 if syntax_iter['target_name'] != None:
                    '''
                    Get target info from [pre-saved] TNS_response
                    '''
                    target_ra  = object_info['ra']
                    target_dec = object_info['dec']

                    target_coords = SkyCoord(target_ra , target_dec ,unit = (u.hourangle,u.deg))

                    target_x_pix, target_y_pix = w1.all_world2pix(target_coords.ra.degree, target_coords.dec.degree, 1)

                    syntax_iter['target_ra'] = target_coords.ra.degree
                    syntax_iter['target_dec']= target_coords.dec.degree

                    syntax_iter['target_x_pix'] = target_x_pix
                    syntax_iter['target_y_pix'] = target_y_pix



                 else:
                    try:

                        target_coords = SkyCoord(headinfo['RA'] , headinfo['DEC'] ,unit = (u.hourangle,u.deg))
                        target_x_pix, target_y_pix = w1.all_world2pix(target_coords.ra.degree, target_coords.dec.degree, 1)

                        syntax_iter['target_ra'] =  target_coords.ra.degree
                        syntax_iter['target_dec']=  target_coords.dec.degree

                        syntax_iter['target_x_pix'] = target_x_pix
                        syntax_iter['target_y_pix'] = target_y_pix

                    except:
                        print( '> NO RA:DEC keywords found - attempting astrometry without <')

#==============================================================================
# FWHM - Using total image source detection
#==============================================================================

            print('\n --- Finding FWHM --- ')

            # get approx fwhm, dataframe of sources used and updated syntax_iter
            # returns fwhm from gaussian fit - and dataframe of sources used
            mean_fwhm,df,syntax_iter = find.fwhm(image,syntax_iter)


            print('FWHM (from all sources[%s]): %s '% (str(len (df)),round(mean_fwhm,3)))
            if round(mean_fwhm,3) > 25:
                print('FWHM Error %s - check image quality' % round(mean_fwhm,3) )
                # raise Exception('Reductions failed - skipping')

            # perform approx aperture photometry on sources
            df = find.phot(image,syntax_iter,df,mean_fwhm)

            '''
            ap_corr needs to be moved to another place

            '''
            ap_corr_base = find.ap_correction(image,syntax_iter,df)

#==============================================================================
# Catalog source detecion
#==============================================================================
            if syntax_iter['do_catalog']:
                print('\n --- Catalog searching/matching --- ')

                print('Searching for viable sources')

                while True:

                    # Search for sources in images that have corrospondong magnityide entry in given catalog
                    specificed_catalog = call_catalog.search(image,headinfo,target_coords,syntax_iter,catalog_syntax_iter,use_filter)

                    # Re-aligns catalog sources with source detection and centroid
                    c,syntax_iter = call_catalog.match(image,headinfo,target_coords,syntax_iter,catalog_syntax_iter,use_filter,specificed_catalog,mean_fwhm)

                    # If  UPWCS keyword, image has been already ran through ASTROMETRY, no need to recheck
                    # c = c.drop(c[c['cp_dist'] > syntax_iter['match_dist']].index)
                    # sigma clip distances - avoid mismatches default sigma = 3
                    sigma_dist =  np.array(sigma_clipped_stats(c['cp_dist']))

                    lower_x_bound = syntax_iter['scale']
                    lower_y_bound = syntax_iter['scale']
                    upper_x_bound = syntax_iter['scale']
                    upper_y_bound = syntax_iter['scale']

                    c = c[c.x_pix < image.shape[1] - upper_x_bound]
                    c = c[c.x_pix > lower_x_bound]
                    c = c[c.y_pix < image.shape[0] - upper_y_bound]
                    c = c[c.y_pix > lower_y_bound]

                    print('Average pixel offset: %s '% round(np.nanmedian(list(sigma_dist)),3))




                    '''
                    Attempt to catch a systematic offset - if wcs and catalog values are off by more than 'offset_param'
                    - recheck wcs if not already done so by autophot

                    '''
                    if np.nanmedian(list(sigma_dist)) <= syntax_iter['offset_param']:
                        break

                    if not syntax_iter['allow_wcs_recheck']:
                        break

                    if 'UPWCS' in headinfo:
                        print('UPWCS found - skipping astrometry')
                        break


                    # Sources from catalog match closly with recentered values, i.e we have a good match
                    # If bad WCS values:
                    #    wipe existing Valeus and run thourgh astrometry.net and
                    #    re-run throgh catalog matching

                    print('Removing and Rechecking WCS')

                    # remove wcs values and update -  will delete keys from header file
                    fit_open = fits.open(fpath,ignore_missing_end = True)
                    fit_header_wcs_clean = removewcs(fit_open,delete_keys = True)
                    fit_header_wcs_clean.writeto(fpath , overwrite = True,output_verify = 'silentfix+ignore')
                    fit_open.close()

                    # Run astromety - return filepath of wcs file
                    astro_check = AstrometryNetLOCAL(fpath,syntax = syntax_iter)

                    new_wcs  = fits.open(astro_check,ignore_missing_end = True)
                    fit_open = updatewcs(fit_open,new_wcs)
                    new_wcs.close()

                    fit_open.writeto(os.getcwd()+'/' + base , overwrite = True,output_verify = 'silentfix+ignore')

                    fpath = os.getcwd()+ '/' + base
                    fit_open.close()

                    fit_open = fits.open(fpath,ignore_missing_end = True)
                    headinfo = getheader(fpath)

                    headinfo['UPWCS'] = ('T', 'CROSS CHECKED WITH ASTROMETRY.NET')

                    fit_open.writeto(fpath , overwrite = True,output_verify = 'silentfix+ignore')
                    fit_open.close()


#==============================================================================
# Aperture Photometry on catalog sources
#=============================================================================

                print('\n --- Photometry --- ')

                c_temp_dict = {}

                dist_list = []

                # Go through each matched catalog source and get its' distance to every other source
                for i in c.index:

                    try:
                        dist = np.sqrt((c.x_pix[i]-np.array(c.x_pix))**2 +( c.y_pix[i]-np.array(c.y_pix))**2)
                        dist = dist[np.where(dist!=0)]
                        # add minimum - used to find isolated sources
                        dist_list.append(np.nanmin(dist))

                    except:
                        dist_list.append(np.nan)
                '''
                loop will attempt to perform psf photometry unless unavailable
                further down into the script
                '''

                #Unless selected - perorm psf over aperture
                if syntax_iter['do_ap_phot']:
                    do_ap = True
                else:
                    do_ap = False



# =============================================================================
# Lets Phot!
# =============================================================================

            '''
            Attempt psf photome
            if it selected in input.yml - default: True
            '''

            if syntax_iter['do_psf_phot']:
                try:

                    print('Using PSF Photometry')

                    '''
                    Build model psf from good sources

                    returns:
                        - residual table:
                        - sigma fit from gaussian fitting
                        - heights/amplitudes of sources used to make psf with x and y locations
                        - updated syntax_iter
                    '''

                    r_table,sigma_fit,psf_heights,syntax_iter = psf.build_r_table(image,df,syntax_iter,mean_fwhm)

                    # if it fails, select aperture photometry and exit this attempt
                    if sigma_fit == None:
                        do_ap=True
                        print('> Sigma fit ERROR <')


                    elif np.any(r_table ==  None):
                        do_ap=True
                        print('PSF fitting unavialable')


                    else:

                        # Fwhm to use for model psf and hereafter
                        mean_fwhm = 2 * np.sqrt(2 * np.log(2)) * np.nanmedian(sigma_fit)

                        print('FWHM PSF: %.3f ' % mean_fwhm)
                        if mean_fwhm > 20:
                            print('*** High FWHM %s***' % round(mean_fwhm,1))
                            raise Exception

                        '''
                        Get data of sources used for psf
                        '''

                        psf_stats,_ = psf.do(psf_heights,r_table,syntax_iter,mean_fwhm)

                        '''
                        Get approximate magnitude of psf to check if target is better than psf
                        very dodgey if target is not specific and target is set to center of image

                        '''

                        approx_psf_mag = float(mag(np.nanmean(psf_stats['psf_counts']),0))
                        print('Approx PSF mag',round(approx_psf_mag,3))

                        '''
                        Beginning fitting model psf to table of match sources from catalog
                        returns:
                            - c_psf: updated dataframe with fitted outputs
                            - model_psf - function of psf model used - see psf.py module
                        '''
                        # image.writeto('/Users/seanbrennan/Desktop/fix.fits')

                        hdul = fits.HDUList()
                        hdul.append(fits.PrimaryHDU())


                        hdul.append(fits.ImageHDU(data=image))

                        hdul.writeto('/Users/seanbrennan/Desktop/fix.fits',overwrite = True)

                        c_psf,model_psf = psf.fit(image,c,r_table,syntax_iter,mean_fwhm)

                        '''
                        Get counts/ magnitudes of sources from psf fitting
                        returns:
                            - c_psf: updated dataframe with psf data
                            - syntax_iter: udpated syntax_iter file
                        '''
                        c_psf,syntax_iter = psf.do(c_psf,r_table,syntax_iter,mean_fwhm)

                        # background flux
                        bkg = c_psf['bkg']/exp_time

                        # source flux
                        ap_sum = c_psf.psf_counts/exp_time

                        # Error in flux
                        ap_sum_err = c_psf.psf_counts_err/exp_time

                        # Signal to noise of source
                        SNR_val = np.array(ap_sum/ap_sum_err)



                except Exception as e:

                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname, exc_tb.tb_lineno,e)

                    # if something goes wrong- do aperture  photometry instead
                    do_ap = True
                    pass

            # add fwhm to syntax file
            syntax_iter['fwhm'] = mean_fwhm

            # Perform aperture photoometry if pre-selected or psf fitting wasn't viable
            if syntax_iter['do_ap_phot'] or do_ap == True:

                print('Using Aperture Photometry')

                # list of tuples of pix coordinates of sources
                positions  = list(zip(np.array(c.x_pix),np.array(c.y_pix)))

                '''
                Aperture photometry model:

                    returns:
                        - ap: total sum of counts wtihitn aperture of radius "radius"
                        - bkg: background sum of counts from annulus of inner/outer radius r_in/r_out
                '''
                ap,bkg = ap_phot(positions,
                                 image,
                                 radius =  syntax_iter['ap_size']    * mean_fwhm,
                                 r_in =    syntax_iter['r_in_size']  * mean_fwhm,
                                 r_out =   syntax_iter['r_out_size'] * mean_fwhm)

                # Background flux from annulus
                bkg = bkg/exp_time

                # Source flux
                ap_sum = ap/exp_time

                # SNR from ccd equation
                # needs to be updated
                SNR_val = SNR(ap_sum,bkg,exp_time,0,syntax_iter['ap_size']* mean_fwhm,gain,0)


            # Adding everything to temporary dataframe - no idea why though
            c_temp_dict['dist'] = dist_list
            c_temp_dict['SNR'] = SNR_val
            c_temp_dict['count_rate_bkg'] = bkg
            c_temp_dict['count_rate_star']= ap_sum

            # add to exisitng dataframe [c]
            c_add = pd.DataFrame.from_dict(c_temp_dict)
            c_add = c_add.set_index(c.index)

            c = pd.concat([c,c_add],axis = 1,sort = False)

            # drop if the counts are negative - account for mismatched source or very faint source
            c = c.drop(c[c['count_rate_star'] < 0.0].index)


            '''
            if iso_cat == True:
                remove sources that have a neighbouring source within a user-defined distance
                given by iso_cat_dist.
            '''
            if syntax_iter['iso_cat']:
                c = c.drop(c[c['dist'] < syntax_iter['iso_cat_dist']].index)



    # =============================================================================
    # Plot of image with sources used and target
    # =============================================================================
            fig_source_check = plt.figure(figsize = set_size(540,aspect = 1))

            try:
                from astropy.visualization.mpl_normalize import ImageNormalize
                from astropy.visualization import  ZScaleInterval, SquaredStretch

                x_pix_sources,y_pix_sources = w1.all_world2pix(c.ra.values,c.dec.values,1)

                norm = ImageNormalize( stretch = SquaredStretch())

                vmin,vmax = (ZScaleInterval(nsamples = 1200)).get_limits(image)

                ax = fig_source_check.add_subplot(111)

                ax.imshow(image,
                          vmin = vmin,
                          vmax = vmax,
                          norm = norm,
                          origin = 'lower',
                          cmap = 'Greys')

                ax.scatter(np.array(c.x_pix),np.array(c.y_pix),
                           marker = '+',
                           color = 'red',

                           label = 'Centroiding [%s]' % str(len(c.x_pix)))

                ax.scatter(x_pix_sources,y_pix_sources,
                           marker = 'o',
                           facecolor = 'None',
                           color = 'blue',

                           label = 'Catalog sources [%s]' % str(len(x_pix_sources)))

                if syntax_iter['target_name'] != 'None':
                    tname = syntax_iter['target_name']
                else:
                    tname = 'Center of Field'

                ax.scatter([target_x_pix],[target_y_pix],
                           marker = 'o',
                           facecolor = 'None',
                           color = 'green',
                           label = 'Target: %s' % tname)
                if not do_ap:

                    ax.scatter(psf_heights.x_pix,psf_heights.y_pix,
                               marker = '*',
                               color = 'orange',
                               label = 'PSF Sources')

                ax.set_xlabel('x pixel')
                ax.set_ylabel('y pixel')

                plt.legend(loc = 'lower left')

                if syntax_iter['save_source_plot']:
                    fig_source_check.savefig(cur_dir + '/' +'source_check_'+str(base.split('.')[0])+'.pdf',format = 'pdf')

                plt.close(fig_source_check)


            except Exception as e:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname1, exc_tb.tb_lineno,e)
                plt.close(fig_source_check)




# =============================================================================
# Aperture correction
# =============================================================================

            if do_ap:
                ap_corr = ap_corr_base
            else:
                ap_corr = 0

# =============================================================================
# Find Zeropoint
# =============================================================================


            if syntax_iter['do_zeropoint']:
                print('\n --- Zeropoint --- ' )
                try:
                    '''
                    initilaise dictionaries for zeropoint and errors
                    '''
                    zp = {}
                    zp_err ={}



                    # get magnitude errors on catalog source from SNR
                    zp_mag_err = sigma_mag_err(c['SNR'])
                    dmag = {
                        'U':['U','B'],
                        'B':['B','V'],
                        'V':['V','R'],
                        'R':['V','R'],
                        'I':['R','I'],

                        'u':['u','g'],
                        'g':['g','r'],
                        'r':['g','r'],
                        'i':['r','i'],
                        'z':['i','z']
                        }


                    '''
                    Add zeropoint [zp_[choosen filter]] and zeropoint error [zp_[choosen filter]_err]  to c dataframe:

                    zeropoint:
                        inputs:
                            catalog magnitude:
                            flux of star: given by count_rate_star
                    '''


                    if syntax['apply_ct_zerpoint']:
                        idx = np.where(np.logical_and(~np.isnan(c[dmag[use_filter][0]]) , ~np.isnan(c[dmag[use_filter][1]])))[0]
                        c[dmag[use_filter][0]+'_'+dmag[use_filter][1]] = c[dmag[use_filter][0]] - c[dmag[use_filter][1]]
                        ct_gradient  = -0.0078
                        dmag = c[dmag[use_filter][0]+'_'+dmag[use_filter][1]]
                    else:
                        ct_gradient = None
                        dmag = None

                    c['zp_'+str(use_filter)] = np.array(zeropoint(c[str(use_filter)], c['count_rate_star'],ct_gradient = ct_gradient,dmag = dmag))

                    # error is magnitude error from catalog and instrumental revovery mangitude error from SNR added in qaudrature
                    c['zp_'+str(use_filter)+'_err'] = np.array(np.sqrt( c[str(use_filter)+'_err']**2 + zp_mag_err**2) )

                    # dataframe has been updated with zeropiint calculations - now to sigma clip to get viable zeropiont
                    zpoint = np.asarray(c['zp_'+str(use_filter)])
                    zpoint_err = np.asarray(c['zp_'+str(use_filter)+'_err'])

                    # remove nan values and apply mask
                    nanmask = np.array(~np.isnan(zpoint))
                    zpoint = zpoint[nanmask]
                    zpoint_err = zpoint_err[nanmask]

                    # get zeropoint values within zp_sigma sigma deviation of mean and applt mask

                    zp_mask = np.array(~sigma_clip(zpoint, sigma = syntax_iter['zp_sigma']).mask)
                    zpoint_clip = zpoint[zp_mask]
                    zpoint_err_clip = zpoint_err[zp_mask]

                    # Get instrumental magnitude for catalog sources from autophot photometry
                    zp_inst_mag = mag(c['count_rate_star'],0)[nanmask]
                    # clip according to zp_mask
                    zp_inst_mag_clip =  zp_inst_mag[zp_mask]
                    '''
                    Get weighted average of zeropoints weighted by their magnitude errors
                    '''
                    zpoint_err_clip[zpoint_err_clip ==0] = 1e-30
                    zpoint_err_clip[np.isnan(zpoint_err_clip)] = 1e-30

                    # return value [zp_wa[0]] and error  [zp_wa[1]]
                    zp_wa =  weighted_avg_and_std(np.array(zpoint_clip),np.sqrt(1/zpoint_err_clip))

                    print('> zp_%s: %s +/- %s < ' % (str(use_filter),round(zp_wa[0],3),round(zp_wa[1],3)))
                    # Adding fwhm and zeropint to headerinfo
                    headinfo['fwhm'] = (round(mean_fwhm,3), 'fwhm w/ autophot')
                    headinfo['zp']   = (round(zp_wa[0],3), 'zp w/ autophot')

                    fits.writeto(fpath,image,headinfo,overwrite = True,output_verify = 'silentfix+ignore')

                    # adding to output list
                    output.update(zp)
                    output.update(zp_err)


                except Exception as e:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname1, exc_tb.tb_lineno,e)
                    print('Zeropint not Found')
                    zp_wa = [np.nan,np.nan]
                    pass

                syntax_iter['zp'] =zp_wa[0]
                syntax_iter['zp_err'] =zp_wa[1]

                zp     = {'zp_'+str(use_filter):zp_wa[0]}
                zp_err = {'zp_'+str(use_filter)+'_err':zp_wa[1]}

                headinfo['ZP'] = (zp_wa[0],'ZP by AUTOPHOT')


                # Instrumental magnitude
                c['mag_inst_'+str(use_filter)] = mag(c['count_rate_star'],0)

                # Error in instrumental Magnitude
                c_SNR_err = sigma_mag_err(c.SNR)
                c['mag_inst_'+str(use_filter)+'_err'] = c_SNR_err

                # Observed magnitude
                c['mag_'+str(use_filter)] = mag(c['count_rate_star'],zp_wa[0])

                # Error in observed magnitude
                c['mag_'+str(use_filter)+'_err'] = np.sqrt(c_SNR_err**2 + zp_wa[1]**2)


    # =============================================================================
    #     Plotting Zeropint hisograms w/ clipping
    # =============================================================================
                try:
                    import matplotlib.gridspec as gridspec
                    from scipy.optimize import curve_fit
                    from matplotlib.ticker import MultipleLocator

                    def gauss(x,a,x0,sigma):
                        return a*np.exp(-(x-x0)**2/(2*sigma**2))

                    fig_zeropoint = plt.figure(figsize = set_size(540,aspect = 1))

                    gs = gridspec.GridSpec(2, 2)

                    ax1 = fig_zeropoint.add_subplot(gs[:-1, :-1])
                    ax2 = fig_zeropoint.add_subplot(gs[-1, :-1])
                    ax3 = fig_zeropoint.add_subplot(gs[:, -1])

                    markers, caps, bars = ax1.errorbar(zpoint,zp_inst_mag,xerr = zpoint_err,
                                 label = 'Before clipping',
                                 marker = 'o',
                                 linestyle="None",
                                 color = 'r',
                                 capsize=1,
                                 capthick=1)

                    [bar.set_alpha(0.3) for bar in bars]
                    [cap.set_alpha(0.3) for cap in caps]


                    ax1.set_ylabel('Instrumental magnitude')
                    ax1.set_xlabel('Zeropoint Magnitude')

                    ax1.invert_yaxis()

                    markers, caps, bars = ax2.errorbar(zpoint_clip,zp_inst_mag_clip,
                                 xerr = zpoint_err_clip,
                                 label = 'After clipping [$%s sigma$]' %  str(int(syntax_iter['zp_sigma'])),
                                 marker = 'o',
                                 linestyle="None",
                                 color = 'blue',
                                 capsize=1,
                                 capthick=1)

                    [bar.set_alpha(0.3) for bar in bars]
                    [cap.set_alpha(0.3) for cap in caps]



                    ax2.set_ylabel('Instrumental Magnitude')
                    ax2.set_xlabel('Zeropoint Magnitude')
                    ax2.invert_yaxis()

                    n, bins, patches = ax3.hist(zpoint_clip,
                             bins = 'auto',
                             color = 'green',
                             label = 'Zeropoint Distribution',
                             density = True)

                    mean = zp_wa[0]
                    sigma = 0.1

                    bins_fix = [(bins[i-1] + bins[i]) / 2  for i in range(1,len(bins))]

                    try:
                        popt,pcov = curve_fit(gauss,bins_fix,n, p0=[np.max(n),mean,sigma])
                        r = np.arange(np.min(bins_fix)-0.1,np.max(bins_fix)+0.1,0.01)
                        ax3.plot(r,gauss(r,*popt),'r',label='Gaussian fit')
                        ax3.xlim(np.min(bins_fix)-0.1,np.max(bins_fix)+0.1)
                    except:
                        pass

                    ax3.set_xlabel('Zeropoint Magnitude')
                    ax3.set_ylabel('Probability')

                    for axes in [ax1,ax2,ax3]:
                        axes.legend(loc = 'upper left',fancybox=True,frameon = False)

                    if syntax_iter['save_zp_plot']:
                        fig_zeropoint.savefig(cur_dir + '/' +'zp_'+str(base.split('.')[0])+'.pdf',format = 'pdf')

                    plt.close(fig_zeropoint)

                except Exception as e:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname1, exc_tb.tb_lineno,e)
                    plt.close(fig_zeropoint)

                    pass

# =============================================================================
# Limiting Magnitude
# =============================================================================
            if syntax_iter['do_mag_lim']:
                print ('\n --- Limiting Magnitude --- ')

                from scipy.stats import binned_statistic

                # Create copy of image to avoid anything being written to original image
                image_copy = image.copy()
                try:

                    # Bin size in magnitudes
                    b_size = 0.5

                    lim_err =  sigma_mag_err(syntax_iter['lim_SNR'])

                    print('SNR limit cutoff: %s - Mag error limit: %s ' % (syntax_iter['lim_SNR'],round(lim_err,3)))

                    # Sort values by filter
                    c_mag_lim = c.sort_values(by = [str(use_filter)])

                    # x values: autophot magnitudes
                    x = c_mag_lim[str(use_filter)].values

                    # y values: absolute differnece between catalog magnitude and autophots
                    y =  c_mag_lim['mag_'+str(use_filter)].values - c_mag_lim[str(use_filter)].values

                    # remove nans
                    idx = np.isnan(x)
                    x = x[~idx]
                    y = y[~idx]

                    if len(x)  <= 1 :
                        print('> Magnitude diagram not found <')
                        mag_limit = np.nan
                        test_mag = np.nan
                    else:

                        fig_magnitude = plt.figure(figsize = set_size(540,aspect = 1),facecolor='w', edgecolor='k')
                        ax = fig_magnitude.add_subplot(111)
                        s, edges, _ = binned_statistic(x,y, statistic='median', bins=np.linspace(np.nanmin(x),np.nanmax(x),int((np.nanmax(x)-np.nanmin(x))/b_size)))
                        bin_centers = np.array(edges[:-1]+np.diff(edges)/2)

                        ax.hlines(s,edges[:-1],edges[1:], color="black")

                        ax.scatter(bin_centers, s, c="green",label = 'Binned median')

                        ax.scatter(x,y,color = 'red',marker = 'o',label = r'$Mag_{err}$',alpha = 0.3)
                        ax.axhline(lim_err,label = r'$SNR_{err}$ cuttoff [ '+str(syntax_iter['lim_SNR'])+' ]: ' + str(round(lim_err,3)),linestyle = '--',color = 'red')

                        ax.axhline(-1*lim_err,linestyle = '--',color = 'red')
                        ax.legend(loc = 'lower left')
                        ax.set_ylim(-1,1)
                        ax.set_ylabel(r'$\Delta$ Mag ($Mat_{cat}$ - $M_{recover}$) ',fontsize = 'x-large')
                        ax.set_xlabel(r'$Mag_{cat}$',fontsize = 'x-large')
                        '''
                        Find magnitude where all over where all proceeding magnitude bins are greater than the error SNR cutoff
                        '''
                        for i in range(len(s)):
                            t = s[i]
                            if abs(t)>lim_err:
                                for j in range(i+1,len(s)):
                                    if abs(s[j]) < lim_err:
                                        break
                                else:
                                    ax.axhline(round(t,3),linestyle = ':',label ='Cutoff [ Mag: '+str(round(bin_centers[i],3))+' ]: '+str(round(t,3)), color = 'blue')
                                    break
                        # Set magnitude limit with inall i'th value here
                        mag_limit = bin_centers[i]

                        if syntax_iter['save_mag_lim_plot']:
                            fig_magnitude.savefig(cur_dir + '/' +'mag_lim_'+str(base.split('.')[0])+'.pdf')

                        plt.close(fig_magnitude)


                except Exception as e:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname1, exc_tb.tb_lineno,e)
                    mag_limit = np.nan
                    test_mag = np.nan
                    try:
                        plt.close(fig_magnitude)
                    except:
                        pass
                    pass

                print('Approx. limiting magnitude: %s ' % str(round(mag_limit,3)))

                if target_x_pix < 0 or target_x_pix> image.shape[1] or target_y_pix < 0 or target_y_pix> image.shape[1] :
                    raise Exception ( ' *** EXITING - Target pixel coordinates outside of image [%s , %s] *** ' % (int(target_x_pix), int(target_y_pix)))

                elif do_ap or syntax_iter['do_ap_phot']:
                    test_mag = mag_limit  - syntax_iter['zp']
                    SNR_mock_target = [np.nan]

                    '''
                    if ap_phot selected - sets magnitude limit to nan <- will fix
                    '''
                    mag_limit = np.nan


                image_copy = image.copy()
                if np.isnan(mag_limit):
                    mag_limit = 20

# =============================================================================
# Do photometry on all sources
# =============================================================================
#            if syntax_iter['do_all_phot']:
#                print('\n --- Perform photometry on all sources in field ---')
#
#                _,df_all,_= find.fwhm(image,syntax_iter,sigma_lvl = syntax_iter['do_all_phot_sigma'],fwhm = mean_fwhm)
#
#                ra_all,dec_all = w1.all_pix2world(df_all.x_pix.values,df_all.y_pix.values,0 )
#
#                df_all['ra'] = ra_all
#                df_all['dec'] = dec_all
#
#                if do_ap:
#
#
#                    positions  = list(zip(df_all.x_pix.values,df_all.y_pix.values))
#
#                    target,target_bkg = ap_phot(positions,
#                                                 image,
#                                                 radius = syntax_iter['ap_size']    * mean_fwhm,
#                                                 r_in   = syntax_iter['r_in_size']  * mean_fwhm,
#                                                 r_out  = syntax_iter['r_out_size'] * mean_fwhm)
#
#                    source_flux = (target/exp_time)
#                    source_bkg = (target_bkg/exp_time)
#
#                    SNR_sources = SNR(source_flux,source_bkg,exp_time,0,syntax_iter['ap_size']* mean_fwhm,gain,0)
#                    df_all['snr'] = SNR_sources
#                    df_all['mag'] = mag(source_flux,zp_wa[0] ) + ap_corr
#
#                    mag_err = sigma_mag_err(SNR_sources)
#                    df_all['mag_err'] =  np.sqrt(mag_err**2 + zp_wa[1]**2)
#
#
#
#                else:
#                    positions = df_all[['x_pix','y_pix']]
#
#                    psf_sources,_ = psf.fit(image,
#                                            positions,
#                                            r_table,
#                                            syntax_iter,
#                                            mean_fwhm)
#
#                    psf_sources_phot,_ = psf.do(psf_sources,
#                                                r_table,
#                                                syntax_iter,
#                                                mean_fwhm)
#
#                    sources_flux = np.array(psf_sources_phot.psf_counts/exp_time)
#
#                    sources_err = np.array(psf_sources_phot.psf_counts_err/exp_time)
#
#                    SNR_sources = np.array(sources_flux/sources_err)
#
#                    mag_err = sigma_mag_err(SNR_sources)
#
#                    df_all['snr'] = SNR_sources
#                    df_all['flux_inst'] = psf_sources['H_psf']
#                    df_all['mag'] = mag(sources_flux,zp_wa[0] ) + ap_corr
#
#                    mag_err = sigma_mag_err(SNR_sources)
#                    df_all['mag_err'] =  np.sqrt(mag_err**2 + zp_wa[1]**2)
#
#                try:
#                    df_all.to_csv(syntax_iter['write_dir']+str(base.split('.')[0])+'.cat',index = False)
#                    print('Photometry of all sources saved as: %s' % base+'.cat')
#                except Exception as e:
#
#
#                    exc_type, exc_obj, exc_tb = sys.exc_info()
#                    fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
#                    print(exc_type, fname1, exc_tb.tb_lineno,e)
#                    sys.stdout = default_stdout
#                    print('Warning - Table could not be saved to csv')
#                    flog.close()
#                    pass

# =============================================================================
# Get Template
# =============================================================================
            subtraction_ready = False
            if syntax_iter['do_subtraction']:

                print('\n --- Doing image subtraction  --- ')
                try:
                    # Pan_starrs template images use some old WCS keywrods, astropy can handle them just for cleaniness
                    warnings.filterwarnings('ignore')

                    ra = target_coords.ra.degree
                    dec = target_coords.dec.degree

                    # Get pixel scale - output in deg
                    image_scale = np.max(wcs.utils.proj_plane_pixel_scales(w1))

                    # size of image in arcseonds
                    size = round(image_scale * 3600 * np.nanmax(image.shape))

                    # Aligning images with astroalihn https://github.com/toros-astro/astroalign

                    template_found = False

                    if syntax['use_user_template']:
                        print('Using user template')
        #                print(syntax['fits_dir'] + '/templates')
                        if not os.path.exists(syntax['fits_dir'] + '/templates'):
                            print('Templates folder not found')
                        else:
                            print('Found user templates')
                            if not os.path.exists(syntax['fits_dir'] + '/templates/' + use_filter + '_template'):
                                print('cannot find User filter template')
                            else:
                                list_dir = os.listdir(syntax['fits_dir'] + '/templates/' + use_filter + '_template')
                                fits_list_dir = [i for i in list_dir if i.split('.')[-1] in ['fits','fts','fit']]

                                if len(fits_list_dir) >1:
                                    print(syntax['fits_dir'] + '/templates/' + use_filter + '_template')
                                    print(os.listdir(syntax['fits_dir'] + '/templates/' + use_filter + '_template'))
                                    print('Check template folder - too many files in fodler')
#                                else:
                                fpath_template = syntax['fits_dir'] + '/templates/' + use_filter + '_template/' + fits_list_dir[0]
#                                    subtraction_ready = True
                                template_found = True

                    if syntax_iter['get_template']:
                        print('Searching for template ...')
#
#                        # Template retrival from 2mass - shite because it's too faint
                        if syntax['catalog'] == '2mass':
                            try:
                                from astroquery.skyview import SkyView

                                # https://astroquery.readthedocs.io/en/latest/skyview/skyview.html
                                hdu  = SkyView.get_images(target_coords,survey = ['2MASS-'+use_filter.upper()],coordinates = 'ICRS',radius = size * u.arcsec)
                                fits.writeto(fpath.replace(fname_ext,'_template')+'no_rot'+fname_ext, hdu[0][0].data, headinfo, overwrite=True,output_verify = 'silentfix+ignore')

                            except Exception as e:
                                    if  syntax_iter['verbose']:
                                        sys.stdout = default_stdout
                                    exc_type, exc_obj, exc_tb = sys.exc_info()
                                    fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                                    print(exc_type, fname1, exc_tb.tb_lineno,e)


                        # Template retrival from panstarrs
                        if syntax['catalog'] == 'pan_starrs' or syntax['catalog'] == 'skymapper':

                            if os.path.isfile(syntax['fits_dir']+'/'+'templates/'+ use_filter + '_template/'+'template_'+use_filter+'_retrieved'+fname_ext):
                                print('Found previously retrived template:'+str(syntax['fits_dir'].split('/')[-1])+'/templates/')
                                template_found = True
                                fpath_template = syntax['fits_dir']+'/'+'templates/'+ use_filter + '_template/'+'template_'+use_filter+'_retrieved'+fname_ext
                            else:

                                pan_starrs_pscale = 0.25 # arcsec per pxiel

                                print('Searching for template on PanSTARRS')

                                fitsurl = get_pstars(float(ra), float(dec), size=int(size/pan_starrs_pscale), filters=use_filter)
                                sys.stdout = default_stdout
                                with fits.open(fitsurl[0],ignore_missing_end = True,lazy_load_hdus = True) as hdu:
                                    try:

                                        hdu.verify('silentfix+ignore')
                                        headinfo_template = hdu[0].header
                                        template_found  = True

                                        sys.stdout = default_stdout

                                        # save templates into original folder under the name template
                                        pathlib.Path(syntax['fits_dir']+'/'+'templates/'+ use_filter + '_template/').mkdir(parents = True, exist_ok=True)
                                        fits.writeto(syntax['fits_dir']+'/'+'templates/'+ use_filter + '_template/'+'template_'+use_filter+'_retrieved'+fname_ext, hdu[0].data,
                                                     headinfo_template,
                                                     overwrite=True,
                                                     output_verify = 'silentfix+ignore')
                                        if os.path.isfile(syntax['fits_dir']+'/'+'templates/'+ use_filter + '_template/'+'template_'+use_filter+'_retrieved'+fname_ext):
                                            print('Retrived template saved in: '+str(syntax['fits_dir'].split('/')[-1])+'/templates/')
                                            fpath_template = syntax['fits_dir']+'/'+'templates/'+ use_filter + '_template/'+'template_'+use_filter+'_retrieved'+fname_ext
                                            template_found = True

                                    except Exception as e:

                                        sys.stdout = default_stdout
                                        exc_type, exc_obj, exc_tb = sys.exc_info()
                                        fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                                        print(exc_type, fname1, exc_tb.tb_lineno,e)

                                sys.stdout = flog




                    if not template_found:
                        print('Template not found')
                    else:
                        with fits.open(fpath_template,ignore_missing_end = True,lazy_load_hdus = True) as hdu:
                            headinfo_template = hdu[0].header
                            try:
                                try:


                                    import astroalign as aa
                                    print('Aligning via Astro Align')
                                    aligned_image, footprint = aa.register(hdu[0].data.astype(float),image)
                                    aligned_image[np.isnan(aligned_image)] = 1e-5
                                    if  syntax_iter['verbose']:
                                        sys.stdout = default_stdout
                                except Exception as e:
                                    print('ASTRO ALIGN failed: %s' % e)

                                    print('Aligning via WCS')
                                    from reproject import reproject_interp
                                    aligned_image, footprint = reproject_interp(hdu[0], headinfo)
                                    aligned_image[np.isnan(aligned_image)] = 1e-5
                                    if syntax_iter['verbose']:
                                        sys.stdout = default_stdout



                                fits.writeto(fpath.replace(fname_ext,'_template')+fname_ext, aligned_image, headinfo_template, overwrite=True,output_verify = 'silentfix+ignore')
                                fpath_template = fpath.replace(fname_ext,'_template')+fname_ext
                                subtraction_ready = True

                            except Exception as e:
                                exc_type, exc_obj, exc_tb = sys.exc_info()
                                fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                                print(exc_type, fname1, exc_tb.tb_lineno,e)
#                                if  syntax_iter['verbose']:
#                                    sys.stdout = default_stdout

                    warnings.filterwarnings("default")
                    if os.path.isfile(fpath.replace(fname_ext,'_template')+fname_ext):
                        print('Template saved as:',os.path.basename(fpath.replace(fname_ext,'_template')))

                except Exception as e:
                    if syntax_iter['verbose']:
                            sys.stdout = default_stdout
                    print('*** Error with Template aquisiton ***')
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname1, exc_tb.tb_lineno,e)

# Check for template image

            if syntax_iter['get_template']:
                if not os.path.isfile(fpath.replace(fname_ext,'_template')+fname_ext):
                    print('Template file NOT found')

                else:
                    fpath_template = fpath.replace(fname_ext,'_template')+fname_ext
                    subtraction_ready = True

# =============================================================================
# Image subtraction using HOTPANTS
# =============================================================================

            if syntax['do_subtraction'] and subtraction_ready:
                print('\n --- Performing image subtraction using HOTPANTS --- ')
                fpath_sub = HOTterPANTS(fpath,fpath_template,syntax_iter)

                fig_sub = plt.figure(figsize = set_size(540,aspect = 1))

                ax = fig_sub.add_subplot(111)

                image_sub_tmp = fits.getdata(fpath_sub)

                norm = ImageNormalize( stretch = SquaredStretch())

                vmin,vmax = (ZScaleInterval(nsamples = 1200)).get_limits(image_sub_tmp)

                ax.imshow(image_sub_tmp,
                          vmin = vmin,
                          vmax = vmax,
                          norm = norm,
                          origin = 'lower',
                          cmap = 'Greys')

                ax.scatter([target_x_pix],[target_y_pix],marker = 'o',facecolor = 'None',
                           color = 'green',linewidth = 0.8,label = 'Target: %s' % tname)

                plt.legend()

                fig_sub.savefig(cur_dir + '/' +os.path.basename(fpath.replace(fname_ext,'_subtraction_QUICKLOOK'))+'.pdf')

                plt.close()

# =============================================================================
# Perform photometry on target
# =============================================================================


            if syntax_iter['do_phot']:

                if syntax_iter['phot_on_sub'] and subtraction_ready:
                    image    = fits.getdata(fpath_sub)
                    print('Target photometry performed on subtracted image')
                else:
                    print('Target photometry performed on original image')

                image_copy  = image.copy()

                print('\n --- Target Photometry --- ')

                target_x_pix_TNS, target_y_pix_TNS = w1.all_world2pix(syntax_iter['target_ra'], syntax_iter['target_dec'], 0)

                close_up = image_copy[int(target_y_pix_TNS)-syntax_iter['scale']: int(target_y_pix_TNS) + syntax_iter['scale'],
                                      int(target_x_pix_TNS)-syntax_iter['scale']: int(target_x_pix_TNS) + syntax_iter['scale']]


                x = np.arange(0,2*syntax_iter['scale'])
                xx,yy= np.meshgrid(x,x)
                sigma = mean_fwhm / (2*np.sqrt(2*np.log(2)))

                try:

                    pars = lmfit.Parameters()
                    pars.add('A',value = np.nanmax(close_up),min = 1e-5)
                    pars.add('x0',value = close_up.shape[0]/2)
                    pars.add('y0',value = close_up.shape[0]/2)
                    pars.add('sky',value = np.nanmedian(close_up))
                    pars.add('sigma',value = np.nanmedian(close_up))

                    def residual(p):
                        p = p.valuesdict()
                        return ( (close_up - gauss_2d((xx,yy),p['x0'],p['y0'],p['sky'],p['A'],sigma).reshape(close_up.shape))).flatten()

                    mini = lmfit.Minimizer(residual, pars)
                    result = mini.minimize(method = 'least_squares')

                    target_x_pix_corr= result.params['x0'].value
                    target_y_pix_corr= result.params['y0'].value

                    target_x_pix =  target_x_pix_corr  - syntax_iter['scale']  + target_x_pix_TNS
                    target_y_pix =  target_y_pix_corr  - syntax_iter['scale']  + target_y_pix_TNS




                except Exception as e:

                    print('- Could not fit gaussian profile to target - using TNS coords')
                    target_x_pix = target_x_pix_TNS
                    target_y_pix = target_y_pix_TNS

                    target_x_pix_corr = syntax_iter['scale']
                    target_y_pix_corr = syntax_iter['scale']

                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname, exc_tb.tb_lineno,e)
                    pass

                close_up = image[int(target_y_pix) - syntax_iter['scale']: int(target_y_pix) + syntax_iter['scale'],
                                 int(target_x_pix) - syntax_iter['scale']: int(target_x_pix) + syntax_iter['scale']]

                mean = np.nanmean(close_up)

                positions  = list(zip([target_x_pix],[target_y_pix]))

                target,target_bkg = ap_phot(positions,
                                             image,
                                             radius = syntax_iter['ap_size']    * mean_fwhm,
                                             r_in   = syntax_iter['r_in_size']  * mean_fwhm,
                                             r_out  = syntax_iter['r_out_size'] * mean_fwhm)

                target = (target/exp_time)[0]
                target_bkg = (target_bkg/exp_time)[0]

                SNR_target = SNR(target,target_bkg,exp_time,0,syntax_iter['ap_size']* mean_fwhm,gain,0)



                # if not do_ap and not syntax_iter['do_ap_phot']:

                    # if mag(target,0) - approx_psf_mag < 0.25:

                    #     print('PSF not applicable')
                    #     print('\n --- Approx tagret mag < Approx PSF mag --- ')
                    #     print(str(round(float(mag(target,0)),3)) ,'<',str(round(float(approx_psf_mag),3)))
                    #     do_ap = True
                    #     ap_corr = ap_corr_base

                if do_ap == True:

                    def residual(p):
                        p = p.valuesdict()
                        return ( (close_up - gauss_2d((xx,yy),p['x0'],p['y0'],p['sky'],p['A'],p['sigma']).reshape(close_up.shape))).flatten()

                    mini = lmfit.Minimizer(residual, pars)
                    result = mini.minimize(method = 'least_squares')
                    target_fwhm = result.params['sigma'].value * (2*np.sqrt(2*np.log(2)))


                if syntax_iter['save_target_plot'] and do_ap:

                    target_x_pix_corr = syntax_iter['scale']
                    target_y_pix_corr = syntax_iter['scale']

                    fig_target = plt.figure(figsize = set_size(540,aspect = 1))

                    ax = fig_target.add_subplot(111)

                    fig_target.suptitle('Target')

                    im = ax.imshow(close_up)

                    ax.scatter(target_x_pix_corr,target_y_pix_corr,label ='Location of target',marker = '*',color = 'red')
                    ax.plot([], [], ' ', label= 'x: '+'{0:.3f}'.format(target_x_pix)+' y: '+'{0:.3f}'.format(target_y_pix))
                    # ax.plot([], [], ' ', label= 'dec '+str(object_info['dec']))

                    circ1 = Circle((target_x_pix_corr,target_y_pix_corr),syntax_iter['ap_size'] * mean_fwhm,label = 'Integration Radius',fill=False,color = 'red')
                    circ2 = Circle((target_x_pix_corr,target_y_pix_corr),syntax_iter['r_in_size'] * mean_fwhm,label = 'Inner Annulus',fill=False,color = 'red',linestyle = '--')
                    circ3 = Circle((target_x_pix_corr,target_y_pix_corr),syntax_iter['r_out_size'] * mean_fwhm,label = 'Outer Annulus',fill=False,color = 'red',linestyle = '--')
                    ax.add_patch(circ1)
                    ax.add_patch(circ2)
                    ax.add_patch(circ3)

                    plt.legend(loc = 'best')

                    divider = make_axes_locatable(ax)

                    cax = divider.append_axes("right", size="5%", pad=0.05)
                    fig_target.colorbar(im, cax=cax)

                    fig_target.savefig(cur_dir + '/' +'target_'+str(base.split('.')[0])+'.pdf')

                    plt.close(fig_target)

                if syntax_iter['do_psf_phot'] and not do_ap:

                    print('Perforoming PSF Photometry on Target')

                    c_target = pd.DataFrame(data = [[target_x_pix,target_y_pix]],columns = ['x_pix','y_pix'])

                    c_psf_target1,_ = psf.fit(image,
                                              c_target,
                                              r_table,
                                              syntax_iter,
                                              mean_fwhm,
                                              save_plot = syntax_iter['save_target_plot'],
                                              show_plot = syntax_iter['plot_target'] ,
                                              fname = str(base.split('.')[0]),
                                              return_fwhm = True)


                    c_psf_target,_ = psf.do(c_psf_target1,
                                            r_table,
                                            syntax_iter,
                                            mean_fwhm)

                    target = np.array(c_psf_target.psf_counts/exp_time)[0]


                    target_bkg = np.array(c_psf_target.bkg/exp_time)[0]


                    target_err = c_psf_target.psf_counts_err/exp_time

                    SNR_target = np.array(target/target_err)[0]

                    target_fwhm = c_psf_target['target_fwhm'].values[0]
                else:
                    print(' > Doing Aperture Photometry on Target <')

# =============================================================================
# Limiting Magnitude
# =============================================================================



                try:
                    syntax_iter['image_radius']
                except:
                    syntax_iter['image_radius'] = 1.3 * mean_fwhm

                if syntax['get_lim_mag_test']:

                    if SNR_target > 10 and abs(target_fwhm- mean_fwhm) < 1:
                            lmag = np.nan
                            output.update({'lmag':lmag})

                            print('SNR = %.f - skipping limiting magnitude' % SNR_target)
                    else:
                        print('SNR = %.f - checking limiting magnitude' % SNR_target)
                        print('Discrepancy in FWHM of %.1f pixels' % abs(target_fwhm - mean_fwhm))

                        expand_scale = np.floor(3*syntax_iter['image_radius'])

                        # cut out an Expanded area around target x and y pixel - needs to be even
                        if (expand_scale - syntax_iter['scale'] )%2 !=0:
                            expand_scale-=1

                        # get close up at bigger scale to allow for injected sources to be well spaced
                        close_up_expand = image_copy[int(target_y_pix - expand_scale): int(target_y_pix + expand_scale),
                                                     int(target_x_pix - expand_scale): int(target_x_pix + expand_scale)]

                        try:
                             model = model_psf
                             r_table = r_table
                        except:
                            model = None
                            r_table = None

                        lmag_inst,low_lmag_inst = limiting_magnitude_prob(syntax_iter,close_up_expand,model,r_table)

                        lmag = lmag_inst + zp_wa[0]

                        # Wost case scenario
                        # low_lmag = lmag_inst + zp_wa[0]

                        output.update({'lmag':lmag})
                        print('Limiting Magnitude: %.1f' % lmag)

                print('\n--- Output ---\n')


# =============================================================================
# Error on target magnitude
# =============================================================================

                # Error from SNR

                SNR_error = sigma_mag_err(SNR_target)

                fit_error  = mag(target,0) - mag(target+target_err,0)

                target_mag_err = SNR_error + fit_error[0]


# =============================================================================
# Output
# =============================================================================

                mag_target={use_filter:mag(target,zp_wa[0]) + ap_corr}
                mag_inst={use_filter+'_inst':mag(target,0)}

                mag_inst_err={use_filter+'_inst_err':mag(target_mag_err,0)}
                zp={'zp_'+use_filter:zp_wa[0]}
                zp_err={'zp_'+use_filter+'_err':zp_wa[1]}
                fwhm_out={'fwhm':mean_fwhm}

                target_fwhm_out = {'target_fwhm':target_fwhm}



                mag_err = {use_filter+'_err':np.sqrt(target_mag_err**2 + zp_wa[1]**2)}
                SNR_dict = {'SNR':SNR_target}

                time_exe = {'time':str(datetime.datetime.now())}

                if mag_target[use_filter] == 0.0:
                    print(' Target not seen - SNR = 0 ')
                    mag_target[use_filter] = [np.nan]

                output.update(mag_inst)
                output.update(mag_inst_err)
                output.update(zp)
                output.update(zp_err)

                output.update(mag_target)
                output.update(mag_err)
                output.update(SNR_dict)
                output.update(time_exe)
                output.update(fwhm_out)
                output.update(target_fwhm_out)

                if do_ap:
                    output.update({'method':'ap'})
                else:
                    output.update({'method':'psf'})

                print('Target counts: %.5f' % target)
                print('Target counts error: %.5f' % target_err)

                print('Target SNR: %.5f' % SNR_target)
                print('Target SNR Error: %.5f' % SNR_error)

                print('Instrumental Magnitude: %.3f' % mag(target,0)[0])
                print('Instrumental Magnitude error: %.3f' % fit_error)
                print('Limiting Magnitude: %.3f' % lmag)
                print('Target Magnitude: %.3f +/- %.3f ' % (mag_target[use_filter][0],
                                                            mag_err[use_filter+'_err']) )

                if lmag < mag_target[use_filter][0] and not np.isnan(lmag) or np.isnan(mag_target[use_filter][0]):
                        print('\n*** Image is magnitude limited ***')


                # Calibration file used in reduction
                c.round(6).to_csv(syntax_iter['write_dir']+'image_calib_'+str(base.split('.')[0])+'_filter_'+str(use_filter)+'.csv',index = False)


                target_output = pd.DataFrame(output,columns=output.keys(), index=[0])
                target_output.to_csv(cur_dir + '/' +'out.csv',index=False)

                for key,val in zp.items():
                        try:
                            zp_wa
                        except:
                            pass

                sys.stdout = default_stdout
                flog.close()




            flog.close()

            print('\n' +'Time Taken [ %s ]: %ss' % (str(os.getpid()),round(time.time() - start)))
            print('Sucess: '+ str(base) + ' - PID: '+str(os.getpid()))


            if syntax['do_all_phot']:
                    print('Performing Photometry on all sources')

                    _,df_all,_ = find.fwhm(image,
                                           syntax_iter,
                                           sigma_lvl = syntax_iter['do_all_phot_sigma'],
                                           fwhm = mean_fwhm)

                    print('Number of sources: %s' % len(df_all))

                    print('Saving PSF subtracted images(s)')

                    syntax_iter['show_residuals'] = True


                    _,_ = psf.fit(image,
                                  df_all,
                                  r_table,
                                  syntax_iter,
                                  mean_fwhm,
                                  return_subtraction_image = True)



            return output,base

        # Parent try/except statement for loop
        except Exception as e:

            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname1, exc_tb.tb_lineno,e)
            sys.stdout = default_stdout
            flog.close()

            print('\n' +'Failure: '+ str(base) + '- PID: '+str(os.getpid()))

            return None,fpath

    except Exception as e:
        print('Fatal Error')

        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname1, exc_tb.tb_lineno,e)

        sys.stdout = default_stdout
        flog.close()
        print('\n' +'Failure: '+str(base) + '- PID: '+str(os.getpid()))

        files = (file for file in os.listdir(os.getcwd())
             if os.path.isfile(os.path.join(os.getcwd(), file)))

        for file in files:
             if fpath == file:
                 break

             if os.path.basename(fpath).split('.')[0] in file:
                   shutil.move(os.path.join(os.getcwd(), file),
                    os.path.join(cur_dir, file))

        try:

            c.to_csv('table_ERROR_'+str(base.split('.')[0])+'_filter_'+str(use_filter)+'.csv',index=False)
        except:
            pass

        return None,fpath



