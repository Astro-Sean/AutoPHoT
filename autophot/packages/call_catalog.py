#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 09:40:29 2018

@author: seanbrennan
"""

def search(image, headinfo, target_coords, syntax, catalog_syntax, filter_):

    """
    Search area around transient/target location in photometric catalogs

    Current catalog (selectable in syntax):

        - Skymapper: Southern Hemisphere
        - Pan Starrs: North of declination -30 degree
        - Apass: All-sky survey
        - 2MASS: JHK all sky survey

    Future:
        - SDSS: Future implemtation
        - Ability to make custom catalog from different surveys


    Input:

        - Image: Numpy 2D array
        - headinfo: astropy.io.fits.header.Header
        - target_coords: astropy.coordinates.sky_coordinate.SkyCoord
        - syntax: dict
        - catalog_syntax: dict
        - filter_: str

    Output:

        - data:  pandas DataFrame

    """

    import warnings
    if not syntax['catalog_warnings'] or syntax['master_warnings']:
            warnings.filterwarnings("ignore")

    import numpy as np
    import os,sys,inspect
    import requests
    import pathlib
    import shutil
    import os.path

    from astropy.table import Table
    from astropy.wcs import wcs
    from astroquery.vizier import Vizier
    from astropy.io.votable import parse_single_table
    from astropy.coordinates import Angle

    try:

        # Get wxs information
        w1 = wcs.WCS(headinfo)

        # Radius arounbd target
        radius  = float(syntax['radius'])

        # Target name, if applicabel
        target = syntax['target_name']

        # Get current directory location,, create directory if needed
        currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        dir_name = currentdir+'/'+'catalog_queries'
        pathlib.Path(dir_name).mkdir(parents = True, exist_ok=True)



        '''
        Getting target Ra and Dec

        - if target is none but a Ra and DEC is given, create new target name

        - if ra and dec not given us center of image as location - for quick reduction of image

        '''

        if target == None:
             if syntax['target_ra'] != None:
                 target = 'target_ra_'+str(round(syntax['target_ra']))+'_dec_'+str(round(syntax['target_dec']))
                 print('New target name: %s' %target)
             else:
                 target = 'target'


        # Search limitation with Pan Starrs rlimited to 0.5deg
        if radius > 0.5 and syntax['catalog'] == 'pan_starrs' :
                    print('WARNING: Search Limitation with PanStarrs API -> Radius = 0.5 [deg] ')
                    radius = 0.5

        # Choosen catalog for input.yml, create directory if needed
        catalog_dir = syntax['catalog']
        pathlib.Path(dir_name + '/' + catalog_dir).mkdir(parents = True, exist_ok=True)

        # Folder for target, create directory if needed
        target_dir =  dir_name + '/' + catalog_dir+'/'+target.lower()
        pathlib.Path(target_dir).mkdir(parents = True, exist_ok=True)

        # Filename of fetchec catalog
        fname = str(target) + '_RAD_' + str(radius)

        # Can force to use certain catalog - untested 03-10-19
        if syntax['force_catalog_csv']:
            print('Using '+syntax['force_catalog_csv_name']+' as catalog')
            fname = str(syntax['force_catalog_csv_name']) + '_RAD_' + str(radius)

        # if catalog is found via it's filename - use this and return data
        if os.path.isfile(target_dir +'/'+ fname+'.csv'):
            print('Catalog found for Target: %s\nCatalog: %s \nFile: %s' % (target,str(catalog_dir).upper(),fname))

            data = Table.read(target_dir+'/'+ fname+'.csv',format = 'csv')
            data = data.to_pandas()

        else:
            # If no previously catalog found
            print('Searching for new catalog: %s ' % syntax['catalog'])

            if syntax['catalog'] in ['gaia']:

                import astropy.units as u
                from astroquery.gaia import Gaia
                import warnings
                warnings.filterwarnings('ignore')

                width = u.Quantity(radius, u.deg)
                height = u.Quantity(radius, u.deg)

                data = Gaia.query_object_async(coordinate=target_coords, width=width, height=height)

                data = data.to_pandas()
                data.to_csv(fname+'.csv', sep=',',index = False)
                print(len(data))
                # Move file to new location - 'catalog queries'
                shutil.move(os.path.join(os.getcwd(), fname+'.csv'),
                            os.path.join(target_dir, fname+'.csv'))

                warnings.filterwarnings('default')

            if syntax['catalog'] in ['apass','2mass']:
                print('Downloading from %s'  % syntax['catalog'] )

                # No row limit
                Vizier.ROW_LIMIT = -1
                catalog_search = Vizier.query_region(target_coords,
                                                     radius = Angle(radius,'deg'),
                                                     catalog = syntax['catalog'])

                # Select first catalog from list
                data = catalog_search[0].to_pandas()
                data.to_csv(fname+'.csv', sep=',',index = False)


                # Move file to new location - 'catalog queries'
                shutil.move(os.path.join(os.getcwd(), fname+'.csv'),
                            os.path.join(target_dir, fname+'.csv'))


            # some catalogs need specific download path using 'requests'
            if syntax['catalog'] in ['pan_starrs','skymapper']:

                mindet=1

                if syntax['catalog'] == 'pan_starrs':


                    server=('https://archive.stsci.edu/'+'panstarrs/search.php')
                    params = {'RA': target_coords.ra.degree, 'DEC': target_coords.dec.degree,
                              'SR': radius, 'max_records': 10000,
                              'outputformat': 'VOTable',
                              'ndetections': ('>%d' % mindet)}

                if syntax['catalog'] == 'skymapper':

                    server=('http://skymapper.anu.edu.au/sm-cone/public/query?')
                    params = {'RA': target_coords.ra.degree, 'DEC': target_coords.dec.degree,
                              'SR': radius,
                              'RESPONSEFORMAT': 'VOTABLE'}

                with open('temp.xml', "wb") as f:

                    print('Downloading from %s'  % syntax['catalog'] )
                    response = requests.get(server,params = params)
                    f.write(response.content)

                # Parse local file into astropy.table object
                data = parse_single_table('temp.xml')

                # Delete temporary file
                os.remove('temp.xml')

                # Convert table to dataframe
                data_table = data.to_table(use_names_over_ids=True)
                data = data_table.to_pandas()

                # invalid entries in panstarrs are -999 - change to nans
                if syntax['catalog'] == 'pan_starrs':
                    data = data.replace(-999,np.nan)

                # No sources in field - temporary fix - will add "check different catalog"
                if len(data) == 0:
                    sys.exit('WARNING: %s No Coverage - Select new catalog' %  syntax['catalog'])

                # Save to csv and move to 'catalog_queries'
                data.to_csv(fname+'.csv',index = False)

                shutil.move(os.path.join(os.getcwd(), fname+'.csv'),
                            os.path.join(target_dir, fname+'.csv'))


        # Add in x and y pixel locatins under wcs
        x_pix,y_pix = w1.wcs_world2pix(data[catalog_syntax['RA']], data[catalog_syntax['DEC']],1)

        data.insert(loc = 5, column = 'x_pix', value = x_pix)
        data.insert(loc = 6, column = 'y_pix', value = y_pix)

        # Remove boundary sources
        data = data[data.x_pix < image.shape[1] - syntax['pix_bound']]
        data = data[data.x_pix > syntax['pix_bound']]
        data = data[data.y_pix < image.shape[0] - syntax['pix_bound']]
        data = data[data.y_pix > syntax['pix_bound']]

        print('Data length: %.f' % len(data))

        warnings.filterwarnings("default")

    except Exception as e:
         exc_type, exc_obj, exc_tb = sys.exc_info()
         fname1 = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
         print(exc_type, fname1, exc_tb.tb_lineno,e)
         data = None

    return data



def match(image, headinfo, target_coords, syntax, catalog_syntax, filter_,data, fwhm):

    """
    Match positions from catalog with locations on image to check for source aboove threshold level given by
    'bkg_level' in syntax


    Input:

        - image: Numpy 2D array
        - headinfo: astropy.io.fits.header.Header
        - target_coords: astropy.coordinates.sky_coordinate.SkyCoord
        - syntax: dict
        - catalog_syntax: dict
        - filter_: str
        - data: pandas DataFrame
        - fwhm: float

    Output:

        - data_new_frame:  pandas DataFrame

    """

    import warnings
    if not syntax['catalog_warnings'] or syntax['master_warnings']:
            warnings.filterwarnings("ignore")

    import numpy as np
    import matplotlib.pyplot as plt

    import os,sys
    import pandas as pd
    import lmfit

    from astropy.stats import sigma_clipped_stats
    from photutils import DAOStarFinder

    from autophot.packages.functions import r_dist,gauss_2d

    x_new_source = []
    y_new_source = []
    x_new_cen    = []
    y_new_cen    = []
    cp_dist      = []
    cat_idx      = []
    non_detections = []
    detections     = []
    fwhm_list = []

    # Remove values that don't have matching value in selected value
    data_update = data[~np.isnan(data[catalog_syntax[filter_]])]

    # Look at most accuarte measurements first based on errors
    try:
        data_update.sort_values(by=[catalog_syntax[filter_+'_err']],inplace = True, na_position='last')
    except:
        data_update[filter_+'_err'] = np.nan * len(data_update)
        pass


    # Grid of close-up mathcing scale of image
    x = np.arange(0,2*syntax['scale'])
    xx,yy= np.meshgrid(x,x)

    k = 0

    # Wiggle room for catalog matching
    dx = syntax['catalog_matching_dx']
    dy = syntax['catalog_matching_dy']

    useable_sources = 0

    try:
        print('\n --- Matching sources to catalog --- ')

        for idx in data_update.index:
            if useable_sources >= syntax['max_catalog_sources']:
                break


            try:

                 # Skip if source location is off the image
                 if data_update.x_pix[idx] <= 0 or data_update.y_pix[idx] <= 0:
                     x_new_cen.append(np.nan)
                     y_new_cen.append(np.nan)
                     cp_dist.append(np.nan)
                     continue

                 # catalog pixel coordinates of source take as an approximate location
                 y = data_update.y_pix[idx]
                 x = data_update.x_pix[idx]

                 # Add index key for original catalog file comparision and matching
                 cat_idx.append(int(idx))

                 # Add x and y pxel location
                 x_new_source.append(x)
                 y_new_source.append(y)

                 # Create cutout image of size (2*syntax['scale'],2*syntax['scale'])
                 close_up = image[int(y)-syntax['scale']: int(y) + syntax['scale'],
                                  int(x)-syntax['scale']: int(x) + syntax['scale']]


                 # Cutout not possible - too close to edge or invalue pixel data i.e. nans of infs
                 if close_up.shape != (2*syntax['scale'],2*syntax['scale']):
                     x_new_cen.append(np.nan)
                     y_new_cen.append(np.nan)
                     cp_dist.append(np.nan)
                     fwhm_list.append(np.nan)

                     continue

                 # Preset pixel error popup skip this source
                 if 1e-5 in close_up or np.isnan(np.min(close_up)):
                     x_new_cen.append(np.nan)
                     y_new_cen.append(np.nan)
                     cp_dist.append(np.nan)
                     fwhm_list.append(np.nan)

                     continue

                 # Get close up image properties
                 mean, median, std = sigma_clipped_stats(close_up,
                                                         sigma = syntax['source_sigma_close_up'],
                                                         maxiters = syntax['iters'])


                 # Preform source detection with threshold set in input.yml
                 daofind = DAOStarFinder(fwhm = fwhm,
                                         threshold = syntax['bkg_level']*std,
                                         roundlo= -1.0, roundhi=1.0,
                                         sharplo = 0.2, sharphi=1.0)


                 sources = daofind(close_up - median)

                 if sources == None:
                     sources = []

                 if len(sources) == 0 :

                     non_detections.append(data_update[catalog_syntax[filter_]].loc[[idx]].values[0])

                     x_new_cen.append(np.nan)
                     y_new_cen.append(np.nan)
                     cp_dist.append(np.nan)
                     fwhm_list.append(np.nan)

                     continue


                 # Approximate location of source
                 xc_guess = np.array(sources['xcentroid'])[0]
                 yc_guess = np.array(sources['ycentroid'])[0]

                 # If more than one source detected in close up
                 # assume source closest to center is desired source
                 if len(sources) >=2:
                     r_vals = r_dist(syntax['scale'],np.array(sources['xcentroid']),
                                     syntax['scale'],np.array(sources['ycentroid']))

                     r_idx = np.argmin(r_vals)

                     # if closest source is too far away from predicted loction - ignore
                     if r_vals[r_idx] > syntax['match_dist']:
                         x_new_cen.append(xc_guess)
                         y_new_cen.append(yc_guess)
                         cp_dist.append(np.nan)
                         fwhm_list.append(np.nan)
                         continue

                     xc_guess = np.array(sources['xcentroid'])[r_idx]
                     yc_guess = np.array(sources['ycentroid'])[r_idx]
                 try:

                     # Fit gaussian to source for centroiding location
                     pars = lmfit.Parameters()

                     pars.add('A', value = np.max(close_up),min = 0)
                     pars.add('x0',value = xc_guess,min = 0.5*close_up.shape[0]-dx,max = 0.5*close_up.shape[0]+dx)
                     pars.add('y0',value = yc_guess,min = 0.5*close_up.shape[0]-dy,max = 0.5*close_up.shape[0]+dy)
                     pars.add('sky',value = np.median(close_up))
                     pars.add('sigma',value = fwhm/2*np.sqrt(2*np.log(2)) )

                     def residual(p):
                        p = p.valuesdict()
                        return (close_up - gauss_2d((xx,yy),p['x0'],p['y0'],p['sky'],p['A'],p['sigma']).reshape(close_up.shape)).flatten()

                     mini = lmfit.Minimizer(residual, pars)
                     result = mini.minimize(method = 'leastsq')

                     xcen = result.params['x0'].value
                     ycen = result.params['y0'].value

                     S    = result.params['sky'].value
                     H    = result.params['A'].value

                     sigma =  result.params['sigma'].value

                 except Exception as e:

                     x_new_cen.append(np.nan)
                     y_new_cen.append(np.nan)
                     cp_dist.append(np.nan)
                     fwhm_list.append(np.nan)


                     exc_type, exc_obj, exc_tb = sys.exc_info()
                     fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                     print(exc_type, fname, exc_tb.tb_lineno,e)

                     continue

                 k+=1

                 # Add new source location accounting for difference in fitted location / expected location
                 x_new_cen.append(xcen -   syntax['scale'] + x)
                 y_new_cen.append(ycen -   syntax['scale'] + y)

                 cp_dist.append(np.sqrt( (xcen - syntax['scale'])**2 + (ycen - syntax['scale'])**2) )

                 detections.append(data_update[catalog_syntax[filter_]].loc[[idx]].values[0])

                 fwhm_list.append(sigma * 2*np.sqrt(2*np.log(2)))

                 useable_sources +=1

                 if syntax['source_plot']:
                     if len(sources) == 1 :

                          fig = plt.figure(figsize = (6,6))
                          ax = fig.add_subplot(111)
                          ax.imshow(close_up)

                          ax.set_title('Source @ x = '+'{0:.3f}'.format(xcen +  x - syntax['scale'])+' : y = '
                                      +'{0:.3f}'.format(ycen +  y - syntax['scale']))

                          small_ap = plt.Circle((xcen,ycen), syntax['ap_size']*headinfo['FWHM'], color='r',fill = False,label = 'Aperture')
                          big_ap = plt.Circle((xcen,ycen), syntax['inf_ap_size']*headinfo['FWHM'], color='b',fill = False,label = 'Aperture Correction')

                          ax.add_artist(small_ap)
                          ax.add_artist(big_ap)

                          ax.plot([],[],' ', label = 'Sky ='+'{0:.3f}'.format(S)+'Height ='+'{0:.3f}'.format(H))
                          ax.scatter(syntax['scale'],syntax['scale'],marker = '+',s = 100 , color = 'r',linewidths=0.01,label = 'Catalog')
                          ax.scatter(xc_guess,yc_guess,marker = '+',s = 100 , color = 'b',linewidths=0.01,label= 'Source detection [closest object to catalog]' )
                          ax.scatter(xcen,ycen,marker = '+',s = 100 , color = 'green',linewidths=0.01,label= 'Least square fit' )

                          ax.legend(loc = 'upper right')
                          plt.tight_layout()
                          plt.close()

            except Exception as e:

                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno,e)

                x_new_cen.append(xc_guess)
                y_new_cen.append(yc_guess)
                cp_dist.append(np.nan)
                fwhm_list.append(np.nan)

                continue

        if syntax['show_nondetect_plot']:

            non_detections = np.array(non_detections)[np.isfinite(non_detections)]
            detections = np.array(detections)[np.isfinite(detections)]


            if len(non_detections) == 0:
                print('All sources detected')

            fig = plt.figure(figsize=(6,8))
            ax = fig.add_subplot(111)
            ax.hist(non_detections,bins = 'auto',
                         align = 'mid',
                         color = 'green',
                         histtype = 'step',
                         label = 'Non-Detection')

            ax.hist(detections,
                         bins = 'auto',
                         align = 'mid',
                         color = 'red',
                         histtype = 'step',
                         label = 'Detection')
            ax.set_title('Non - Detections')
            ax.set_xlabel('Magnitude')
            ax.set_ylabel('Binned Occurance')
            ax.legend(loc = 'best')
            plt.show()


        frame_data = [np.array(cat_idx).astype(int),
                      np.array(data_update[catalog_syntax['RA']]),
#                      np.array(data_update[catalog_syntax['RA_err']]),
                      np.array(data_update[catalog_syntax['DEC']]),
                      # np.array(data_update[catalog_syntax['DEC_err']]),
                      np.array(x_new_source),
                      np.array(y_new_source),
                      np.array(x_new_cen),
                      np.array(y_new_cen),
                      np.array(cp_dist),
                      np.array(fwhm_list)
                      ]

        frame_cols = ['cat_idx',
                      'ra',
#                      'ra_err',
                      'dec',
#                      'dec_err',
                      'x_pix_source',
                      'y_pix_source',
                      'x_pix',
                      'y_pix',
                      'cp_dist',
                      'fwhm'
                      ]
        # filters = ['J','H','K','U', 'u','v', 'B', 'g', 'V', 'r', 'R', 'i', 'I', 'z']

        data_new_frame = pd.DataFrame(frame_data).T
        data_new_frame.columns = frame_cols
        data_new_frame.set_index('cat_idx',inplace = True)

#         for col in filters:
#             try:
#                 if catalog_syntax[col] in data_update.columns:
#                     data_new_frame[col]        = data_update[catalog_syntax[col]]
#                     data_new_frame[col+'_err'] = data_update[catalog_syntax[col+'_err']]



#             except:

# #                exc_type, exc_obj, exc_tb = sys.exc_info()
# #                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
# #                print(exc_type, fname, exc_tb.tb_lineno,e)

#                 x_new_cen.append(np.nan)
#                 y_new_cen.append(np.nan)
#                 cp_dist.append(np.nan)
#                 fwhm_list.append(np.nan)
#                 pass

        data_new_frame[filter_]        = data_update[catalog_syntax[filter_]]
        data_new_frame[filter_+'_err'] = data_update[catalog_syntax[filter_+'_err']]

        data_new_frame = data_new_frame[~np.isnan(data_new_frame['x_pix'])]
        data_new_frame = data_new_frame[~np.isnan(data_new_frame['y_pix'])]

        print('Number of sources from catalog: ',len(data_new_frame))

        warnings.filterwarnings("default")

    except Exception as e:

                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno,e)
                cp_dist.append( np.nan )

    return data_new_frame,syntax
