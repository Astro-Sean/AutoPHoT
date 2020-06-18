#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Jan 29 12:24:07 2019

@author: seanbrennan
"""
import matplotlib.pyplot as plt

plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
plt.rcParams['font.size'] = 7
plt.rcParams['figure.subplot.hspace']= 0.1
plt.rcParams['figure.subplot.wspace']= 0.1
plt.rcParams['figure.dpi']= 100
plt.rcParams['axes.titlesize'] = 7
plt.rcParams['axes.labelsize'] = 7

plt.rcParams['lines.linewidth'] = 1
plt.rcParams['savefig.format'] = 'pdf'
plt.rcParams['lines.markersize'] = 3

plt.rcParams['legend.markerscale'] = 2
plt.rcParams['legend.fontsize'] = 7

plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['xtick.minor.width'] = 0.35
plt.rcParams['xtick.direction'] = 'in'

plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['ytick.minor.width'] = 0.35
plt.rcParams['ytick.direction'] = 'in'

import numpy as np
from scipy.optimize import least_squares
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def gauss_2d(image, x0,y0, sky , A, sigma):
    import numpy as np
    (x,y) = image
    a = (x-x0)**2
    b = (y-y0)**2
    c = (2*sigma**2)
    d =  A*np.exp( -(a+b)/c)
    e =  d + sky
    return  e.ravel()

def pixel_correction(x,m):
    diff =  x - int(x)
    if diff/m >= 0.5 or  diff ==0:
        return np.ceil(x) - 0.5
    if diff/m < 0.5:
        return np.floor(x)

# Whatever pixel (0,0 TL),coordinate on the image to the corrosponding corrdinate in center of pixel
def array_correction(x):
    diff =  float(x) % 1
    if diff >= 0.5:
        return np.ceil(x)
    if diff < 0.5:
        return np.floor(x)

def norm(array):
    import numpy as np
    norm_array = (array - np.min(array))/(np.nanmax(array)-np.min(array))
#    norm_array = array/np.nanmax(array)
    return norm_array

def renorm(array,lb,up):
    s = up - lb
    n =  (array - np.min(array))/(np.nanmax(array)-np.min(array))
    m = (s * n) + lb
    return m

def find_2d_int_percent(count_percent,fwhm):

    from scipy.integrate import dblquad
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))

    A = 1 / (2 * np.pi * sigma **2)
    gauss = lambda y,x: A * np.exp( -1 * (((y)**2+(x)**2)/(2*sigma**2)))
    fit = lambda x0: (count_percent - dblquad(gauss, -1*x0,x0,lambda x: -1*x0,lambda x: x0)[0])

    r = least_squares(fit,x0 = 3)
    return r.x[0]

def get_an(image,xc,yc,r_in = None,r_out = None):
    from photutils import CircularAnnulus

    if r_in == None:
        r_in = 0.95
    if r_out == None:
        r_in = 0.95

    ap_an = CircularAnnulus((xc,yc),r_in = r_in ,r_out =  r_out)
    mask_annulus = ap_an.to_mask(method='subpixel',subpixels=5 )
    mask_annulus = mask_annulus[0].to_image(shape=((image.shape)))
    mask_annulus /= mask_annulus.max()
    an =  (abs(image * mask_annulus))

    return an

def scale_roll(x,xc,m):
    dx = (x - xc)
    if m !=1:
        shift = int(round(dx *m))
    else:
        shift = int(dx *m)
    return shift

def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)



# =============================================================================
# Point spread function by Me :)
# =============================================================================

def build_r_table(base_image,selected_sources,syntax,fwhm):

    '''
    Build psf model:
        Fits gaussian to bright sources and removes fit gaussian to leave residual.
        Residual signifies how the model gaussian deviates from the true psf of the image
        Build up averaged table of these residual images to constuct a psf
        with a gaussian core with this residual table
    '''

    import numpy as np
    from astropy.stats import sigma_clipped_stats
    from photutils import DAOStarFinder


    import pandas as pd
    import lmfit
    import logging

    import sys,os
    import matplotlib.pyplot as plt
    import pathlib
    from astropy.stats import sigma_clip
    from astropy.stats import SigmaClip
    from photutils import Background2D, MedianBackground

    from autophot.packages.functions import r_dist
    from autophot.packages.uncertain import SNR
    from autophot.packages.aperture  import ap_phot

    try:

        logger = logging.getLogger(__name__)

        image = base_image.copy()

        fitting_radius = int(np.ceil(fwhm))

        regriding_size = int(syntax['regrid_size'])
        m = regriding_size

        sigma_fit = []

        construction_sources = []

        if regriding_size % 2 > 0:
            logger.info('regrid size must be even adding 1')
            regriding_size += 1

        # Residual Table in extended format
        residual_table = np.zeros((2 * syntax['scale'] * regriding_size, 2 * syntax['scale']*regriding_size))

        selected_sources['dist'] =  r_dist(syntax['target_x_pix'],selected_sources.x_pix,
                                           syntax['target_y_pix'],selected_sources.y_pix)

        lower_bkg_source = sigma_clip(selected_sources['median'],
                                      sigma=3,
                                      maxiters=3,
                                      masked = True)

        selected_sources = selected_sources[~lower_bkg_source.mask]

        flux_idx = [i for i in selected_sources.flux_ap.sort_values(ascending = False).index]
        flux_idx = [i for i in selected_sources.dist.sort_values(ascending = True).index]

        sources_used = 1
        n = 0
        failsafe = 0
        psf_mag = []
        image_radius_lst = []

        sigma  = fwhm/(2*np.sqrt(2*np.log(2)))

        # level
        level = syntax['lim_SNR']

        logger.info('Limiting threshold: %s sigma' % level)


        while sources_used <= syntax['psf_source_no']:

            if failsafe>25:
                logger.info('PSF - Failed to build psf')
                residual_table=None
                sigma_fit = fwhm


            if n >= len(flux_idx):
                if sources_used  >= syntax['min_psf_source_no']:
                    logger.info('Using worst case scenario number of sources')
                    break
                logger.info('PSF - Ran out of sources')
                residual_table=None
                sigma_fit = fwhm
                break
            try:

                idx = flux_idx[n]

                n+=1

                psf_image = image[int(selected_sources.y_pix[idx])-syntax['scale']: int(selected_sources.y_pix[idx]) + syntax['scale'],
                                  int(selected_sources.x_pix[idx])-syntax['scale']: int(selected_sources.x_pix[idx]) + syntax['scale']]
                if len(psf_image) == 0:
                    logger.info('PSF image ERROR')
                    continue
                if np.min(psf_image) == np.nan:
                    continue


                mean, median , std = sigma_clipped_stats(psf_image, sigma = syntax['source_sigma_close_up'], maxiters = syntax['iters'])

                daofind = DAOStarFinder(fwhm=np.floor(fwhm),
                                     threshold = syntax['bkg_level']*std,
                                     roundlo = -1.0, roundhi = 1.0,
                                     sharplo =  0.2, sharphi = 1.0)


                sources = daofind(psf_image - median)

                if sources is None:
                    sources = []

                if len(sources) > 1:

                    dist = [list(r_dist(
                            sources['xcentroid'][i]
                            ,sources['xcentroid']
                            ,sources['ycentroid'][i]
                            ,sources['ycentroid']) for i in range(len(sources)))]

                    dist = np.array(list(set(np.array(dist).flatten())))

                    if all(dist < 2):
                        pass
                    else:
                        continue

                psf_image = image[int(selected_sources.y_pix[idx])-syntax['scale']: int(selected_sources.y_pix[idx]) + syntax['scale'],
                                  int(selected_sources.x_pix[idx])-syntax['scale']: int(selected_sources.x_pix[idx]) + syntax['scale']]



                if syntax['psf_bkg_surface']:


                    sigma_clip = SigmaClip(sigma=syntax['Lim_SNR'])
                    bkg_estimator = MedianBackground()
                    bkg = Background2D(psf_image, (3, 3),
                                       filter_size=(5, 5),
                                       sigma_clip=sigma_clip,
                                       bkg_estimator=bkg_estimator)


                    surface = bkg.background

                    # bkg_median = np.nanmedian(bkg.background_median)
                    psf_image_bkg_free = psf_image- surface

                if syntax['psf_bkg_poly']:


                    # Fit 2D surface to close_up
                    from astropy.modeling import models, fitting
                    surface_function_init = models.Polynomial2D(degree=syntax['psf_bkg_poly_degree'])

                    fit_surface = fitting.LevMarLSQFitter()


                    x = np.arange(0,psf_image.shape[0])
                    y = np.arange(0,psf_image.shape[1])
                    xx,yy= np.meshgrid(x,y)


                    with warnings.catch_warnings():
                        # Ignore model linearity warning from the fitter
                        warnings.simplefilter('ignore')
                        surface_fit = fit_surface(surface_function_init, xx, yy, psf_image)


                    surface = surface_fit(xx,yy)

                    # bkg_median = np.nanmedian(surface)
                    psf_image_bkg_free = psf_image - surface


                if syntax['psf_bkg_local']:

                    pos = list(zip([psf_image.shape[0]/2],[psf_image.shape[0]/2]))

                    ap , bkg = ap_phot(pos ,
                                       image,
                                       radius = syntax['ap_size'] * fwhm,
                                       r_in   = syntax['r_in_size'] * fwhm,
                                       r_out  = syntax['r_out_size'] * fwhm)

                    psf_image_bkg_free = psf_image - ( np.ones(psf_image.shape) * bkg)

                    # bkg_median = bkg[0]

                x = np.arange(0,2*syntax['scale'])
                xx,yy= np.meshgrid(x,x)

                pars = lmfit.Parameters()

                pars.add('A',value = np.nanmax(psf_image_bkg_free),min = 0)
                pars.add('x0',value = psf_image_bkg_free.shape[0]/2,min = 0,max = psf_image_bkg_free.shape[0])
                pars.add('y0',value = psf_image_bkg_free.shape[0]/2,min = 0,max = psf_image_bkg_free.shape[0])
                pars.add('sigma',value = sigma )
                pars.add('sky',value = np.nanmedian(psf_image_bkg_free))

                def residual(p):
                    p = p.valuesdict()
                    return (psf_image_bkg_free - gauss_2d((xx,yy),p['x0'],p['y0'],p['sky'],p['A'],p['sigma']).reshape(psf_image_bkg_free.shape)).flatten()

                mini = lmfit.Minimizer(residual, pars)
                result = mini.minimize(method = 'leastsq')
                xc = result.params['x0'].value
                yc = result.params['y0'].value

                ap_range = np.arange(0.1,syntax['scale']/fwhm,1/25)
                ap_sum = []

                for nm in ap_range:

                    ap,bkg = ap_phot( [(xc,yc)] ,
                                 psf_image_bkg_free,
                                 radius = nm * fwhm,
                                 r_in = syntax['r_in_size'] * fwhm,
                                 r_out = syntax['r_out_size'] * fwhm)

                    ap_sum.append(ap)

                ap_sum = ap_sum/np.nanmax(ap_sum)

                radius = ap_range[np.argmax(ap_sum>=syntax['norm_count_sum'])]

                image_radius = radius * fwhm
                image_radius_lst.append(image_radius)

                # global pixel coorindates base on bn gaussian fit
                xc_global = xc - syntax['scale'] + int(selected_sources.x_pix[idx])
                yc_global = yc - syntax['scale'] + int(selected_sources.y_pix[idx])

                # recenter image absed on location of best fit x and y
                psf_image = image[int(yc_global)-syntax['scale']: int(yc_global) + syntax['scale'],
                                  int(xc_global)-syntax['scale']: int(xc_global) + syntax['scale']]

                if syntax['psf_bkg_surface']:


                    sigma_clip = SigmaClip(sigma=syntax['Lim_SNR'])
                    bkg_estimator = MedianBackground()
                    bkg = Background2D(psf_image, (3, 3),
                                       filter_size=(5, 5),
                                       sigma_clip=sigma_clip,
                                       bkg_estimator=bkg_estimator)


                    surface = bkg.background

                    # bkg_median = np.nanmedian(bkg.background_median)
                    psf_image_bkg_free = psf_image- surface

                if syntax['psf_bkg_poly']:


                    # Fit 2D surface to close_up
                    from astropy.modeling import models, fitting
                    surface_function_init = models.Polynomial2D(degree=syntax['psf_bkg_poly_degree'])

                    fit_surface = fitting.LevMarLSQFitter()


                    x = np.arange(0,psf_image.shape[0])
                    y = np.arange(0,psf_image.shape[1])
                    xx,yy= np.meshgrid(x,y)

                    with warnings.catch_warnings():
                        # Ignore model linearity warning from the fitter
                        warnings.simplefilter('ignore')

                        surface_fit = fit_surface(surface_function_init, xx, yy, psf_image)


                        surface = surface_fit(xx,yy)

                        # bkg_median = np.nanmedian(surface)
                        psf_image_bkg_free = psf_image - surface


                if syntax['psf_bkg_local']:

                    pos = list(zip([psf_image.shape[0]/2],[psf_image.shape[0]/2]))

                    ap , bkg = ap_phot(pos ,
                                       image,
                                       radius = syntax['ap_size'] * fwhm,
                                       r_in   = syntax['r_in_size'] * fwhm,
                                       r_out  = syntax['r_out_size'] * fwhm)

                    psf_image_bkg_free = psf_image - ( np.ones(psf_image.shape) * bkg)

                    # bkg_median = bkg[0]

                psf_image_slice = psf_image_bkg_free[int(psf_image_bkg_free.shape[0]/2 - fitting_radius):int(psf_image_bkg_free.shape[0]/2 + fitting_radius) ,
                                            int(psf_image_bkg_free.shape[0]/2 - fitting_radius):int(psf_image_bkg_free.shape[0]/2 + fitting_radius) ]


                x_slice = np.arange(0,2*fitting_radius)
                xx_sl,yy_sl= np.meshgrid(x_slice,x_slice)

                pars = lmfit.Parameters()
                pars.add('A',value = np.nanmax(psf_image_slice),min = 1e-5)
                pars.add('x0',value = psf_image_slice.shape[0]/2)
                pars.add('y0',value = psf_image_slice.shape[0]/2)
                pars.add('sigma',value = sigma )

                def residual(p):
                    p = p.valuesdict()
                    return (psf_image_slice - gauss_2d((xx_sl,yy_sl),p['x0'],p['y0'],0,p['A'],p['sigma']).reshape(psf_image_slice.shape)).flatten()

                mini = lmfit.Minimizer(residual, pars,nan_policy = 'omit')

                result = mini.minimize(method = 'leastsq')

                positions  = list(zip([xc_global ],[yc_global ]))

                psf_counts,psf_bkg = ap_phot(positions,
                                             image,
                                             radius = syntax['ap_size']    * fwhm,
                                             r_in   = syntax['r_in_size']  * fwhm,
                                             r_out  = syntax['r_out_size'] * fwhm)

                psf_SNR = SNR(psf_counts,psf_bkg,syntax['exp_time'],0,syntax['ap_size']* fwhm,syntax['gain'],0)[0]

                if psf_SNR < syntax['construction_SNR']:
                    logger.debug('PSF constuction source too low: %s' % int(psf_SNR))
                    continue
                else:
                    logger.debug('PSF construction source SNR: %s' % int(psf_SNR))
                    # print('\rPSF source %d / %d :: SNR: %d' % (int(psf_SNR)),end = '')
                    pass

                xc = result.params['x0'].value
                yc = result.params['y0'].value

                H = result.params['A'].value
                sigma = abs(result.params['sigma'].value)

                xc_os = ( xc - fitting_radius ) + syntax['scale']
                yc_os = ( yc - fitting_radius ) + syntax['scale']

                residual = psf_image_bkg_free - gauss_2d((xx,yy),xc_os,yc_os,0,H,sigma).reshape(psf_image_bkg_free.shape)

                residual /= H

                psf_mag.append(-2.5*np.log10(H))

                residual_regrid = np.repeat(np.repeat(residual, regriding_size, axis=0), regriding_size, axis=1)

                x_roll = scale_roll(fitting_radius,xc,regriding_size)
                y_roll = scale_roll(fitting_radius,yc,regriding_size)

                residual_roll = np.roll(np.roll(residual_regrid,y_roll,axis=0),x_roll,axis = 1)

                if syntax['r_table_shift_check']:
                    fig, ax = plt.subplots(nrows = 2,ncols = 2, figsize = (20,20))

                    def drawArrow(n , A, B):
                        n.arrow(A[0], A[1], B[0] - A[0], B[1] - A[1],
                          head_width=0.5, length_includes_head=True)

                    im = {}

                    im[0] = ax[0,0].imshow(psf_image_bkg_free)
                    ax[0,0].set_title('Source:'+str(n))
                    ax[0,0].axvline(xc_os,linestyle = ':',color = 'red')
                    ax[0,0].axhline(yc_os,label = 'Best fit',linestyle = ':',color = 'red')
                    ax[0,0].axvline(psf_image_bkg_free.shape[0]/2,color = 'black',linestyle = ':')
                    ax[0,0].axhline(psf_image_bkg_free.shape[0]/2,color = 'black',label = 'Center of image',linestyle = ':')


                    im[1] = ax[0,1].imshow(residual)
                    ax[0,1].set_title('Residuals from source @ Normal scale')
                    ax[0,1].axvline(xc_os,linestyle = ':',color = 'red')
                    ax[0,1].axhline(yc_os,label = 'Best fit',linestyle = ':',color = 'red')
                    ax[0,1].axvline(psf_image_bkg_free.shape[0]/2,color = 'black',linestyle = ':')
                    ax[0,1].axhline(psf_image_bkg_free.shape[0]/2,color = 'black',label = 'Center of image',linestyle = ':')

                    im[2] = ax[1,0].imshow(residual_regrid)
                    ax[1,0].set_title('Residuals  shift s.t xc,yc are at center')
                    ax[1,0].scatter(array_correction(xc_os*m),array_correction(yc_os*m),marker = 'x',color = 'red',label = 'Location of best fit before roll with regridding')
                    ax[1,0].scatter(residual_regrid.shape[0]/2,residual_regrid.shape[0]/2,marker = 'o',facecolors='none', edgecolors='black',label = 'Image center')
                    ax[1,0].set_xticks(np.arange(-.5, residual_regrid.shape[0], 1), minor=True)
                    ax[1,0].set_yticks(np.arange(-.5, residual_regrid.shape[0], 1), minor=True)
                    ax[1,0].grid(which='minor', color='black', linestyle='-', linewidth = 1 , alpha = 0.1)

                    im[3] = ax[1,1].imshow(residual_roll)
                    ax[1,1].set_title('Regrid Residuals after roll')
                    ax[1,1].scatter(array_correction(x_roll +xc_os*m),array_correction(y_roll +yc_os*m),marker = 'x',color = 'red',label = 'Location of best fit plus roll adjustment ')
                    ax[1,1].scatter(residual_regrid.shape[0]/2,residual_regrid.shape[0]/2,marker = 'o',facecolors='none', edgecolors='black',label = 'Image center')
                    ax[1,1].plot([None],[None],label = 'Roll: xr '+str(x_roll)+' yr '+str(y_roll))
                    ax[1,1].set_xticks(np.arange(-.5, residual_regrid.shape[0], 1), minor=True);
                    ax[1,1].set_yticks(np.arange(-.5, residual_regrid.shape[0], 1), minor=True)
                    ax[1,1].grid(which='minor', color='black', linestyle='-', linewidth = 1 , alpha = 0.1)


                    from mpl_toolkits.axes_grid1 import make_axes_locatable

                    for i in range(len(ax.reshape(-1))):
                        na = ax.reshape(-1)[i]

                        na.legend(loc = 'best')

                        divider = make_axes_locatable(na)
                        cax = divider.append_axes("right", size="5%", pad=0.05)
                        fig.colorbar(im[i], cax=cax)

                    pathlib.Path(syntax['write_dir']+'residual_shift_check/').mkdir(parents = True, exist_ok=True)

                    import time

                    if os.path.exists(syntax['write_dir']+'residual_shift_check/'+'result.png'):
                        plt.savefig(syntax['write_dir']+'residual_shift_check/'+'result_{}.png'.format(int(time.time())))
                    else:
                        plt.savefig(syntax['write_dir']+'residual_shift_check/'+'result.png')


                    plt.close(fig)

                residual_table += residual_roll
                construction_sources.append([xc_global,yc_global,H/syntax['exp_time']])
                logger.debug('Residual table updated: %d / %d ' % (sources_used,syntax['psf_source_no']))

                print('\rResidual table updated: %d / %d ' % (sources_used,syntax['psf_source_no']) ,end = '')
                sources_used +=1

                sigma_fit.append(sigma)

            except Exception as e:
                logger.exception(e)

                logger.error('trying another source')
                failsafe+=1
                n+=1
                try:
                    plt.close(fig)
                except:
                    pass

                continue
        print('')

        if sources_used < syntax['min_psf_source_no']:
            logger.warning('BUILDING PSF: Not enough useable sources found')
            return None,None,construction_sources.append([np.nan]*5),syntax

        logger.debug('PSF Successful')
        residual_table/= sources_used

        residual_table  = rebin(residual_table,(2*syntax['scale'],2*syntax['scale']))

        construction_sources = pd.DataFrame(construction_sources)


        '''
        Black magic time -- Removing negative counts in PSF model and updating residual table

        '''
        x_rebin = np.arange(0,residual_table.shape[0])

        xx_rebin,yy_rebin = np.meshgrid(x_rebin,x_rebin)
        psf_sigma = np.mean(sigma_fit)

        psf_core = gauss_2d((xx_rebin,yy_rebin),residual_table.shape[0]/2,residual_table.shape[0]/2,0,1,psf_sigma).reshape(residual_table.shape)

        tmp_psf = psf_core + residual_table

        tmp_psf[tmp_psf > 0] = 0

        residual_table -= tmp_psf

        image_radius_lst = np.array(image_radius_lst)

        syntax['image_radius'] = image_radius_lst.mean()
        logger.info('Image_radius [pix] : %s +/- %s' % (round(image_radius_lst.mean(),3) ,round(image_radius_lst.std(),3)))

        # for finding error on psf, not implemented yet
        construction_sources.columns = ['x_pix','y_pix','H_psf']
        construction_sources['H_psf_err'] = [0]*len(construction_sources)

    except Exception as e:
                logger.exception('BUILDING PSF: ',e)
                raise Exception



    return residual_table,sigma_fit,construction_sources,syntax

'''
Once PSF is built, Do fitting
'''

def fit(image,
        sources,
        residual_table,
        syntax,
        fwhm,
        save_plot = False,
        show_plot = False,
        rm_bkg_val = True,
        hold_pos = False,
        return_fwhm = False,
        return_subtraction_image = False,
        fname = None

        ):

    import numpy as np
    import pandas as pd
    import pathlib
    import lmfit
    import logging
    # import sys,os
    import matplotlib.pyplot as plt
    from astropy.stats import SigmaClip
    from photutils import Background2D, MedianBackground




    from autophot.packages.aperture  import ap_phot

    logger = logging.getLogger(__name__)


    fitting_radius = int(np.ceil(1.5*fwhm))
    regriding_size = int(syntax['regrid_size'])

    sources = sources
    residual_table = residual_table

    def build_psf(xc, yc, sky, H, r_table, slice_scale = None,pad_shape = None):

        '''
        Slice_scale = only return "slice scaled" image
        '''

        try:

            psf_shape = r_table.shape


            if pad_shape != None:
                psf_shape = pad_shape


            if pad_shape != None:
                top =    int((pad_shape[0] - r_table.shape[0])/2)
                bottom=  int((pad_shape[0] - r_table.shape[0])/2)
                left =   int((pad_shape[1] - r_table.shape[1])/2)
                right =  int((pad_shape[1] - r_table.shape[1])/2)

                r_table = np.pad(r_table, [(top, bottom), (left, right)], mode='constant', constant_values=0)

            x_rebin = np.arange(0,psf_shape[0])
            y_rebin = np.arange(0,psf_shape[1])

            xx_rebin,yy_rebin = np.meshgrid(x_rebin,y_rebin)

            sigma  = fwhm/(2*np.sqrt(2*np.log(2)))

            core = gauss_2d((xx_rebin,yy_rebin),xc,yc,sky,H,sigma).reshape(psf_shape)

            residual_rebinned = np.repeat(np.repeat(r_table, regriding_size, axis=0), regriding_size, axis=1)

            x_roll = scale_roll(xc,int(r_table.shape[0]/2),regriding_size)
            y_roll = scale_roll(yc,int(r_table.shape[1]/2),regriding_size)

            residual_roll = np.roll(np.roll(residual_rebinned,y_roll,axis=0),x_roll,axis = 1)

            residual = rebin(residual_roll,psf_shape)

            psf =  (sky  + (H * residual)) + core

            if np.isnan(np.min(psf)):
                logger.info(sky,H,np.min(residual),np.min(core))

            psf[np.isnan(psf)] = 0
            psf[psf<0] = 0

            if slice_scale != None:
                psf = psf[int ( 0.5 * r_table.shape[0] - slice_scale): int(0.5*r_table.shape[0] + slice_scale),
                          int ( 0.5 * r_table.shape[0] - slice_scale): int(0.5*r_table.shape[0] + slice_scale)]

        except Exception as e:
            logger.exception(e)
            psf = np.nan

        return psf


    psf_params = []

    x = np.arange(0,2*syntax['scale'])

    xx,yy= np.meshgrid(x,x)

    if hold_pos:
        dx = 1e-3
        dy = 1e-3
    else:
        dx = syntax['dx']
        dy = syntax['dy']


    lower_x_bound = syntax['scale']
    lower_y_bound = syntax['scale']
    upper_x_bound = syntax['scale']
    upper_y_bound = syntax['scale']


    if return_subtraction_image:

        from astropy.visualization.mpl_normalize import ImageNormalize

        from astropy.visualization import  ZScaleInterval, SquaredStretch


        norm = ImageNormalize( stretch = SquaredStretch())

        vmin,vmax = (ZScaleInterval(nsamples = 1500)).get_limits(image)

    '''
    Known issue - for poor images, some sources may be too close to boundary, remove this
    '''
    if not return_fwhm:
        logger.info('Image cutout size: (%.f,%.f) (%.f,%.f)' % ((lower_x_bound,upper_x_bound,lower_y_bound,upper_y_bound)))

    sources = sources[sources.x_pix < image.shape[1] - upper_x_bound]
    sources = sources[sources.x_pix > lower_x_bound]
    sources = sources[sources.y_pix < image.shape[0] - upper_y_bound]
    sources = sources[sources.y_pix > lower_y_bound]

    logger.info('Fitting PSF to %d sources' % len(sources))

    for n  in range(len(sources.index)):
        if return_fwhm:
            print('\rFitting PSF to source: %d / %d' % (n+1,len(sources)), end = '',flush=False)
        try:

            idx = list(sources.index)[n]


            source =   image[int(sources.y_pix[idx])-lower_y_bound: int(sources.y_pix[idx]) + upper_y_bound,
                             int(sources.x_pix[idx])-lower_x_bound: int(sources.x_pix[idx]) + upper_x_bound]

            if source.shape != (int(2*syntax['scale']),int(2*syntax['scale'])):

                bkg_median = np.nan
                H = np.nan
                H_psf_err = np.nan
                psf_params.append((idx,bkg_median,H,H_psf_err))
                continue

            xc = syntax['scale']
            yc = syntax['scale']

            xc_global =  sources.x_pix[idx]
            yc_global =  sources.y_pix[idx]

            if not rm_bkg_val:
                source_bkg_free = source
                bkg_median = 0
            else:

                try:
                    if syntax['psf_bkg_surface']:

                        sigma_clip = SigmaClip(sigma=syntax['Lim_SNR'])
                        bkg_estimator = MedianBackground()
                        bkg = Background2D(source, (3, 3), filter_size=(5, 5),
                                           sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)


                        surface = bkg.background

                        bkg_median = np.nanmedian(bkg.background_median)
                        source_bkg_free = source - surface

                    if syntax['psf_bkg_poly']:


                        # Fit 2D surface to close_up
                        from astropy.modeling import models, fitting
                        surface_function_init = models.Polynomial2D(degree=syntax['psf_bkg_poly_degree'])

                        fit_surface = fitting.LevMarLSQFitter()


                        x = np.arange(0,source.shape[0])
                        y = np.arange(0,source.shape[0])
                        xx,yy= np.meshgrid(x,y)

                        with warnings.catch_warnings():
                            # Ignore model linearity warning from the fitter
                            warnings.simplefilter('ignore')
                            surface_fit = fit_surface(surface_function_init, xx, yy, source)

                        surface = surface_fit(xx,yy)

                        bkg_median = np.nanmedian(surface)
                        source_bkg_free = source - surface


                    if syntax['psf_bkg_local']:

                        pos = list(zip([xc_global],[yc_global]))

                        ap , bkg = ap_phot(pos ,
                                           image,
                                           radius = syntax['ap_size'] * fwhm,
                                           r_in   = syntax['r_in_size'] * fwhm,
                                           r_out  = syntax['r_out_size'] * fwhm)

                        source_bkg_free = source - ( np.ones(source.shape) * bkg)

                        bkg_median = bkg[0]

                except Exception as e:
                    logger.exception(e)


            source = source_bkg_free[int(source_bkg_free.shape[1]/2 - fitting_radius):int(source_bkg_free.shape[1]/2 + fitting_radius) ,
                                     int(source_bkg_free.shape[0]/2 - fitting_radius):int(source_bkg_free.shape[0]/2 + fitting_radius) ]

            if source.shape != (int(2*fitting_radius),int(2*fitting_radius)):
                bkg_median = np.nan
                H = np.nan
                H_psf_err = np.nan
                psf_params.append((idx,bkg_median,H,H_psf_err))
                continue

            if np.sum(np.isnan(source)) == len(source):
                logger.warning('all nan image')
                continue


            if hold_pos:
                dx = 1e-3
                dy = 1e-3
            else:
                dx = syntax['dx']
                dy = syntax['dy']

            x_slice = np.arange(0,2*fitting_radius)
            xx_sl,yy_sl= np.meshgrid(x_slice,x_slice)

            if return_fwhm:
                logger.info('Fitting gaussian to source to get FWHM')

                pars = lmfit.Parameters()
                pars.add('A',value = np.nanmax(source),min = 1e-5)
                pars.add('x0',value = source.shape[0]/2)
                pars.add('y0',value = source.shape[0]/2)
                pars.add('sigma',value = fwhm / 2*np.sqrt(2 * np.log(2)) )

                def residual(p):
                    p = p.valuesdict()
                    return (source - gauss_2d((xx_sl,yy_sl),p['x0'],p['y0'],0,p['A'],p['sigma']).reshape(source.shape)).flatten()

                mini = lmfit.Minimizer(residual, pars,nan_policy = 'omit')

                result = mini.minimize(method = 'leastsq')

                source_fwhm = result.params['sigma'].value  * 2*np.sqrt(2 * np.log(2))

                logger.info('Target FWHM: %.3f' % source_fwhm)

            try:

                pars = lmfit.Parameters()

                pars.add('A', value = np.nanmax(source)*0.75,min = 1e-5)
                pars.add('x0',value = 0.5*residual_table.shape[0],min = 0.5*residual_table.shape[0]-dx,max = 0.5*residual_table.shape[0]+dx)
                pars.add('y0',value = 0.5*residual_table.shape[1],min = 0.5*residual_table.shape[1]-dy,max = 0.5*residual_table.shape[1]+dy)

                def residual(p):
                    p = p.valuesdict()
                    res = ((source - build_psf(p['x0'],p['y0'],0,p['A'],residual_table,slice_scale = source.shape[0]/2)))
                    return res.flatten()

                mini = lmfit.Minimizer(residual,pars,nan_policy = 'omit',scale_covar=True)

                result = mini.minimize(method = 'leastsq')

                xc = result.params['x0'].value

                yc = result.params['y0'].value

                H = result.params['A'].value

                H_psf_err = result.params['A'].stderr



            except Exception as e:
                logger.exception(e)
                continue

            best = 2
            cf68 = [1,3]

            # This is temporary
            H_psf_err = 0

            '''
            Work in progress getting error output to work
            '''

            if syntax['use_confinterval']:
                try:
                    error_params = lmfit.conf_interval(mini, result,sigmas=[1, 2],verbose=False)

                    H_psf_err_low = error_params['A'][cf68[0]][1] - error_params['A'][best][1]
                    H_psf_err_high = error_params['A'][best][1] - error_params['A'][cf68[1]][1]
                    H_psf_err = [H_psf_err_low,H_psf_err_high]

                except:
                    syntax['use_covarience'] = True
                    pass

            if syntax['use_covarience']:

                H_psf_err = result.params['A'].stderr

            else:
                logger.warning('Error not computed')
                H_psf_err = 0

            psf_params.append((idx,bkg_median,H,H_psf_err))

            if return_subtraction_image:

                try:

                    image_section = image[int(yc_global -  syntax['scale']): int(yc_global + syntax['scale']),
                                          int(xc_global -  syntax['scale']): int(xc_global + syntax['scale'])].copy()

                    image[int(yc_global  - syntax['scale']): int(yc_global +  syntax['scale']),
                          int(xc_global  - syntax['scale']): int(xc_global +  syntax['scale'])] =  image_section - build_psf(xc , yc, 0, H, residual_table)

                    image_section_subtraction = image_section - build_psf(xc , yc, 0, H, residual_table)

                    fig, ax1, = plt.subplots()

                    ax_before = ax1.inset_axes([0.95, 0.70, 0.4, 0.25])
                    ax_after  = ax1.inset_axes([0.95, 0.20, 0.4, 0.25])

                    ax1.imshow(image,
                               vmin = vmin,
                               vmax = vmax,
                               norm = norm,
                               origin = 'lower',
                               cmap = 'gist_heat',
                               interpolation = 'nearest')

                    ax1.set_xlim(0,image.shape[0])
                    ax1.set_ylim(0,image.shape[1])

                    ax1.scatter(xc_global,
                                yc_global,
                                marker = 'o',
                                facecolor = 'None',
                                color = 'green',
                                s = 25)


                    ax_before.imshow(image_section,
                                      vmin = vmin,
                                      vmax = vmax,
                                      norm = norm,
                                     origin = 'lower',
                                     cmap = 'gist_heat',
                                     interpolation = 'nearest')

                    ax_after.imshow(image_section_subtraction,
                                     vmin = vmin,
                                     vmax = vmax,
                                     norm = norm,
                                    origin = 'lower',
                                    cmap = 'gist_heat',
                                    interpolation = 'nearest')

                    ax_after.axis('off')
                    ax_before.axis('off')
                    ax1.axis('off')

                    ax_after.set_title('After')
                    ax_before.set_title('Before')

                    logger.info('Image %s / %s saved' % (str(idx),str(len(sources.index))))
                    plt.close(fig)

                except Exception as e:
                    logger.exception(e)
                    plt.close('all')
                    pass





            if syntax['show_residuals'] or show_plot == True or save_plot == True:
                fig_source_residual = plt.figure(figsize = (14,8),facecolor='w', edgecolor='k')
                gs = fig_source_residual.add_gridspec(2, 5)
                try:

                    from matplotlib.patches import Circle

                    from mpl_toolkits.axes_grid1 import make_axes_locatable
                    from mpl_toolkits import mplot3d


                    ax1 = fig_source_residual.add_subplot(gs[0, 0])
                    ax2 = fig_source_residual.add_subplot(gs[1,0], sharex=ax1)
                    ax3 = fig_source_residual.add_subplot(gs[:,-4:-2], projection='3d')
                    ax4 = fig_source_residual.add_subplot(gs[:,-2:], projection='3d')

                    ax1.imshow(source_bkg_free)
                    ax1.scatter(xc,yc,label = 'Best fit',
                                marker = '*',
                                color = 'red',
                                alpha = 0.5)

                    ax2.imshow(source_bkg_free - build_psf(xc,yc,0,H,residual_table))
                    ax2.scatter(xc,yc,label = 'Best fit',
                                marker = '*',
                                color = 'red',
                                alpha = 0.5)


                    x  = np.linspace(0, 2*syntax['scale'], 2*syntax['scale'])
                    y  = np.linspace(0, 2*syntax['scale'], 2*syntax['scale'])

                    X, Y = np.meshgrid(x, y)

                    ax3.plot_surface(X, Y,source_bkg_free,
                                     cmap='viridis',rstride=1, cstride=1)


                    ax4.plot_surface(X, Y,source_bkg_free - build_psf( xc,yc,0,H ,residual_table),
                                     cmap='viridis',
                                     rstride=1,
                                     cstride=1,
                                     edgecolor='none')

                    ax3.set_zlabel("Amplitude")
                    ax4.set_zlabel("Amplitude")

                    ax1.text(1.1, 0.5,'Image',
                             horizontalalignment='center',
                             verticalalignment='center',
                             rotation = 270,
                             transform = ax1.transAxes)

                    ax2.text(1.1, 0.5,'Residual Image',
                             horizontalalignment='center',
                             verticalalignment='center',
                             rotation = 270,
                             transform = ax2.transAxes)

                    for ax in [ax1,ax2,ax3,ax4]:
                        ax.set_xlabel('X pixel')
                        ax.set_ylabel('Y pixel')

                    plt.tight_layout()

                    if save_plot == True:
                        fig_source_residual.savefig(syntax['write_dir']+'target_psf_'+fname+'.pdf',
                                                    format = 'pdf')

                        logger.info('Image %s / %s saved' % (str(n+1),str(len(sources.index)) ))
                    else:

                        pathlib.Path(syntax['write_dir']+'/'+'psf_subtractions/').mkdir(parents = True, exist_ok=True)

                        plt.savefig(syntax['write_dir']+'psf_subtractions/'+'psf_subtraction_{}.png'.format(int(n)))

                    plt.close(fig_source_residual)

                except Exception as e:
                    logger.exception(e)
                    plt.close('all')


        except Exception as e:
             logger.exception(e)


             bkg_median = np.nan
             H = np.nan
             H_psf_err = np.nan
             psf_params.append((idx,bkg_median,H,H_psf_err))
             plt.close(fig_source_residual)
             continue

    new_df =  pd.DataFrame(psf_params,columns = ('idx','bkg','H_psf','H_psf_err'),index = sources.index)

    if return_fwhm:
        new_df['target_fwhm'] = round(source_fwhm,3)
    print('')

    return pd.concat([sources,new_df],axis = 1),build_psf






'''
Perform PSF calculations
'''

def do(df,residual,syntax,fwhm):

    try:
        from photutils import CircularAperture
        from photutils import aperture_photometry
        from scipy.integrate import dblquad
        import sys,os
        import logging


        logger = logging.getLogger(__name__)

        xc = syntax['scale']
        yc = syntax['scale']

        # Integration radius
        int_scale = syntax['image_radius']
#        int_scale = syntax['ap_size'] * fwhm

        int_range_x = [xc - int_scale , xc + int_scale]
        int_range_y = [yc - int_scale , yc + int_scale]

        sigma  = fwhm/(2*np.sqrt(2*np.log(2)))

        # Core Gaussian component with height 1 and sigma value sigma
        core= lambda y, x: gauss_2d((x,y),syntax['scale'],syntax['scale'],0,1,sigma)
        core_int = dblquad(core, int_range_y[0],int_range_y[1],lambda x:int_range_x[0],lambda x:int_range_x[1])[0]

        # Aperture Photometry over residual
        apertures = CircularAperture((syntax['scale'],syntax['scale']), r=int_scale)
        phot_table = aperture_photometry(residual, apertures,method='subpixel',subpixels=5)

        phot_table['aperture_sum'].info.format = '%.8g'
        residual_int = phot_table[0]

        # Counts from core compoent on PSF
        syntax['c_counts'] = float(core_int)

        # Counts from Residual component of PSF
        syntax['r_counts'] = float(residual_int['aperture_sum'])

        # Counts in a PSF with fwhm 2 sqrt(2 ln 2) * sigma and height 1
        sudo_psf = core_int+float(residual_int['aperture_sum'])

        psf_int     = df.H_psf     * sudo_psf
        psf_int_err = df.H_psf_err * sudo_psf

        df['psf_counts'] = psf_int.values
        df['psf_counts_err'] = psf_int_err.values

    except Exception as e:
        exc_type, exc_oba, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        logger.info(exc_type, fname, exc_tb.tb_lineno,e)
        df = np.nan

    return df,syntax

