#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 14:10:14 2020

@author: seanbrennan
"""

'''

AutoPHoT Limiting Magnitude Module
'''

import matplotlib.pyplot as plt

plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
plt.rcParams['font.size'] = 7
plt.rcParams['figure.subplot.hspace']= 0.2
plt.rcParams['figure.subplot.wspace']= 0.2
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
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['xtick.minor.width'] = 0.35
plt.rcParams['xtick.direction'] = 'in'

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
            height meanltiple of width

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

def gauss_2d(image, x0,y0, sky , A, sigma):
    import numpy as np
    (x,y) = image
    a = (x-x0)**2
    b = (y-y0)**2
    c = (2*sigma**2)
    d =  A*np.exp( -(a+b)/c)
    e =  d + sky
    return  e

def limiting_magnitude_prob(syntax,image,model =None,r_table=None):

    '''
    syntax - dict
        dDtionary of input paramters
    image - np.array
        Image of region of interest with target in center of image
    model - function
        - psf function from autophot
    '''
    try:

        from photutils import CircularAperture
        import matplotlib.pyplot as plt
        import numpy as np
        import matplotlib.gridspec as gridspec
        import random
        from astropy.stats import SigmaClip
        from photutils import Background2D, MedianBackground
        from scipy.optimize import curve_fit
        import warnings

        from autophot.packages.functions import r_dist

        import logging

        logger = logging.getLogger(__name__)

        # level for detection - Rule of thumb ~ 5 is a good detection level
        level = syntax['lim_SNR']

        # Lower_level
        low_level = 3

        logger.info('Limiting threshold: %d sigma' % level)

        if syntax['psf_bkg_surface']:



            sigma_clip = SigmaClip(sigma=syntax['lim_SNR'])
            bkg_estimator = MedianBackground()
            bkg = Background2D(image, (3, 3), filter_size=(4, 4),
                               sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)


            surface = bkg.background

            image_no_surface = image - surface

        if syntax['psf_bkg_poly']:

            from astropy.modeling import models, fitting
            surface_function_init = models.Polynomial2D(degree=syntax['psf_bkg_poly_degree'])

            fit_surface = fitting.LevMarLSQFitter()

            x = np.arange(0,image.shape[0])
            y = np.arange(0,image.shape[0])
            xx,yy= np.meshgrid(x,y)


            with warnings.catch_warnings():

                # Ignore model linearity warning from the fitter
                warnings.simplefilter('ignore')
                surface_fit = fit_surface(surface_function_init, xx, yy,image)


            surface = surface_fit(xx,yy)
            image_no_surface = image - surface

        if syntax['psf_bkg_local']:

            surface = np.ones(image.shape) * np.nanmedian(image)
            image_no_surface = image - np.nanmedian(image)


        # "size" of source
        source_size = 1.3 * syntax['image_radius']

        # Mask out target region
        mask_ap  = CircularAperture([image.shape[0]/2,image.shape[1]/2],r = source_size)

        mask = mask_ap.to_mask(method='center')
        logging.info('Number of pixels in star: %d' % np.sum(mask.to_image(image.shape)))

        # Mask out center region
        mask_image  = (image_no_surface ) * (1-mask.to_image(image.shape))

        number_of_points = 1000

        fake_points = {}
        for i in range(number_of_points):
            fake_points[i] = []
            for j in range(int(np.sum(mask.to_image(image.shape)))):

                xp_ran = random.randint(0,int(image.shape[0]-1))
                yp_ran = random.randint(0,int(image.shape[1]-1))

                rp = r_dist(image.shape[0]/2,xp_ran,image.shape[0]/2,yp_ran)

                if rp < source_size:

                    r_ran = random.randint(int(source_size  - rp),int(image.shape[0]/2 - rp - 1 ))
                    theta_ran = random.uniform(0,2*np.pi)
                    xp_ran = r_ran*np.cos(theta_ran) + xp_ran
                    yp_ran = r_ran*np.sin(theta_ran) + yp_ran
                fake_points[i].append((int(xp_ran),int(yp_ran)))

        fake_sum = {}
        for i in range(len(fake_points.keys())):
            fake_sum[i] = []
            list_tmp = []
            for j in fake_points[i]:

                list_tmp.append(mask_image[j[0]][j[1]])

            fake_sum[i] = list_tmp

        fake_mags = {}

        for f in fake_sum.keys():

            fake_mags[f] = np.sum(fake_sum[f])

        hist, bins = np.histogram(list(fake_mags.values()),bins = len(list(fake_mags.values())),density = True)
        center = (bins[:-1] + bins[1:]) / 2

        sigma = np.nanstd(list(fake_mags.values()))
        mean = np.nanmean(list(fake_mags.values()))
        A = np.nanmax(hist)

        def gauss(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))

        popt,pcov = curve_fit(gauss,center,hist,p0=[A,mean,sigma],absolute_sigma=True )

        mean = popt[1]
        std  = abs(popt[2])

        logging.info('Mean: %s - std: %s' % (round(mean,3),round(std,3)))
        logging.info('Detection at %s std: %s' % (level,round(mean + level*std,3)))

        limiting_mag_figure = plt.figure(figsize = set_size(540,aspect = 1))
        gs = gridspec.GridSpec(2, 2)
        ax0 = limiting_mag_figure.add_subplot(gs[:, :-1])

        ax1 = limiting_mag_figure.add_subplot(gs[-1, -1])
        ax2 = limiting_mag_figure.add_subplot(gs[:-1, -1])

        line_kwargs = dict(alpha=0.5,color='black',ls = '--')

        # the histogram of the data
        n, bins, patches = ax0.hist(list(fake_mags.values()),
                                    density=True,
                                    bins = 30,
                                    facecolor='blue',
                                    alpha=1,
                                    label = 'Pseudo-Flux Distribution')

        ax0.axvline(mean,**line_kwargs)
        ax0.axvline(mean + 1*std,**line_kwargs)
        ax0.text(mean + 1*std,np.max(n),r'$1\sigma$')
        ax0.axvline(mean + 2*std,**line_kwargs)
        ax0.text(mean + 2*std,np.max(n),r'$2\sigma$')
        ax0.axvline(mean +  level*std,**line_kwargs)
        ax0.text(mean + level*std,np.max(n),r'$'+str(level)+r'\sigma$')

        x_fit = np.linspace(ax0.get_xlim()[0], ax0.get_xlim()[1], 1000)

        ax0.plot(x_fit, gauss(x_fit,*popt),label = 'Gaussian Fit',color = 'red')

        ax0.ticklabel_format(axis='y', style='sci',scilimits = (-2,0))
        ax0.yaxis.major.formatter._useMathText = True

        ax0.set_xlabel('Pseudo-Flux')
        ax0.set_ylabel('Normalised Probability')

        im2 = ax2.imshow(image-surface,origin='lower')
        plt.colorbar(im2,ax=ax2)
        ax2.set_title('Image - Surface')


        x = random.sample(range(0,int(image.shape[0])), int(image.shape[0])//3)
        y = random.sample(range(0,int(image.shape[1])), int(image.shape[0])//3)

        counts = abs(mean + level*std)

        flux  = counts / syntax['exp_time']

        mag_level = -2.5*np.log10(flux)


        # Worse case scenario

        counts = abs(mean + low_level*std)

        low_flux  = counts / syntax['exp_time']

        low_mag_level = -2.5*np.log10(low_flux)

        try:

            model_label = 'PSF'

            def mag2image(m):
                return (syntax['exp_time']/(syntax['c_counts']+syntax['r_counts']))*(10**(m/-2.5))

            def input_model(x,y,f):

                return model(x, y,0,f,r_table, pad_shape = image.shape,slice_scale = None)

            fake_sources = input_model(x[0], y[0],mag2image(mag_level))
            ax1.scatter(x[0],y[0],marker = 'o',s=150, facecolors='none', edgecolors='r',alpha = 0.5)

        except:
            logging.info('PSF model not available - Using Gaussian')
            model_label = 'Gaussian'
            x_grid = np.arange(0,image.shape[0])
            xx,yy= np.meshgrid(x_grid,x_grid)

            def mag2image(flux):
                sigma = syntax['fwhm'] / 2*np.sqrt(2*np.log(2))
                Amplitude =  flux/np.sqrt(np.pi*sigma)
                return Amplitude

            def input_model(x,y,f):

                return gauss_2d((xx,yy),x,y,0,f,syntax['fwhm'] / 2*np.sqrt(2*np.log(2)))

            fake_sources = input_model(x[0], y[0],mag2image(counts))
            ax1.scatter(x[0],y[0],marker = 'o',s=150, facecolors='none', edgecolors='r',alpha = 0.5)


        try:
            for i in range(1,len(x)):
                fake_sources+= input_model(x[i], y[i],mag2image(mag_level))
                ax1.scatter(x[i],y[i],marker = 'o',s=150, facecolors='none', edgecolors='r',alpha = 0.5)

            im1=ax1.imshow(image - surface + fake_sources, origin='lower')
            ax1.set_title('Fake Sources [%s]' % model_label)

        except Exception as e:
            logging.exception(e)
            im1=ax1.imshow(image - surface , origin='lower')
            ax1.set_title('[ERROR] Fake Sources [%s]' % model_label)


        plt.colorbar(im1,ax=ax1)


        ax0.legend(loc = 'upper left',frameon = False)

        limiting_mag_figure.savefig(syntax['write_dir']+'limiting_mag_porb.pdf',
                                        format = 'pdf')
        plt.close(limiting_mag_figure)

    except Exception as e:
        logging.exception(e)


    return mag_level,low_mag_level