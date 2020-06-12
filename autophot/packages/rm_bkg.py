#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def subtract_bkg(source,
                 syntax,
                 xc=0,
                 yc=0):

    import logging
    import warnings
    import numpy as np
    from astropy.stats import SigmaClip
    from photutils import Background2D, MedianBackground
    from autophot.packages.aperture import ap_phot
    from astropy.modeling import models, fitting

    logger = logging.getLogger(__name__)

    try:
        if syntax['psf_bkg_surface']:

            sigma_clip = SigmaClip(sigma=syntax['Lim_SNR'])
            bkg_estimator = MedianBackground()
            bkg = Background2D(source,
                               (3, 3),
                               filter_size=(5, 5),
                               sigma_clip=sigma_clip,
                               bkg_estimator=bkg_estimator)


            surface = bkg.background

            bkg_median = np.nanmedian(bkg.background_median)

            source_bkg_free = source - surface

        if syntax['psf_bkg_poly']:


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

            pos = list(zip([xc],[yc]))

            ap,bkg = ap_phot(pos,
                             source,
                             radius = syntax['ap_size'] * syntax['fwhm'],
                             r_in   = syntax['r_in_size'] * syntax['fwhm'],
                             r_out  = syntax['r_out_size'] * syntax['fwhm'])

            source_bkg_free = source - ( np.ones(source.shape) * bkg)

            bkg_median = bkg[0]

    except Exception as e:
        logger.exception(e)

    return source_bkg_free,bkg_median

