#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 16:49:01 2018

@author: seanbrennan
"""

from photutils import CircularAperture
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground,BkgZoomInterpolator



def rm_surface(close_up,syntax,center = None,radius = None):
    # Add in fwhm info here!
    
    if center == None:
        center = [syntax['scale'],syntax['scale']]
    
    if radius == None:
        radius = syntax['fwhm_guess']

    
    aperture = CircularAperture((center[0],center[1]),
                                r=radius)

    mask = aperture.to_mask(method='center')
    

    
#    mask = masks[0]
    
    image = mask.to_image(shape=((close_up.shape[1],close_up.shape[0])))
    
    
    sigma_clip = SigmaClip(sigma=3., maxiters=5)
    
    bkg_estimator = MedianBackground()
    
    bkg_interpolator = BkgZoomInterpolator()
    
    bkg1 = Background2D(close_up,
                        box_size =  (2, 2), 
                        filter_size=(2, 2),
                        sigma_clip=sigma_clip,
                        bkg_estimator=bkg_estimator,
                        mask = image,
                        interpolator = bkg_interpolator,
                        )
    back3 = bkg1.background
    
    return close_up-back3, back3, bkg1