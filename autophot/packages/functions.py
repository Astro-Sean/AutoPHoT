def getheader(fpath):

    '''
    Attempt dynamic  header retrival

    Look for sci image header assuming it is first in list of header
    extensions
    '''

    # Need to rename this function
    from astropy.io.fits import getheader
    from astropy.io import fits

    try:
        try:
            with fits.open(fpath,ignore_missing_end=True) as hdul:
                # Need to verify/fix or can get errors down the line
                hdul.verify('silentfix+ignore')
                headinfo = hdul
            try:
                # try for Telescop keyword
                headinfo['Telescop']
            except:
                raise Exception
        except:
            # try to find 'sci' extension
            headinfo = getheader(fpath,'sci')
    except:
        with fits.open(fpath) as hdul:
            hdul.verify('silentfix+ignore')
            headinfo = hdul

    try:
        # If header is a list
        if isinstance(headinfo,list):
            # If list length contains multiple headers, concat them
            if len(headinfo)>1:
                headinfo_list = headinfo[0].header
                for i in range(1,len(headinfo)):
                    headinfo_list += headinfo[i].header
                headinfo = headinfo_list
            else:
                # is length of list choose this a header
                headinfo = headinfo[0].header
    except Exception as e:
        print(e)

    return headinfo

def getimage(fpath):

    from astropy.io import fits

    try:
        with fits.open(fpath,ignore_missing_end=True) as hdul:
            hdul.verify('silentfix+ignore')
            image = hdul[1].data
    except:
        image = fits.getdata(fpath)

    return image


def zeropoint(mag, counts, ct_gradient = None, dmag = None, airmass = None):

    """
    Calculate zeropint using:

        mag = -2.5 * log10 (counts) + ct(dmag) + zp + airmass
    """
    import numpy as np

    zp_list = [mag]

    mag_inst = -2.5 * np.log10(counts)
    zp_list.append(-1 * mag_inst)


    if type(counts) is np.array:
        counts[np.where(counts < 0.0)] = np.nan

    if ct_gradient != None and all(dmag) != None:
        ct = ct_gradient * dmag
        zp_list.append(-1 * ct)


    zp = sum(zp_list)
    return zp

def mag(counts, zp,ct_gradient = False,dmag = False,airmass = None):

    try:
        """
        Calculate zeropint using:

            mag = - 2.5 * log10 (counts) + ct(dmag) + zp + airmass
        """

        import numpy as np
        import sys
        import os


#        try:
#            counts = np.array(counts)
#    #        counts[np.isnan(counts)] = -1
#    #        counts[counts<0] = np.nan
#        except:
#            pass

        # Iitial list with zeropoint
        if isinstance(counts,float):
            counts = [counts]

        mag_inst = np.array([-2.5*np.log10(c)+zp if c > 0.1 else 0. for c in counts ])


    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno,e)

    return mag_inst

def gauss_fwhm(sigma):
   import numpy as np
   fwhm = 2*np.sqrt(2*np.log(2)) * sigma
   return fwhm

def moffat_fwhm(gamma, beta):
   import numpy as np
   x = np.sqrt((2**(1/beta))-1)
   y = 2*gamma*x
   return y

def pix_dist(x1,x2,y1,y2):
    import numpy as np
    z1 = (x1-x2)**2
    z2 = (y1-y2)**2
    z3 = np.sqrt(z1+z2)
    return z3

def gauss_2d(image, x0,y0, sky , A, sigma):
    import numpy as np
    (x,y) = image
    a = (x-x0)**2
    b = (y-y0)**2
    c = (2*sigma**2)
    d =  A*abs(np.exp( -(a+b)/c))
    e =  d + sky
    return  e.ravel()

def r_dist(x1,x2,y1,y2) :
    import numpy as np

    a = (x1-x2)**2
    b = (y1-y2)**2

    r = np.sqrt(a+b)

    return np.array(r)

def weighted_avg_and_std(values, weights):
    import numpy as np
    import math

    mask = ~np.ma.masked_array(values, np.isnan(values)).mask

    values = np.array(values)[mask]
    weights = np.array(weights)[mask]

    if len(values) == 0 or len(weights) == 0:
        values  = np.array(np.nan)
        weights = np.array(np.nan)


    average = np.average(values, weights=weights)

    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))


def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)


def weighted_median(data, weights):
    import numpy as np

    """
    Args:
      data (list or numpy.array): data
      weights (list or numpy.array): weights
    """
    data, weights = np.array(data).squeeze(), np.array(weights).squeeze()
    s_data, s_weights = map(np.array, zip(*sorted(zip(data, weights))))
    midpoint = 0.5 * sum(s_weights)
    if any(weights > midpoint):
        w_median = (data[weights == np.max(weights)])[0]
    else:
        cs_weights = np.cumsum(s_weights)
        idx = np.where(cs_weights <= midpoint)[0][-1]
        if cs_weights[idx] == midpoint:
            w_median = np.mean(s_data[idx:idx+2])
        else:
            w_median = s_data[idx+1]
    return w_median



