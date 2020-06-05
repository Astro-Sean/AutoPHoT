<p align="center">
  <img src=https://github.com/Astro-Sean/autophot/blob/master/logo.jpg>
</p>

[![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/version.svg)](https://anaconda.org/astro-sean/autophot) [![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/latest_release_date.svg)](https://anaconda.org/astro-sean/autophot) [![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/license.svg)](https://anaconda.org/astro-sean/autophot) [![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/downloads.svg)](https://anaconda.org/astro-sean/autophot )[![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/installer/conda.svg)](https://conda.anaconda.org/astro-sean)


## Description

A core part of my PhD is the development and publishing of AutoPHoT, an
automated photometry pipeline that is optimized for the measurements of
transients. As well as being used in my own group, AutoPHoT will be
accessible to the wider transient community through a web interface. As
part of the pipeline, image astrometry, metadata and quality will be verified
and calibrated onto a homogeneous standard system.

Massive stars and their resultant supernovae have played a crucial role in the
Universe. Many massive stars appear to undergo violent episodes of mass
loss, prior to gravitational collapse. The mechanisms that trigger these
mass loss episodes, known as 'supernova imposters’, are still poorly
understood. It is a focus of my group [sn.ie](http://sn.ie/) to understand these
mechanisms, as well as their connection to core-collapse supernova and
the physical mechanism that power them.

A first step in the investigating these objects is to measure the change
in brightness over time (photometry). A large amount of time is spent
by astronomers working through large quantities of datasets, performing
photometric measurements and calibrating these measurements to a
standard photometric system.

Project is currently in closed development although feedback is appreciated.
email: sean.brennan2@ucdconnect.ie [subject: github autophot]

## Developer notes

 Colour terms still in testing/development

## Installation

* Install via conda:

```python
conda install -c astro-sean autophot
```

* Code relies on [Astrometry.net](https://arxiv.org/abs/0910.2233) by Dustin Lang to solve for world coordinates system (WCS). Code can be downloaded/installed [here](http://astrometry.net/doc/readme.html) and [here](http://astrometry.net/doc/build.html#build.)
Once installed, locate the solve-field executable [default location: /usr/local/astrometry/bin/solve-field] and update (if needed) 'solve field exe loc' in input.yml.

* Image subtraction uses [HOTPANTS](http://www.ascl.net/1504.004) by Andy Becker - HOTPANTS can be found [here](https://github.com/acbecker/)
Known error with installation if installing on MacOS - if upon installation you get 'malloc.h' file not found, replace

```c
#include <malloc.h>
```
with
 ```c
 #if !defined(  MACH  )
 #include <malloc.h>
 #endif
 #if defined(  MACH  )
 #include <stdlib.h>
 #endif
```
to every .c file.
Once installed, locate the hotpants executable and update 'hotpants exe loc' in input.yml.
Future work will include pip install/upgrade installation method.
## Usage
* Main execution is from call.py
* Main settings are set in input.yml file
* To use file directory set out via 'fits dir' in input.yml:
```python
python call.py
```
* Will create separate folder title [directory name] output in parent directory.
* For single file:
```python
python call.py -f [filepath]
```
* Will save output in folder with name that of the fits file in same directory as fits file
* To use current directory:
```python
python call.py -d
```
* Will create separate folder title [directory name] output in parent directory.
*  to hide print statements
```python
python call.py -v
```
* Currently setup to reduce ACAM images from the WHT:
```python
python call.py -r
```
will reduce images from ACAM with bias and flats located in calib/
* If you want to add flats for ACAM images, keep same file notation i.e nflat_acam_[filter].fits
* Will also crop ACAM images to focus of circular illuminated region - if 'INSTRUMNE' is ACAM.
 Cannot currently handle input of lists using call.py -f [filepath] method, will be included in future development. If you want to run multiple files in a single instance, move to a folder and run autophot -d in that directory or change fits_dir in input.yml and run autophot.
## Output
* Cosmic ray cleaned image with updated header - [ original filename ] + APT:
	* Fwhm
	* Zeropoint
	* WCS values (if needed)
* Template image from Panstarrs (if possible) - [ original filename ] + template
	* Image from the [Panstarrs image cutout server](https://ps1images.stsci.edu/cgi-bin/ps1cutouts) in same filter as input image aligned astrometrically use wcs values.
* Subtracted image (if possible) - [ original filename ] +  subtraction
	Subtracted image used original image and template image using HOTPANTs.
* Astrometric data - [ original filename ] +  APT astrometry:
	* if wcs needed to be corrected/updated, The script will call astrometry.net and write a .txt file with astrometric outputs.
* Target photometric data - target output.csv
	* filepath
	* telescope
	* mjd
	* zeropoint [ zp [] ]
	* zeropoint error [ zp [] err ]
	* limiting magnitude [ lmag ]
	* instrumental magnitude of target [ [] inst ]
	* calibrated magnitude of target [ [] ]
	* calibrated magnitude of target err[ [] err ]
	* Signal to noise of target [SNR]
	* Photometric method used (ap/psf) [ method ]
* if user selected:
	* do all phot
		* Photometric data on sources in field (all sources and calibration sources) with sigma level set by do all phot sigma saved to .cat file.
	* save mag lim plot:
		* Comparison between catalog magnitude and recovered magnitude for sources in field.
	* save source plot:
		* Image of sources used for zeropoint (circle) with their catalog position (cross) convert to pixel values using wcs values.
	* save target plot:
		* Image of target:
			* if PSF used: shows source residuals after extraction.
			* if AP used: shows image and annuli.
	* save zp plot:
		* Graph of zeropoints found from catalog sources with sigma clipping.
  * get_template:
    * get photometric template [currently only from panstarrs server]
  * get_mag_lim:
    * Perform limiting magnitude script on local area around target via atifical star injected [requires PSF model to have been built]
  * phot_on_sub:
    * Perform photometry on subtracted image [requries template image]
## Road map
* Allow for user supplied template image
* Currently only able to search for images using PanSTARRs server - will extend to use SDSS query.
* Include [HEALPix](https://healpix.sourceforge.io/) to better optimize catalog selection, currently user defined.
* Add logging module inside of changing stdout in python.

## Version History
* 0.1
    * Initial Release 21-04-19:
* 0.2
	* Re-release 13-06-19:
		* New Features:
			* initial parsers integration.
			* Ability to download template images for pan starrs server.
			* Return photometry of all sources in field default @ 25 sigma.
* 0.3
	* Release 23-07-19:
		* New Features:
			* Introduced HOTPANTS.
			* ACAM reduction.
			* if number of fits files is >500 user will be asked if they want to continue.
			* added 2mass into catalog query.
		* Bugs Fixed:
			* Mitigated 'Tee object has no attribute issatty' error from reproject interp.
			* fixed remove wcs.
			* updated write yaml to not continue looking if filter keyword returns the word 'clear'.
* 0.4
  * Release 17-10-19:
    * Introduced limiting magnitudes
    * Introduced color terms
    * Adjusted image aligning and subtraction
