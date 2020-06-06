<p align="center">
  <img src=https://github.com/Astro-Sean/autophot/blob/master/logo.jpg>
</p>


<div align="center">

[![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/version.svg)](https://anaconda.org/astro-sean/autophot) [![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/latest_release_date.svg)](https://anaconda.org/astro-sean/autophot) [![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/license.svg)](https://anaconda.org/astro-sean/autophot) [![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/downloads.svg)](https://anaconda.org/astro-sean/autophot ) [![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/installer/conda.svg)](https://conda.anaconda.org/astro-sean)

</div>

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

## Road map

* Currently only able to search for images using PanSTARRs server - will extend to use SDSS query.
* Include [HEALPix](https://healpix.sourceforge.io/) to better optimize catalog selection, currently user defined.
* Add logging module inside of changing stdout in python.

## Version History
* 0.0.1
    Initial upload
