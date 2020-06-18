<p align="center">
  <img src=https://github.com/Astro-Sean/autophot/blob/master/logo.png width="75%">
</p>

<div align="center">

[![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/version.svg)](https://anaconda.org/astro-sean/autophot) [![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/latest_release_date.svg)](https://anaconda.org/astro-sean/autophot) [![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/license.svg)](https://anaconda.org/astro-sean/autophot) [![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/downloads.svg)](https://anaconda.org/astro-sean/autophot ) [![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/installer/conda.svg)](https://conda.anaconda.org/astro-sean)

</div>

## Introduction

The Automated Photometry of Transients (AutoPhOT) pipeline allows for rapid and automatic analysis of image data for transient events.

The novel pipeline is built from the ground up, based on Python3 and makes use of Astropy packages. AutoPhOT is able to handle homogenised data from different telescopes and applies techniques such as image calibration, image subtraction and novel PSF fitting in an automated and intelligent way.

**Project is currently in beta stage. Feedback is appreciated.
email: sean.brennan2@ucdconnect.ie [subject: github autophot]**

## Developer notes

 * Colour terms still in testing/development and **not** in current build
 * Currently astrometry.net, HOTPANTS and astoalign needs to be installed by USER.
 * PSF selection will look for bright isolated sources, however this can lead to sources being selected near the image boundaries which can cause large errors, especially in IR bands. Future update with annulus selection procedure and spacial variant PSF.

## Installation

**Some users reporting installation hanging during solving environment - currently looking into this - as always creating a fresh conda environment solves the issue**

* Some packages require conda-forge in channel list:

```bash
conda config --add channels conda-forge.
```
as well as the astropy channel:

```bash
conda config --add channels astropy
```

* Install AutoPhOT via conda istall:

```bash
conda install -c astro-sean autophot
```

* Image alignment can use [astroalign](https://www.sciencedirect.com/science/article/pii/S221313372030038X) over WCS alignment from Astropy (using [reproject_interp](https://reproject.readthedocs.io/en/stable/api/reproject.reproject_interp.html)). Install via

```bash
pip install astroalign
```

* Code relies on [Astrometry.net](https://arxiv.org/abs/0910.2233) by Dustin Lang to solve for WCS. Code can be downloaded/installed [here](http://astrometry.net/doc/readme.html) and [here](http://astrometry.net/doc/build.html#build.).
Once installed, locate the solve-field executable [default location: /usr/local/astrometry/bin/solve-field] and update (if needed) 'solve_field_exe_loc' in syntax (see [here](https://github.com/Astro-Sean/autophot/blob/master/autophot_example.ipynb)). **If the user trusts their WCS this step can be ignore as Astrometry.net is not used.**


* Image subtraction a local instance of [HOTPANTS](http://www.ascl.net/1504.004) by Andy Becker - HOTPANTS can be found [here](https://github.com/acbecker/). Once installed, locate the hotpants executable and update 'hotpants_exe_loc' in syntax see [here](https://github.com/Astro-Sean/autophot/blob/master/autophot_example.ipynb). **If the user has no need for image subtraction this step can be ignored.**


**Known error with installation of HOTPANTS**

if installing on MacOS - if upon installation you get 'malloc.h' file not found, replace

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

## Usage

* Check out series of Jupyter Notebooks for AutoPhOT [here](https://github.com/Astro-Sean/autophot/tree/master/example_notebooks)

## Referencing

* As this code is very much still under development and being tested, you should be very cautious before using results from AutoPhOT in a publication - there are papers in prep that are using AutoPhOT light curves, but please speak to us first.

## Road map

* Awaiting user feedback


## Testing

* If you experience errors with a particular file, the most effective means of debug is to share the file with a developer for diagnostic. Once bugs have been addressed all files will be deleted. **All shared data will be kept confidential**.

## Flowcahrt

Basic operation of AutoPhOT for complete photometric calibration of a transient as shown [here](https://github.com/Astro-Sean/autophot/blob/master/autophot_example.ipynb).

<p align="center">
  <img src=https://github.com/Astro-Sean/autophot/blob/master/flowchart.png>
</p>

## Version History
* 0.1
    - Initial upload
* 0.1.1
    - Fixed error opening fits image
    - Adjusted directory creating sequence
    - added convergence parameters to world_2_pix in *main.py*
    - added version to *autophot/__init__.py*
