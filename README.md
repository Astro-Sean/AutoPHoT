<p align="center">
  <img src=https://github.com/Astro-Sean/autophot/blob/master/logo.jpg>
</p>

<div align="center">

[![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/version.svg)](https://anaconda.org/astro-sean/autophot) [![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/latest_release_date.svg)](https://anaconda.org/astro-sean/autophot) [![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/license.svg)](https://anaconda.org/astro-sean/autophot) [![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/downloads.svg)](https://anaconda.org/astro-sean/autophot ) [![Anaconda-Server Badge](https://anaconda.org/astro-sean/autophot/badges/installer/conda.svg)](https://conda.anaconda.org/astro-sean)

</div>

## Introduction

The Automated Photometry of Transients (AutoPhOT) pipeline allows for rapid and automatic analysis of image data for transient events.

The pipeline is built from the ground up based on Python 3 and makes great us of astropy packages. AutoPhOT is able to handle homogenised data from different telescopes and applies techniques such as image calibration, image subtraction and novel PSF fitting in an automated and intelligent way.

**Project is currently in beta stage. Feedback is appreciated.
email: sean.brennan2@ucdconnect.ie [subject: github autophot]**

## Developer notes

 * Colour terms still in testing/development
 * Currently astrometry.net, HOTPANTS and astoalign needs to be installed by USER.

## Installation

* Install via conda:

```python
conda install -c astro-sean autophot
```

* Code relies on [Astrometry.net](https://arxiv.org/abs/0910.2233) by Dustin Lang to solve for world coordinates system (WCS). Code can be downloaded/installed [here](http://astrometry.net/doc/readme.html) and [here](http://astrometry.net/doc/build.html#build.)
Once installed, locate the solve-field executable [default location: /usr/local/astrometry/bin/solve-field] and update (if needed) 'solve field exe loc' in syntax (see [here](https://github.com/Astro-Sean/autophot/blob/master/autophot_example.ipynb)).

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
Once installed, locate the hotpants executable and update 'hotpants exe loc' in syntax see [here](https://github.com/Astro-Sean/autophot/blob/master/autophot_example.ipynb).

## Usage

For quick use see [here](https://github.com/Astro-Sean/autophot/blob/master/autophot_example.ipynb)

for more detailed explanation see here (work in progress)

## Road map

* update

## Flowcahrt

Basic operation of AutoPhOT for complete Photometric calibration of a given target

<p align="center">
  <img src=https://github.com/Astro-Sean/autophot/blob/master/flowchart.png>
</p>

## Version History

* 0.0.1
    Initial upload
