# spec1d
A pipeline to reduce 1d dispersion spectrograph data. 


# Setup
Setup is as easy as,
```bash
pip install setup.py
```
There is no pypi for this yet. 

# Spupnic Usage

## Step 1 - create the log file

It is nesseary to re-make the log file from the fits files in the directory since sometime the log-files aren't made correctly at the 1.9-m at SAAO. First, cd into the directory containitng the fits files and then call spupnic.

```bash
cd ~/Documents/NGTS_monotransit_WG/SAAO/saao-1.9/0130
spupnic --create_log
```

This will make a file called log.fits in the same directory which will be used for the reduction of data. From the log file, you can easily see which frames are science, calibration, dome flats or biases. For now, do not worry if you are using a CuAr lamp, they all go through the CuNe flag for now at least. 

## Step 2 - list science frames 
```bash
(base) sam@Samuels-MBP test % spupnic --list_science
                     FILE                      HJD                   TARGET                 ARC-LAMP                  EXPTYPE
            a1341017.fits       2458850.3728416804             TIC259592689                      OFF                  SCIENCE
```

These are the spectra you can reduce. 

## Step 3 - Reduction of science images

To reduce each science image, we can use this,

```bash
spupnic --extract_science a1341017.fits --grating gr5 --gratingangle -4 [--reject_cosmics]
```
This will load the science image, trace it, subtract the sky spectra, search for the nearest 2 calibrations, find wavlength solutions for them by doing a distorted alignment of reference spectra and re-fitting peaks, average the wavelength solution if required and re-inteprolate the science spectra 

![alt text](https://github.com/samgill844/SpUpNIC/blob/master/images/fig1.png)


```bash
Averageing 2 calibrations
717.394 nm    208.46273186337336      207.8453878844751      [-0.6173439788982762]
703.241 nm    494.1781977925755      493.56004811897407      [-0.6181496736014083]
692.947 nm    694.5655700850274      693.9460883013087      [-0.6194817837186974]
671.704 nm    1092.8190418555323      1092.1784909380106      [-0.6405509175217503]
667.828 nm    1163.731655382685      1163.0918223639956      [-0.6398330186893872]
659.895 nm    1307.4731657003194      1306.824030664928      [-0.6491350353915095]
653.288 nm    1425.946404415987      1425.2941772162753      [-0.6522271997116604]
650.653 nm    1472.9110233082913      1472.2483716000277      [-0.6626517082636383]
640.225 nm    1657.459241984377      1656.7986628509263      [-0.6605791334507103]
638.299 nm    1691.3626551538723      1690.6795491346786      [-0.6831060191937013]
633.443 nm    1776.7321094173828      1776.0630613542303      [-0.6690480631525588]
626.65 nm    1895.9620870699257      1895.2871650726054      [-0.6749219973203253]
```

This will save a *_spectra.fits file in your current directory which the first extension is wavelength, flux, flux_err. Other extensions include 

ext0 - wavelength, flux, flux_err
ext1 - pixel x axis
ext2 - raw sky flux
ext3 - lower sky flux 
ext4 - uper sky flux 
ext5 - average sky flux which is subtraced from the stellar spectra.
ext6 - the traced calibration spectra for the first reference lamp
ext7 - the traced calibration spectra for the second reference lamp.


## Cosmic ray rejection

When processing long exposures, cosmic ray are a pain. We use astroscrappy to remove them from the image by interpolating from around the affected pixels. To use this feature, add the flag
```bash
----reject_cosmics
```
and it should work for you.

# Process entire nights

To make things easy, we can process entire nights easily
```bash
spupnic --process_folder --grating gr5 --gratingangle -4 --reject_cosmics --threads 12
```