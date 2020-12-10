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

This will make a file called log.fits in the same directory which will be used for the reduction of data. From the log file, you can easily see which frames are science, calibration, dome flats or biases. 

## Step 2 - Bias subtraction and Flat-field corrections. 

Coming soon. For now, we don't use them and everything works pretty well. 

## Step 3 - Reduction of CuNe calibration spectra. 

We need to reduce the calibration spectra so that we can accurately calibrate the wavelength axis of the science frames. To do this, we need to pass a CuNe fits file name long with a path to the reference spectra. This code automatically extracts information about the grating and angle (e.g. gr5 and -4) from the fits headers. If for whatever this isnt right, you'll need to fix them by exliplicitly stating the grating and the grating angle (I would do this for sanity).

```bash
spupnic --extract_CuNe a2101095.fits --grating gr7 --gratingangle 17
```
This will create 3 new files in our directory with the prefix a1471096 - two plots and a calibratation file. The first plot shows the CCD image along with the lines identified and a trace of the slit. The second plot shows the extracted spectrum and the associated lines marked up. In the output you should see something like this


![alt text](https://github.com/samgill844/SpUpNIC/blob/master/images/fig1.png)

![alt text](https://github.com/samgill844/SpUpNIC/blob/master/images/fig2.png)


```bash
Calibration:
	1788.5604183193573 -> 811.5311 nm
	1751.2539996291457 -> 800.61567 nm
	1646.944371808619 -> 772.42072 nm
	1614.7128689877975 -> 763.5106 nm
	1568.5214421584471 -> 750.38691 nm
	1524.0965075784022 -> 738.39805 nm
	1410.6944539856065 -> 706.72181 nm
	1374.4206152596084 -> 696.54307 nm
	606.4163648145163 -> 476.48646 nm
	436.84098086506714 -> 427.75282 nm
	410.99737650954694 -> 419.8317 nm
```
These tell you which lines have been sucessfuly matched with the reference spectra. If there aren't any or very few, check the plots - something has gone wrong. This information is saved in a .npy extension which can be then used to calibrate nearby science frames. 

We can cheat here, and instead of giving it every single CuNe file seperately, we can search the log file and get it to reduce every CuNe frame in the current directory. To do this, we can do something like:

```bash
spupnic --grating gr5 --gratingangle -4 --extract_all_CuNe --threads 12
```

This will reduce all the CuNe frames using 12 multiprocessing threads.

## Step 4 - Reduction science frames. 

To reduce a science frame, we simply pass:
```bash
spupnic --grating gr5 --gratingangle -4 --extract_science a1341001.fits
```

This will look for calibration files within 1 hour of the science frame. If possiblem, it will take th eaverage of the calibration file before and after the exposure. If not, the default is simply to choose the nearest calibration file. If none are within 1 hour, it will fail unless the code is modified. 

![alt text](https://github.com/samgill844/SpUpNIC/blob/master/images/fig3.png)

![alt text](https://github.com/samgill844/SpUpNIC/blob/master/images/fig4.png)


Like the CuNe files, we can reduce these en-mass by passing:

```bash
spupnic --grating gr5 --gratingangle -4 --extract_all_science --threads 12
```

## Cosmic ray rejection

When processing long exposures, cosmic ray are a pain. We use astroscrappy to remove them from the image by interpolating from around the affected pixels. To use this feature, add the flag
```bash
----reject_cosmics
```
and it should work for you.