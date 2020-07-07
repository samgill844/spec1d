# SpUpNIC
A pipeline to reduce spupnic data, extract spectra and perform RV measurements. 


# Setup
There is no official setup for this - I might make one in future. For now, clone this repository and make an alias like this for simplicity. 

```bash
alias spupnic="python [path_to_spupnic]/SpUpNIC/spupnic.py"
```

# Usage

## Step 1 - create the log file

It is nesseary to re-make the log file from the fits files in the directory since sometime the log-files aren't made correctly at the 1.9-m at SAAO. First, cd into the directory containitng the fits files and then call spupnic.

```bash
cd ~/Documents/NGTS_monotransit_WG/SAAO/saao-1.9/0130
python ~/Software/SpUpNIC/spupnic.py --create_log
```

This will make a file called log.fits in the same directory which will be used for the reduction of data. From the log file, you can easily see which frames are science, calibration, dome flats or biases. 

## Step 2 - Bias subtraction and Flat-field corrections. 

Coming soon. For now, we don't use them and everything works pretty well. 

## Step 3 - Reduction of CuNe calibration spectra. 

We need to reduce the calibration spectra so that we can accurately calibrate the wavelength axis of the science frames. To do this, we need to pass a CuNe fits file name long with a path to the reference spectra. This code automatically extracts information about the grating and angle (e.g. gr5 and -4) from the fits headers. If for whatever this isnt right, you'll need to fix them.

```bash
python ~/Software/SpUpNIC/spupnic.py --extract_CuNe a1471096.fits --path_to_ref_spectra ~/Software/SpUpNIC/CuNe_ref_spectra
```
This will create 3 new files in our directory with the prefix a1471096 - two plots and a calibratation file. The first plot shows the CCD image along with the lines identified and a trace of the slit. The second plot shows the extracted spectrum and the associated lines marked up. In the output you should see something like this


![alt text](https://github.com/samgill844/SpUpNIC/blob/master/images/fig1.png)

![alt text](https://github.com/samgill844/SpUpNIC/blob/master/images/fig2.png)


```bash
Calibration:
	182.56421197218626 -> 621.728 nm
	268.83513631590824 -> 626.65 nm
	387.961762118894 -> 633.443 nm
	473.2948346762896 -> 638.299 nm
	507.18606197989385 -> 640.225 nm
	691.8217835455781 -> 650.653 nm
	738.8112541217184 -> 653.288 nm
	857.3613794904414 -> 659.895 nm
	1001.2105936018824 -> 667.828 nm
	1072.1899364353528 -> 671.704 nm
	1470.9127530198862 -> 692.947 nm
	1671.6209164671052 -> 703.241 nm
	1957.88402644295 -> 717.394 nm
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

Like the CuNe files, we can reduce these en-mass by passing:

```bash
spupnic --grating gr5 --gratingangle -4 --extract_all_science --threads 12
```