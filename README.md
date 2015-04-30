## Double Yolk Galaxies
Code to find double yolk galaxies and make mock images to run through the filter for testing. This code was used to find merging galaxies described in [this paper](http://arxiv.org/abs/1406.2327). The code breaks easily into two segments one to detect and filter the peaks and one to take images to make mock catalogs.

###Detecting the Peaks

####peak_filter.py
The first part of the code actually implements the ring filter which first smooths the image, then subtracts the smooth image (making a high-pass filter) and finally runs a detection algorithm to detect the separate peaks. Here's what that looks like for a well-chosen sample image. 

![alt text](https://raw.githubusercontent.com/clackner2007/double-yolks/master/ringfilt_img_readme.png "Example of Peak Filter")

The options for peak_filter.py are all command-line options as follows:
```
input list	this is the list of images you want to run. An example is in test_data/input_test. 
-p 	path for the output files, two ascii text lists are generated
-i 	path to image input files. This will be combined with the filenames given in the input list
-d set to the redshift offset that's been applied to the images if you are running simulation tests
-j set to the number of processors you want to use on your computer 
-r the inner radius of the ring in units of the FHWM of the point spread function
-w the width of the ring in pixels
-f the FWHM of the point spread function in pixels.
```
To run, it's just `python peak_filter.py input_list_file [options]`. Don't rely on the default options as they aren't likely to make sense for your data.

The input data must include: a galaxy ID number (unique), a filename for the image (not the path, unless different images have different paths), the redshift of the object, the apparent magnitude and the coordinates. The names of those columns must be as in input_test. The code will accept a fits table instead of an ASCII one, and you can include extra columns, but they will be ignored and not propogated. If the FWHM isn't included on the command-line, it needs to be included as a column called FWHM in the input file. This allows for different PSFs for each object.

`peak_filter.py` creates two files, the first is a list of galaxies call gal_list. This looks almost like the input list, but includes a column `N_PEAKS` which is the number of peaks found for that image. The second file is peak_list which has one line for each measured peak. The entries for each peak are:

1. ID of the galaxy
2. DELTA_Z: redshift offset, almost always zero
3. PEAK_ID: number of the peak in the galaxy
4. PEAK_NPIX: number of pixels above the noise threshold detected as part of each peak
5. PEAK_B-A: the ellipticity of the peak (in the form of an axis ratio). Later, we will eliminate extremely elongated peaks as noise
6. PEAK_X0_PIX/PEAK_Y0_PIX: the luminosity-weighted central position of the peak
7. PEAK_FLUX_FILTERED/PEAK_FLUX: the flux of the peak on the filtered and unfiltered image. This simply counts up the flux in the pixels and isn't a great measurement. The FILTERED_FLUX version is only used internally.

At this point, there are a lot of junk peaks which have to be removed.

####clean_peaks.py
`peak_filter.py` is agnostic about what types of peaks it keeps and it lets through a lot of junk. `clean_peaks.py` filters out a lot of the junk. Because the cleaning is separate from the detecting, it can be run many times with different configuration parameters on the same set of input data from `peak_filter.py`. This code has 2 configuration files. One, `configParams.py` is more global and used by the mock image tests as well. The other, `clean_config.ini` only applies to clean_peaks.py. These will have to be modified depending on your data.

To run: `clean_peaks.py input_list_file gal_list peak_list param_file [options]`. The arguments/options are:
```
inputfile: input file table, as used in peak_filter.py
gal_list: galaxy listing from peak_filter ( called gal_list)
peak_list: peak listing from peak_filter (called peak_list)
param_file: config parameter file (*.ini) for cleaning peaks, clean_config.ini is provided
-p: give path to output
-x: output peak coordinates in pixels in clean_pairs.txt instead of ra/dec
-l: make images of peaks, and put them in output/path/imgs folder along with imgs.html file to easily open them (this step is SLOW)
-i: path to input images, only needed for plotting
-e: make eps plots instead of png (default)
```

The outputs of `clean_peaks.py` are 4 text files

1. `all_peaks_sources.txt`: This contains all the sources and all their peaks before cleaning. Each line contains one source and the number of peaks as well as the separation of the peaks from the brightest sources and the fluxes of each peak. Note that the rows have different lengths.
2. `cleaned_peaks_sources.txt`: This contains only the cleaned peaks, including sources with just one peak. The format is the same as for `all_peaks_sources.txt`
3. `cleaned_peaks.txt`: This contains the positions of the peaks for each source. You can output it either in pixels or ra/dec. There's one source per line, so again, the lines have variable length.
4. `clean_pairs.txt`: This is a pair catalog only. It takes `cleaned_peaks_sources` and only returns galaxies with at least 2 clean peaks. For galaxies with more, it notes that in `N_PEAKS`, but only reports the brightest 2 peaks. 

#####configParams.py
This contains the cosmology parameters, as well as image parameters, like the magnitude zeropoint, the pixelscale (arcseconds/pixel) and the size of the cutout images in pixels. These are assumed to be global for a project.

#####clean_config.ini
This file is an example configuration file for `clean_peaks.py`. It's basically a list of cuts to be made on peaks. The values used here are arbitrary, but some of them can be set by looking at peak finding in mock galaxies and minimizing the contamination/maximizing the completeness. The parameters are as follows:

1. `ecut` this sets the minimum ratio of the ellipticities (as computed from the second moments of the flux distribution). It selects *against* long skinny peaks, which are usually artificial.
2. `dist_cut_1` inner separation limit in kpc. This sets the minimum allowed separation. If you want a sample that's similarly complete at all redshifts, set this to the kpc-scale of the resolution limit at the highest redshift.
3. `dist_cut_2` outer separation limit in kpc. This is driven by the image size (in kpc) at the lowest redshift. Making it larger will introduce more spurios line-of-sight pairs.
4. `flux_cut` this is the minimum allowed flux ratio of the peaks compared to the brightest peak. It's supposed to select major mergers, but doesn't work really well for that in practice. It is needed to eliminate non-galaxy (spiral structure, etc.) peaks
5. `totflux_cut` this is the minimum allowed flux of the peak compared to the total source flux. Again, it's there to remove contamination and the exact value is driven by simulations.
6. `cent_dist_cut` this is the maxiumum allowed offset from the center of the image for any of the peaks in kpc. Basically, we want to be looking at the object in the center of the frame, not ones that happen to be off to the side.
7. `imsize_x` 1st dimension of image in arcseconds (could be removed...I think)
8. `imsize_y` 2nd dimension of image in arcseconds
9. `pearsonr_cut` this removes sources with 3+ peaks in which the peaks are spatially-aligned as we don't expect that to be the case in reality, but it does occur for barred spirals or edge-on disks.


###Testing with Mock Mergers
The second part of this code makes mock merger galaxies by coadding two real galaxy postage stamps with a spatial offset. This mimicks real merging systems, but it doesn't create any of the known structure that exists in mergers (tidal tails, starbursts, shells, rings etc.). Many of these features are low surface brightness, making them unimportant for the peak filter. Once the mock images are generated, they can be run through the peaking finding and cleaning to test the completeness of the peak finder. It's useful to play around with the cleaning parameters at this stage to see how they effect the completeness and contamination.

####make_merger_stamps.py
This is the code that makes the mock merger images. It can use the same images you ran the `peak_filter.py` on, or use entirely different images. Each merger image consists of 2 randomly selected galaxies. The galaxies are required to have similar redshifts (within 0.03 by default) and their morphologies can be constrained. One galaxy image is then offset by a fixed amount in a random direction. The appropriate amount of noise is added to the blank pixels and then the two galaxy images are added together. This means the mock images are sqrt(2) times noisier than the input images.

To run: `make_merger_stamps.py input_list_file param_file [options]`. The arguments/options are:
```
inputfile: input file table, same format as in peak_filter.py
param_file: config parameter file (*.ini) for making mocks, make_mocks.ini is provided clean_config.ini is provided
-o: output path for listing files and images
-i: path to input images, these are the images that will be coadded
```

The outputs of `make_merger_stamps.py` are the coadded images (FITS format) in output/path/imgs. The files names are of the form `GALID1_GALID2_sepkpc.fits`, which gives the ids of the two galaxies in the image and their separation in kiloparsecs (assuming redshift of the middle galaxy, but the redshifts are similar). The code also outputs two text files into the directory specified by the `-o` option:

1. `input_peakfilter_MAG.txt` This is the input file needed to run the peak-finder on the mock images. It contains the new id numbers (starting with 0), the image filenames, the redshift (of the 'middle' galaxy) and the magnitude obtained by summing the fluxes of the two input galaxies. `MAG` in the filename corresponds to the limiting magnitude used to generate the mocks and given in the configuration file.
2. `simulatedSample_MAG.dat` This file gives more details about the mock merger images that can be compared to the output of the peak-finding code above. It includes the id numbers of the mock merger, the input galaxies, the magnitudes of the galaxies and the magnitude of the merger (the sum of the galaxy fluxes), the flux ratio of the galaxies, the redshifts of the galaxies, the morphology (zest) parameters of the galaxies, the separation of the galaxies in the mock image, and the position of the galaxies in pixels in the mock image. The first galaxy is always in the center of the image while the second is offset in a random direction. The morphology parameters can be given in the input list for the mock generation. If they aren't the code just puts a 1 here and doesn't use it for anything.

#####make_mocks.ini
This is an example of the configuration parameter file used by `make_merger_stamps.py`. The available parameters are:

1. `list_stamps`
  1. `ngal` is the number of galaxy pairs to coadd. This times the number of offsets will give the number of mock merger images.
  2. `mag_limit` is the limiting magnitude to use for the input sample. It also shows up in the output table filenames. This maybe useful if you don't want to coadd galaxies which are extremely faint.
2. `zlims` give the redshift limits applied to the input galaxies. Only galaxies between 0.2<z<1.1 will be used.
3. `offsets` lists the offsets to be used in kpc. The number of parameters here is arbitrary, but the names should be of the form `oNUMBER`. The code uses the input galaxy redshifts to convert these to offsets in arcseconds for each galaxy pair. Mock images are made using each offset for each of the `ngal` pairs of galaxies.
4. `morph_class` lists the possible morphology classes. The code expects integer labels here. Originally, this morphology label referred to the ZEST morphology class of the galaxy, but the label could mean anything. In order to use this, the input galaxy list needs to have a column called ZEST with the integer label for each galaxy. By placing values in this list, the code will only find matching galaxies which have a ZEST parameter given in the list. For example if the `morph_class` values are `[1,2]`, the code will only find matches in which at least one galaxy is of ZEST type 1 or 2. If there is no ZEST column in the data, the code uses a default value of 1, so make sure that the configuration always includes at list the value 1 if you don't plan on using it.


####test_mock_recovery.py
After running the mocks through the peak filtering and cleaning code, it's time to look at the outputs and see what the recovery fraction is and how the measured separations and flux ratios compare with the real thing. The following code is really just a basic set of tools to get started with this, far from the final word. `test_mock_recovery.py` basically shows examples of tests and plots one might want to do and `mock_recovery.py` contains some the machinery needed to do these tests. 

The call signature for `test_mock_recovery.py` is:
```
inputfolder: the folder containing the cleaned*.txt sources from cleaned_peaks.py run on the mocks
origfile: file listing properties of the mock mergers. This is generated by make_merger_stamps and is simulatedSample*.dat
-o: output path where the plots and images will be written
-e: make eps plots instead of png plots
-m NIMAGES: make plots of the mocks with the detections shown (like in clean_peaks) for the first NIMAGES plots, default is to do nothing
```
One important concept in this code is that of the `restrict` sample. This is a subsample of the total set of mock mergers. In particular, this sample only contains mock mergers for which the input separations and flux ratios are within certain ranges. For example, we know that the peak finder can't detect mergers at very small separations, but the mock catalog may contain such mergers. In order to look at the completeness in systems for which we think the peak finding *will* work, we need to start with a restricted input sample. One could imagine putting all kinds of restrctions (only blue galaxies, only early-type, etc.) here. `test_mock_recovery.py` outputs a figure called completeness_scatter with 6 subplots:

1. The completeness as a function of redshift. This shows all the detected mock mergers as a function of redshift. The thick cyan line is the completeness in the restricted sample, while the thin one is in the total sample. The restricted sample completeness is usually better if you are restricting to sample where the pairs are easier to measure. The completeness overplots the total redshift distribution of the sample
2. The completeness as a function of separation. In this plot, any restrctions on separation are ignored, but other restrctions (flux ratio, by default) are still in effect. The lower portion of the panel shows the fractional difference in the measured and input separations. The separations are in kpc.
3. The completeness as a function of input flux ratio. The lower panel shows the fractional difference in the measured and input flux ratios.
4. The measured versus input flux ratio.
5. As with panel (1), but the completeness as a function of magnitude of the merging system.
6. As with panel (4), but for the input and measured separations.
These figures give some idea how complete the merger sample is and how well different parameters are measured. All of these plots are based on the cleaned outputs. This allows you to test the cleaning by using different cleaning parameters for the mocks and remaking these plots.

In addition to the completeness, the contamination is also important. The contamination can be divided into two main types: contamination from non-merging galaxies (spiral arms, star formation etc.) and contamination from merging galaxies with separations or flux ratios outside the 'desired range'. The latter type is important because is depends on redshift. At low redshift, it's much easier to detect mergers at small separations, so the contamination from mergers at small separations increases at low-z. `test_mock_recovery.py` does a few simple things to look at the contamination, but there are many more possible tests here. First it prints a set of numbers:
1. the total number of mergers (in the main and restricted samples)
2. the number of pairs measured by the peak filtering and cleaning (length of `clean_pairs.txt` file)
3. the number of measured pairs that correspond to real galaxies. The correspondence is determined by the position of the peak; if the peaks lines up with the galaxies in the mock merger image, then it's a real pair
4. the number of measured pairs that don't correspond to the mock merger.
5. out of the restricted sample, the number of real pairs (3 on this list) that have separations and/or flux ratios outside the limits set by the `restrict` parameter. For instance, if our sample is restricted to mergers separated by at least 2 kpc, detected mergers separated by 1 kpc are a contaminant.

There is one contamination plot the code makes: `peak_tot_flux_ratio.png`. This is a histogram of the peak-to-total flux ratios for peaks corresponding to the galaxies and peaks which are something else (unmasked sources, star-formation, etc.). Setting the total_flux_cut in the cleaning will control the contamination from extra peaks and this plot just shows where best to put that cut.

#####mock_recovery.py
This code has a lot of the work behind the plots in `test_mock_recovery.py` and should be imported when writing new tests of mock galaxies. The classes and functions it has are:

* `readMergers`: function to read the merger data from the `clean*.txt` outputs and the `simulatedSample*.dat` input file. This function takes the output data and reorganizes it into a list of `Merger` objects. These merger objects contain information about the peaks detected in each merger, as well as the galaxies that went into the mock merger. There is an option to clean the peaks, but because the peaks are read from the clean outputs, it's not exercised.
* `plotComplete`: function to plot the completeness as a function of a given property. The property can either by intrinsic (redshift) or measured (flux ratio). Examples of using this function are in `test_mock_recovery.py`. 
* `Merger`: class to contain a single mock merger and it's properties. This class holds the basic properties of the mock merger, i.e. magnitude, redshift, input separation of galaxies, input flux ratio, etc. It also contains the list of peaks detected and the properties (size, position, flux, ellipticity) of the peaks. Finally, and most importantly, it identifies which, if any of the peaks, are associated with the input galaxies that made the merger. The association is based on position; the real galaxy centroid and the peak center have to be within 8 pixels (see `Mergers.assignRealPeaks`). The peaks are then divided into `MeasP1`, `MeasP2` and `extrapeaks`. `Merger.isdbl` is set to true if both the input galaxies are detected as peaks (`MeasP1` and `MeasP2`). These attributes and their related methods are the part of the `Merger` class most useful to computing the completeness and the contamination.
