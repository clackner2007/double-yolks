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

`peak_filter.py` creates two files, the first is a list of galaxies call gal_list. This looks almost like the input list, but includes a column `N_PEAKS` which is the number of peaks found for that image. The second file is peak_list which has one line for each measured peak including the parameters (size, position, flux, etc.). At this point, there are a lot of junk peaks which have to be removed.

####clean_peaks.py
`peak_filter.py` is agnostic about what types of peaks it keeps and it lets through a lot of junk. `clean_peaks.py` filters out a lot of the junk. This piece of code has 2 configuration files. One, `configParams.py` is more global and used by the mock image tests as well. The other, `clean_config.ini` only applies to clean_peaks.py. These will have to be modified depending on your data.

To run: `clean_peaks.py input_list_file gal_list peak_list [options]`. The arguments/options are:
```
inputfile: input file table, as used in peak_filter.py
gal_list: galaxy listing from peak_filter ( called gal_list)
peak_list: peak listing from peak_filter (called peak_list)
param_file: config parameter file (*.ini) for cleaning peaks, clean_config.ini is provided
-p: give path to output
-x: output peak coordinates in pixels in clean_pairs.txt instead of ra/dec
-l: make images of peaks, and put them in output/path/imgs folder along with imgs.html file to easily open them
-i: path to input images, only needed for plotting
-e: make eps plots instead of png (default)
```

The out of `clean_peaks.py` are 4 text files

1. `all_peaks_sources.txt`: This contains all the sources and all their peaks before cleaning. Each line contains one source and the number of peaks as well as the separation of the peaks from the brightest sources and the fluxes of each peak. Note that the rows have different lengths.
2. `cleaned_peaks_sources.txt`: This contains only the cleaned peaks, including sources with just one peak. The format is the same as for `all_peaks_sources.txt`
3. `cleaned_peaks.txt`: This contains the positions of the peaks for each source. You can output it either in pixels or ra/dec. There's one source per line, so again, the lines have variable length.
4. `clean_pairs.txt`: This is a pair catalog only. It takes `cleaned_peaks_sources` and only returns galaxies with at least 2 clean peaks. For galaxies with more, it notes that in `N_PEAKS`, but only reports the brightest 2 peaks. 

#####configParams.py
This contains the cosmology parameters, as well as image parameters, like the magnitude zeropoint, the pixelscale (arcseconds/pixel) and the size of the cutout images in pixels. These are assumed to be global for a project.

#####clean_config.ini
This file contains configuration parameters relating to cleaning up the peak list. It's basically a list of cuts to be made on peaks. The values used here are arbitrary, but some of them can be set by looking at peak finding in mock galaxies and minimizing the contamination/maximizing the completeness. The parameters are as follows:

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
The second part of this code makes mock merger galaxies by coadding two real galaxy postage stamps with an offset. There are also tools to examine the completeness and contamination.

####make_merger_stamps.py

#####make_mocks.ini

####test_mock_recovery.py

#####mock_recovery.py


