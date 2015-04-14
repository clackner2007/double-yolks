# Double Yolk Galaxies
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
-j set to the number of processors you want to use (python multiprocessing)
-r the inner radius of the ring in units of the FHWM of the point spread function
-w the width of the ring in pixels
-f the FWHM of the point spread function in pixels.
```
The input data must include: a galaxy ID number (unique), a filename for the image (not the path, unless different images have different paths), the redshift of the object, the apparent magnitude and the coordinates. The names of those columns must be as in input_test. The code will accept a fits table instead of an ASCII one, and you can include extra columns, but they will be ignored and not propogated. If the FWHM isn't included on the command-line, it needs to be included as a column called FWHM in the input file. This allows for different PSFs for each object.

`peak_filter.py` creates two files, the first is a list of galaxies call gal_list. This looks almost like the input list, but includes a column `N_PEAKS` which is the number of peaks found for that image. The second file is peak_list which has one line for each measured peak including the parameters (size, position, flux, etc.). At this point, there are a lot of junk peaks which have to be removed.

####clean_peaks.py
