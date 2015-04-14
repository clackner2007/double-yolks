# Double Yolk Galaxies
Code to find double yolk galaxies and make mock images to run through the filter for testing. This code was used to find merging galaxies described in [this paper](http://arxiv.org/abs/1406.2327). The code breaks easily into two segments one to detect and filter the peaks and one to take images to make mock catalogs.

###Detecting the Peaks

####peak_filter.py
The first part of the code actually implements the ring filter which first smooths the image, then subtracts the smooth image (making a high-pass filter) and finally runs a detection algorithm to detect the separate peaks. Here's what that looks like for a well-chosen sample image. 

The options for peak_filter.py are all command-line options as follows:
* input list: this is the list of images you want to run. An example is in test_data/input_test. The input data must include: a galaxy ID number (unique), a filename for the image (not the path, unless different images have different paths), the redshift of the object, the apparent magnitude and the coordinates. The names of those columns must be as in input_test. The code will accept a fits table instead of an ASCII one, and you can include extra columns, but they will be ignored and not propogated.
* -p path for the output files, two are generated
* -i path to image input files. This will be combined with the filenames given in the input list
* -d
* -j
* -r
* -w
* -f
