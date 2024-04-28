# ILMT-astrometry
This pipeline performs precise astrometry on the ILMT images, and updates the wcs information in the files.
Download a raw TDI frame from the ILMT cloud, and run the shell script 'run_astrometry.sh' in the same directory and pass the fits frame as the argument.

e.g., sh run_astrometry.sh *.fits

Following modules are required for running this pipeline:
novas,
astroquery,
sep,
astropy,
numpy,
lmfit etc.
