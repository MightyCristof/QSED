# qsed
# 
# Fast SED modeling with low-resolution galaxy+AGN templates.
# 
# Christopher M. Carroll
# 
#
# Steps to run SED modeling
#
# 1. Add directory to !PATH using the your idl_startup file or similar.
# 
#    If you'd like to run without permanently adding to your !PATH (i.e., case-by-case basis),
#    move to the "qsed/" directory and run the following from the command line:
# 
#         .r scripts/add_dir2path
# 
# 
# 2. Bandpass information is provided for all instruments in the "bandpass" directory, and 
#    is saved under "bandpass.sav". 
#    
#    Adding another instrument is a two step process:
#    First, add a subdirectory within "bandpass/" and add the bandpass files to the new 
#    folder. Make sure the bandpass files are labeled in such a way that they are in 
#    ascending order of increasing wavelength.
#         Ex: WISE  -- RSR-W1.txt, RSR-W2.txt, ... (3.4-micron, 4.6-micron, ...)
#         GALEX -- 01-GALEX.FUV.dat, 02-GALEX.NUV.dat (135Ð175 nm, 175Ð280 nm)
# 
#    Second, modify the script "MAGFLUX.PRO" to include the new instrument with associated
#    zero points for magnitude-flux conversion.
# 
#    Then run the following from the command line:
# 
#         merge_bp
# 
# 
# 3. The galaxy+AGN templates provided are ready to go. The current setup uses the 
#    Assef+2010 elliptical, irregular, and AGN template, combined with the Kirkpatrick+2015 
#    star-forming template.
#    
#    To add a different templates, modify the script "PREP_RAW_TEMPLATES.PRO", and run 
#    the following from the command line:
#    
#         prep_raw_templates
#    
# 
# 4. This SED modeling using a grid-based approach to find the best-fit model per source. 
#    The SED grid models are stored in the "models/" directory. A sample grid ("galtemp.sav") 
#    has been provided that uses GALEX/SDSS/UKIDSS/2MASS/WISE, and covers a parameter space
#    of 0 < E(B-V)_AGN < 50 with spacing ÆE(B-V) = 0.05 dex, and 0.00 < z < 3.00 with 
#    spacing Æz = 0.01.
#    
#    To create a different "galtemp", use the procedure "create_galtemp.pro", 
#    and input your preferred color excess, redshift, and instruments. Run the following 
#    from the command line:
#    
#         load_vars,'templates/comp.sav','_comp'
#         load_vars,'bandpass/bandpass.sav','_bp'
#         create_galtemp,<OUTPUT_FILE_NAME>,<EBV_ARRAY>,<Z_ARR>,<ARRAY_OF_INSTRUMENTS>
# 
#    EX: To create the provided "galtemp.sav", I run the command as:
#         ebv = 10.^(dindgen(55)/20-1)-0.1d
#         red = dindgen(3000)/1000.
#         inst = ['galex','sdss','uk','twom','wise']
#         create_galtemp,'galtemp.sav',ebv,red,inst
# 
# 
# 5. Prep photometry for SED modeling. The procedure "READ_CARROL2021.PRO" is specific to 
#    my own dataset for Carroll+2021. I've provided it for reference, as well as a more 
#    generalized procedure, "READ_PHOTOMETRY.PRO". The procedure assumes the data is 
#    contained in a .FITS file and in the "data/" directory. It also assumes the names of 
#    data columns passed to it (i.e., "RA", "DEC", etc.). You will need to change those if 
#    your data is not formatted the same.
#    
#    A test data set has been provided with 50 sources, 49 of which should pass the criteria
#    in "READ_PHOTOMETRY.PRO". Run the following from the command line:
#    
#         read_photometry,'test_data.fits'
#         
#    This should produce the file "test_data_READY.sav" which contains the photometry prepped
#    and ready for SED modeling.
#    
# 
# 6. Run the SED modeling! The SED modeling resamples the photometric errors on each source
#    a number of times specified on input. The procedure is pretty fast, so you can test the
#    speed with your own dataset to decide the number of iterations you need.
#    
#    Run the following from the command line:
# 
#         load_vars,'models/galtemp.sav','_galtemp'
#         load_vars,'templates/comp.sav','_comp'
#         qsed_resamp,'test_data_READY.sav','galtemp.sav','comp.sav',1000,/verbose
# 
#    For very large datasets, the /VERBOSE keyword which will print an alert to screen 
#    every time the procedure completes 10% of the sources.











