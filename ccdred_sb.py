"""A module for performing CCD data reduction steps. 

`ccdred_sb` is modeled after IRAF's `imred.ccdred` package. It contains a set
of functions for creating dark frames, creating a flats, and processing object
image files.

Author: Shane Burns
Date: 7/26/16
Modified: 10/6/21, 9:21 AM
"""
import numpy as np
from astropy.io import fits
from sys import exit
import datetime as dt

def get_image(file_path, list_file_name):
    """Get image data from a collection of FITS files listed in 
    `list_file_name`. Gets *only* the image in the primary `hdu`.
    
    Parameters
    ----------
    file_path, list_file_name : str
    
    Return
    ------
    image_data_list : list of numpy arrays
        List where each element is a numpy array containing the image
        data from one of the files in `list_file_name`.
    """
    list_path = list_file_name
    
    file_names = open(list_path,'r')
    
    # Get the image data and header from each file in the file list.
    image_data_list = []
    image_head_list = []
    for fname in file_names:
        try:
            fname = fname.strip()
            fpath = file_path+fname
            print('Reading %s' % fname)
            hdulist = fits.open(fpath)
            image_data_list.append(hdulist[0].data)
            image_head_list.append(hdulist[0].header)
            hdulist.close()
        except:
            print("ERROR: Can't open file %s" % fname)
    file_names.close()
    
    return image_data_list,image_head_list

def zerocombine(imdir='images/',input_list='zero_list.txt',output='Zero'):
    """Like IRAF's `zerocombine`. It currently combines the frames by
    averaging all of the bias frames listed in the `input_list`. It writes
    out a FITS file with the name given by the `output` parameter. Prints the
    `mean` and `sigma` (standard deviation) of the output file's pixel values.
    
    Parameters
    ----------
    imdir : str
    	Path to image directory (must include trailing '/')
    input_list : str
        Name of file containing list of filenames
    output : str
        Name for the output FITS file.

    Returns
    -------
    fout_name : str
        Name with `.fits` extension of created file
    """
    # Get the current path and open the file-list file 
    image_data_list,image_head_list = get_image(imdir,input_list)
    
    # Find the mean bias image 
    zero_data = np.mean(image_data_list,axis=0)
    
    # Write the master bias FITS file    
    fout_name = output+'.fits'
    fout_path = imdir+fout_name
    
    hdu = fits.PrimaryHDU(zero_data)
    # This assumes that the `gain` and `rdnoise` of all images is the same
    try:
        hdu.header['gain'] = image_head_list[0]['gain']
    except:
        print("Warning: 'GAIN' header keyword not found.")

    try:
        hdu.header['rdnoise'] = image_head_list[0]['rdnoise']
    except:
        print("Warning: 'RDNOISE' header keyword not found.")
        
    hdu.header['ccdtype'] = 'zero'        
    utc = dt.datetime.utcnow()
    time_stamp = utc.strftime('%b %d, %Y at %H:%M:%S UTC')
    hdu.header['comment'] = 'Master bias frame created ' + time_stamp
    hdulist = fits.HDUList([hdu])
    try:
        hdulist.writeto(fout_path)
        print('zerocombine output file name: %s' % fout_name)

        # Compute return values    
        pix_mean = np.mean(zero_data)
        pix_sigma = np.std(zero_data)
        print("    mean pixel value = %f" % pix_mean)
        print("    standard deviation = %f" % pix_sigma)  
        return fout_name
    except:
        print('\n*** Failed to write %s. ***' % fout_name)
        print('*** The file may already exist. ***')
        return None
    
#Begin**********************
def darkcombine(imdir='images/',input_list='dark_list.txt',output='Dark'):
    """Like IRAF's `darkcombine`. Creates dark current images and then it
    currently combines the frames by averaging all of the dark current
    frames listed in the `input_list`. It writes out a FITS file with the name
    given by the `output` parameter. 
    
    Parameters
    ----------
    imdir : str
    	Path to image directory (must include trailing '/')
    input_list : str
        Name of file containing list of filenames
    output : str
        Name for the output FITS file.

    Returns
    -------
    fout_name : str
        Name with `.fits` extension of created file
    """
    # Get the current path and open the file-list file 
    image_data_list,image_head_list = get_image(imdir,input_list)
    
    # Divide each data file by the exposure time for that file to 
    # create dark current frames.
    Ndarks = len(image_data_list)
    dark_frames=[]
    for idx in range(Ndarks):
        exposure_time = image_head_list[idx]['exptime']
        print('the exposure time is %g' % exposure_time)
        dark_current = image_data_list[idx]/exposure_time
        dark_frames.append(dark_current)
    dark_data = np.mean(dark_frames,axis=0)
    print('dark_data type %s'%type(dark_data))
            
    # Write the master bias FITS file    
    fout_name = output+'.fits'
    fout_path = imdir+fout_name
    
    hdu = fits.PrimaryHDU(dark_data)
    # This assumes that the `gain` and `rdnoise` of all images is the same
    try:
        hdu.header['gain'] = image_head_list[0]['gain']
    except:
        print("Warning: 'GAIN' header keyword not found.")

    try:
        hdu.header['rdnoise'] = image_head_list[0]['rdnoise']
    except:
        print("Warning: 'RDNOISE' header keyword not found.")
        
    hdu.header['ccdtype'] = 'dark'
    hdu.header['exptime'] = 1        
    utc = dt.datetime.utcnow()
    time_stamp = utc.strftime('%b %d, %Y at %H:%M:%S UTC')
    hdu.header['comment'] = 'Master dark current frame created ' + time_stamp
    hdulist = fits.HDUList([hdu])
    try:
        hdulist.writeto(fout_path)
        print('darkcombine output file name: %s' % fout_name)

        # Compute return values    
        pix_mean = np.mean(dark_data)
        pix_sigma = np.std(dark_data)
        print("    mean pixel value = %f" % pix_mean)
        print("    standard deviation = %f" % pix_sigma)  
        return fout_name
    except:
        print('\n*** Failed to write %s. ***' % fout_name)
        print('*** The file may already exist. ***')
        return None
#end************************

def flatcombine(imdir='images/',input_list='flat_list.txt',output='Flat',zero_sub = True,
                zero_file = 'zero.fits'):
    """Like IRAF's `flatcombine`. 
    
    It currently works by: 
    
    1. Subtract `zero_file` from each file in `input_list` if `zero_sub = True`
    2. Compute median pixel value for each file in `input_list`
    3. Create the master efficiency map file with name specified by `output`
    
    Parameters
    ----------
    imdir : str
    	Path to image directory (must include trailing '/')
    input_list : str
        Name of file containing list of filenames
    output : str
        Name for the output FITS file.
    zero_sub : boolean
        Flag for bias subtraction
    zero_file : str
        Name of bias file for bias subraction        
        
    Returns
    -------
    fout_name : str
        Name with `.fits` extension of created file
    """

    # Get `zero_file` data if needed
    if zero_sub:    
        zero_path = imdir + zero_file
        hdu_list = fits.open(zero_path)
        zero_img = hdu_list[0].data
        hdu_list.close()

    # Get the current path and open the file-list file 
    image_data_list,image_head_list = get_image(imdir,input_list)

    # Create list of arrays for processed flats
    flat_data_list = []        

    # Process individual flats and load into `flat_data_list`
    for image_data in image_data_list:
        if zero_sub:        
            temp_data = image_data - zero_img
        else:
            temp_data = image_data
        med = np.median(temp_data)
        print('Flat image meadian = %f' % med)
        temp_data = temp_data/med
        
        flat_data_list.append(temp_data)

    # Create master flat    
    flat_data = np.median(flat_data_list,axis=0)
        
    hdu = fits.PrimaryHDU(flat_data)
    # This assumes that the `filter`, `gain` and `rdnoise` of all images 
    # are the same.
    try:
        ast_filter = image_head_list[0]['filter']       
        hdu.header['filter'] = ast_filter
        fout_name = output+'_'+ast_filter+'.fits'
    except:
        hdu.header['filter'] = 'none'
        fout_name = output+'.fits'
        
    try:
        hdu.header['gain'] = image_head_list[0]['gain']
    except:
        print("Warning: 'GAIN' header keyword not found.")

    try:
        hdu.header['rdnoise'] = image_head_list[0]['rdnoise']
    except:
        print("Warning: 'RDNOISE' header keyword not found.")
        
    hdu.header['ccdtype'] = 'flat'
    utc = dt.datetime.utcnow()
    time_stamp = utc.strftime('%b %d, %Y at %H:%M:%S UTC')
    hdu.header['comment'] = 'Master flat frame created ' + time_stamp
    hdulist = fits.HDUList([hdu])
    
    # Write the output FITS file    
    fout_path = imdir+fout_name
    
    try:
        hdulist.writeto(fout_path)   
        print('flatcombine output file name: %s' % fout_name)
        return fout_name
    except:
        print('\n*** Failed to write %s. ***' % fout_name)
        print('*** The file may already exist. ***')
        return None

def ccdproc(imdir='images/',input_list='object_list.txt', output_list = None,
            zero_file='Zero.fits', flat_file='Flat.fits',
            dark_correction=False,dark_file='Dark.fits'):
    """Like IRAF's `ccdproc`. 
    
    It currently works by: 
    
    1. Subtract `zero_file` from each file in the `input_list` 
    2. Divide by `flat_file` from each file in the `input_list`
    3. Write processed FITS files with name specified by `output_list`
    
    Parameters
    ----------
    input_list : str
        Name of file containing list of input file names
    output_list : str
        Name of file containing output file names. If `output_list = None` the
        filenames are the same as those in `input_list` with `_p` appended.
    zero_file : str
        Name of bias correction file
    flat_file : str
        Name of flat correction file       
    dark_correction : boolean
        applies dark correction if True
    dark_file : str
        dark current image file name
        
    Returns
    -------
    fout_list : list of str
        List of names of the processed images.
    """
    print('Running ccdproc')
    
    # Get get object images
    image_data_list,image_header_list = get_image(imdir,input_list)

    # Get calibration images
    zero_path = imdir + zero_file
    hdu_list = fits.open(zero_path)
    zero_img = hdu_list[0].data
    hdu_list.close()

    flat_path = imdir + flat_file
    hdu_list = fits.open(flat_path)
    flat_img = hdu_list[0].data
    hdu_list.close()

    # Create output filenames
    out_name_list = []
    if (output_list == None):
        list_path = input_list
        file_names = open(list_path,'r')
        for name in file_names:
            parts = name.split(sep='.')
            out_name = parts[0]+'_p.fits'
            out_name_list.append(out_name)
    else:
        list_path = imdir + output_list
        file_names = open(list_path,'r')
        for name in file_names:
            out_name = name.strip()
            out_name_list.append(out_name)
    file_names.close()
    
    N_files = len(image_data_list)
    if (len(out_name_list) != N_files):
        print('ERROR!')
        print('Number of file names does not equal the number of images')
        exit
        
    for idx in range(N_files):
        # Process images
        exposure_time = image_header_list[idx]['exptime']
        if dark_correction:
            dark_path = imdir + dark_file
            hdu_list = fits.open(dark_path)
            dark_img = hdu_list[0].data
            hdu_list.close()

            background_img = dark_img*exposure_time + zero_img
            print('Using dark current image: %s' % dark_file)
            print('Using zero image: %s' % zero_file)
        else:
            background_img = zero_img
            print('Using zero image: %s' % zero_file)
        image_data = (image_data_list[idx] - background_img)/flat_img
        print('Using flat image: %s' % flat_file)
        hdu = fits.PrimaryHDU(image_data)

        # Get input file header information
        hdu.header['object'] = image_header_list[idx]['object']
        hdu.header['exptime'] = exposure_time

        try:
            hdu.header['gain'] = image_header_list[idx]['gain']
        except:
            print("Warning: 'GAIN' header keyword not found.")
    
        try:
            hdu.header['rdnoise'] = image_header_list[idx]['rdnoise']
        except:
            print("Warning: 'RDNOISE' header keyword not found.")
            
        try:
            ast_filter = image_header_list[idx]['filter']
            hdu.header['filter'] = ast_filter
        except:
            hdu.header['filter'] = 'none'

        # Update information
        hdu.header['ccdtype'] = 'object'
        utc = dt.datetime.utcnow()
        time_stamp = utc.strftime('%b %d, %Y at %H:%M:%S UTC')
        hdu.header['comment'] = 'Processed ' + time_stamp
        hdulist = fits.HDUList([hdu])
        fout_path = imdir + out_name_list[idx]
        try:
            hdulist.writeto(fout_path)
        except:
            print('\n*** Failed to write %s. ***' % fout_path)
            print('*** The file may already exist. ***')
    return out_name_list
    
def main():
    """Print help if the module is run as a script"""
    print('This is a module for performing CCD data reduction steps.')
    print('Import the module and use python help for documentation.')

if __name__ == '__main__':
    main()