# One of my larger projects using python at university, this project required using the astropy library in order to calibrate images to find Neptune's moons
# This required importing all the images into python where I calibrated the images using the follwing method
# 1. Remove the mean flux value from the overscan region from every pixel in the row 2. Subtract the dark field image to remove CCD artefacts
# 3. Create a master flat by combining the bias and dark subtracted flat fields and normalise this image.
# 4. Flat field each image (subract the master flat from the science images)




from astropy.io import fits
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import mode


master_Flat = []
dark_10_exp_time = fits.getheader("QUESTdata/darks/dark_10.C22.fits")['EXPTIME']
dark_180_exp_time = fits.getheader("QUESTdata/darks/dark_180.C22.fits")['EXPTIME']
dark_10 = fits.getdata("QUESTdata/darks/dark_10.C22.fits")
dark_180 = fits.getdata("QUESTdata/darks/dark_180.C22.fits")


def dark_subtract(dark_data,file_str,Type):
    overscan_width = 40
    shape = np.shape(dark_data)
    interim_Array = np.zeros((shape[0],(shape[1]-overscan_width)))
    file_data  = fits.getdata(file_str)
    for j in range(0,shape[0]):
        #coursework says to use mean, lecture notes say to use median
        #Tried mean and final result too bright - values out of range of original
        #Tried median as per lecture notes and this files look more acceptable
        mean = np.median(file_data[j][-overscan_width:])
        if Type == "Science":
            for l in range(0,shape[1]-overscan_width):
                interim_Array[j,l] = ((file_data[j,l] - mean - dark_data[j,l])/master_Flat[j,l])
        else:
            for l in range(0,shape[1]-overscan_width):
                interim_Array[j,l] = file_data[j,l] - mean - dark_data[j,l]
    
    return interim_Array

def filewriting(data,new_file_name):
    file = fits.PrimaryHDU(data)
    file.header.add_comment(new_file_name)
    fits.writeto(new_file_name,file.data,file.header, overwrite = True)

def getDarkFile(file_exp_time):
    #Get header
    if dark_10_exp_time+5 >= file_exp_time and dark_10_exp_time-5 <= file_exp_time:
        current_dark = dark_10
        print("Dark 10")
    elif dark_180_exp_time+5 >= file_exp_time and dark_180_exp_time-5 <= file_exp_time:
        current_dark = dark_180
        print("Dark 180")
    else:
        print("error")   
    
    return current_dark

def file_calibration(file_shift,new_data):
    calibrated_array = np.roll(new_data,file_shift[0],axis=0)
    calibrated_array = np.roll(calibrated_array,file_shift[1],axis=1)
    
    return calibrated_array



main_Array = []
new_sci = []


flat_file_directory = "QUESTdata/flats/"
flat_files = os.listdir(flat_file_directory)
for i in range(0,len(flat_files)):
    current_file = flat_files[i]
    if (current_file[:1] == ".") or (current_file == "flat.list"):
        print("Exclude file")
    else:
        file_loc = flat_file_directory+current_file
        print(file_loc)
        dark_file_data = getDarkFile(file_exp_time=fits.getheader(file_loc)['EXPTIME'])
        result=dark_subtract(dark_data=dark_file_data, file_str=file_loc, Type="Flat")
        main_Array.append(result)

master_Flat = np.median(main_Array,axis=0)/np.median(np.median(main_Array,axis=0))

#plot histogram
plt.hist(master_Flat)

filewriting(data=master_Flat, new_file_name="masterFlats.fits")

print("Master Flat created!")


sci_file_directory = "QUESTdata/science/"
sci_files = os.listdir(sci_file_directory)
for i in range(0,len(sci_files)):
    current_file = sci_files[i]
    if (current_file[:1] != "."):
        file_loc = sci_file_directory+current_file
        print(file_loc)
        dark_file_data = getDarkFile(file_exp_time=fits.getheader(file_loc)['EXPTIME'])
        new_sci = dark_subtract(dark_data=dark_file_data, file_str=file_loc, Type="Science")
        mode_value = mode(new_sci,axis=None)
        new_sci = new_sci - mode_value[0]

        filewriting(data=new_sci, new_file_name="new "+current_file)
        
print("subtracted science files!")


#new_files = os.listdir()
#print(new_files)
filewriting(data=file_calibration(file_shift=[1317-1291,237-244],new_data=fits.getdata("new 20130911020246s.C22.fits")), new_file_name="calibrated 20130911020246s.C22.fits")
filewriting(data=file_calibration(file_shift=[1317-1259,237-251],new_data=fits.getdata("new 20130911040543s.C22.fits")), new_file_name="calibrated 20130911040543s.C22.fits")
print("calibrated science files!")
 
