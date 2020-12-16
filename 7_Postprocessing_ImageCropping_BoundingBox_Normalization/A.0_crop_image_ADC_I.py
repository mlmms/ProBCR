# Aim: Crop image for normalization purposes
# Status: properly working
# MMS 2019

#INPUT:
# 1) original image
# 2) original label
#OUTPUT:
# 1) resampled label (reference: image)
# 2) cropped image (result of ' Image * Resampled Label ' )

# The resample of the Label is more important in the case of DWI images and ADC, where the Labels were obtained after
# the Registration process. In the T2ax case, it would not be necessary, because the labels were created upon the original T2 images

import glob, os
import numpy as np
import SimpleITK as sitk
from shutil import copyfile
#from numpy import * #all functions will be loaded into the local namespace

multiply=sitk.MultiplyImageFilter()

src = "/Users/monicasilva/Documents/1_Original_data/Final_Cohort_2020_99p/preRadiomics/ADC*Mask_i/input"
os.chdir(src)

images = glob.glob("ADC*.nii")
images.sort()
print("images = ", images)

masks = glob.glob("*label.nii")
masks.sort()
print("masks = ", masks)

dst = "/Users/monicasilva/Documents/1_Original_data/Final_Cohort_2020_99p/preRadiomics/ADC*Mask_i"

### Reading images to numpy arrays ###
for (i, j) in zip(images, masks):

    a = i
    f = j
    i = os.path.join(src, i)
    j = os.path.join(src, j)
    img = sitk.ReadImage(i)
    mask = sitk.ReadImage(j)

    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(img)
    resampler.SetInterpolator(sitk.sitkNearestNeighbor)
    mask = resampler.Execute(mask)
    filename_resampled_mask = "res_" + f

    #OUTPUT 1)
    sitk.WriteImage(mask, os.path.join(dst, filename_resampled_mask))  # save final, and resampled, (DWI) label
    print("Finished resampling of ", a)


    print(i)
    print(img.GetPixelIDTypeAsString())
    print("Image size: " + str(img.GetSize()))
    print("Original Spacing: ", img.GetSpacing())
    print("Original Origin: ", img.GetOrigin())
    print("Direction: ", img.GetDirection())

    print(j)
    print(mask.GetPixelIDTypeAsString())
    print("Image size: " + str(mask.GetSize()))
    print("Original Spacing: ", mask.GetSpacing())
    print("Original Origin: ", mask.GetOrigin())
    print("Direction: ", mask.GetDirection())

    img = sitk.Cast(img, sitk.sitkFloat64)
    mask = sitk.Cast(mask, sitk.sitkFloat64)

    # OUTPUT 2):
    crp_img = multiply.Execute(img, mask)

    print(crp_img.GetPixelIDTypeAsString())
    print("Image size: " + str(crp_img.GetSize()))
    print("Original Spacing: ", crp_img.GetSpacing())
    print("Original Origin: ", crp_img.GetOrigin())
    print("Direction: ", crp_img.GetDirection())


    sitk.WriteImage(crp_img, os.path.join(dst, 'crp' + a))
    #copyfile(os.path.join(src,f), os.path.join(dst,f))
    print("Finished cropping of ", a)




# Crop a bit outside the segmentation (rectangular)
# 1. Calculate mask dimensions (a x b)
# 2. Create rectangle (a x b)
# 3. Multiply image with rectangle
