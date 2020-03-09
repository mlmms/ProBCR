# N4 BIAS FIELD CORRECTION CODE #
import SimpleITK as sitk

# import matplotlib
# matplotlib.use("TkAgg")
# from matplotlib import pyplot as plt

# Read image
# img_file = '/Users/monicasilva/Documents/Deidentification/ProBCR_T2WI_Data/P003.nii'
img_file = '/Users/monicasilva/Documents/2_Preprocessing/T2ax_nii_files/P060_T2ax.nii'
img = sitk.ReadImage(img_file)

# Get image dimensions
print(img.GetPixelIDTypeAsString())
print("Image size: " + str(img.GetSize()))
print("Original Spacing: ", img.GetSpacing())
print("Original Origin: ", img.GetOrigin())
print("Direction: ", img.GetDirection())

# Required to have a “real” pixel type: sitkFloat32 or sitkFloat64.
image = sitk.Cast(img, sitk.sitkFloat32)

corrected_img = sitk.N4BiasFieldCorrection(image)

sitk.WriteImage(corrected_img, "NB_P060_T2ax.nii")
print("Finished N4 Bias Field Correction...")
