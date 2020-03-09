# N4 BIAS FIELD CORRECTION #
import SimpleITK as sitk
import os

path = "/Users/monicasilva/Documents/2_Preprocessing/2.1_Preprocessing_T2ax/T2ax_nifty/T2ax_nii_files"
output = "/Users/monicasilva/Documents/2_Preprocessing/2.1_Preprocessing_T2ax/T2ax_nifty/T2ax_BiasFieldCorrected"

list = os.listdir(path)
list.sort()

for file in list:
    if file.endswith("033_T2ax.nii"):
        print(file)
        img = sitk.ReadImage(os.path.join(path, file))

        image = sitk.Cast(img, sitk.sitkFloat32)  # Required to have a “real” pixel type: sitkFloat32 or sitkFloat64.
        corrected_img = sitk.N4BiasFieldCorrection(image)

        sitk.WriteImage(corrected_img, os.path.join(output, 'NB_' + file))
        print("Finished N4 Bias Field Correction of ", file)

print("Finished N4 for all .nii files")
