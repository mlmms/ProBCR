# N4 BIAS FIELD CORRECTION CODE #
import SimpleITK as sitk
import os

path = "ProBCR_T2WI_Data"
output = "T2_BiasFieldCorrected_Cluster"

list = os.listdir(path)
list.sort()

for file in os.listdir(path):
    if file.endswith(".nii"):
        print(file)
        img = sitk.ReadImage(os.path.join(path, file))

        image = sitk.Cast(img, sitk.sitkFloat32)  # Required to have a “real” pixel type: sitkFloat32 or sitkFloat64.
        corrected_img = sitk.N4BiasFieldCorrection(image)

        sitk.WriteImage(corrected_img, os.path.join(output, 'NB_' + file))
        print("Finished N4 Bias Field Correction of ", file)

print("Finished N4 for all .nii files")
