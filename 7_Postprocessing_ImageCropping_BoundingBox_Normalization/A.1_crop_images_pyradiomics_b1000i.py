# DWI, b0 radiomics, for multiparametric MRI, for PCa * MMS, 2019
# Analysis I

import os
import radiomics
from radiomics import featureextractor
from radiomics import imageoperations
import SimpleITK as sitk
import glob

# 0_Setting up data

#1st analysis: JUST PROSTATE, no surrounding tissues
src = "/Users/monicasilva/Documents/1_Original_data/Final_Cohort_2020_99p/pre_PyRadiomics/6.2_DWIb1000/Analysis_II(withbackground)/3_BB_input"
dst = "/Users/monicasilva/Documents/1_Original_data/Final_Cohort_2020_99p/pre_PyRadiomics/6.2_DWIb1000/Analysis_II(withbackground)/4_BB_output"

os.chdir(src)

list_images = sorted(glob.glob("*b1000_f.nii"))
cases = list_images
for i in range(len(list_images)):
    cases[i] = list_images[i][0:4]

list_images = sorted(glob.glob("*b1000_f.nii"))
list_labels = sorted(glob.glob("*label.nii"))

for (i, j) in zip(list_images, list_labels):
    case = i[0:4]

    case1 = case + str("_crpBB.nii")
    case2 = case + str("_crpBB_mask.nii")
    imagePath = os.path.join(src, i)
    labelPath = os.path.join(src, j)
    #paramPath = "/Users/monicasilva/PycharmProjects/ProBCR/radiomics_kit/data.yml"

    image = sitk.ReadImage(imagePath)
    label = sitk.ReadImage(labelPath)
    label = sitk.Cast(label, sitk.sitkUInt8)
    #radiomics.setVerbosity(level=10)

    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(image)
    resampler.SetInterpolator(sitk.sitkNearestNeighbor)
    label2 = resampler.Execute(label)

    out = radiomics.imageoperations.checkMask(image, label2)
    boundingBox_mask = out[0]

    out2 = radiomics.imageoperations.cropToTumorMask(image, label2, boundingBox_mask)
    crp_image = out2[0]
    crp_mask = out2[1]
    #OUTPUTS:
    sitk.WriteImage(crp_image, os.path.join(dst, case1))
    sitk.WriteImage(crp_mask, os.path.join(dst, case2))

    #ADC map: no normalization!!!!
    norm_crp_img = radiomics.imageoperations.normalizeImage(crp_image)
    case3 = case + str("_crpBB_norm.nii")
    #OUTPUT:
    sitk.WriteImage(norm_crp_img, os.path.join(dst, case3))

    print("Finished case ", case)
