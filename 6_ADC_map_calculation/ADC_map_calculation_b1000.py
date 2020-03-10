#This script calculates the ADC map based on 2 DWI image series.
# b0 and b1000

# Given the undefining of the equation for the value of zero,
# we will chose an arbitrary range of say 0.1 to 0.9 is chosen instead
# of the traditional fitting of the data within unity (concept of "normalization).

import glob, os
import SimpleITK as sitk
import math

src = "/Users/monicasilva/Documents/1_Original_data/Final_Cohort_2020_99p/ADC_maps_calculated/input"
os.chdir(src)

b_low_images = glob.glob("*b0.nii")
b_low_images.sort()
print("b_low_images = ", b_low_images)

b_high_images = glob.glob("*b1000_f.nii")
b_high_images.sort()
print("b_high_images = ", b_high_images)

dst_adcs = "/Users/monicasilva/Documents/1_Original_data/Final_Cohort_2020_99p/ADC_maps_calculated"

# Calculate the ADC map

for (i, j) in zip(b_low_images, b_high_images):

    case_i = i
    case_j = j
    i = os.path.join(src, i)  # file path of i
    j = os.path.join(src, j)  # file path of j

    b_low_image = sitk.ReadImage(i)
    b_high_image = sitk.ReadImage(j)
    adc_map = b_low_image
    b_low_image = sitk.Cast(b_low_image, sitk.sitkFloat32)
    b_high_image = sitk.Cast(b_high_image, sitk.sitkFloat32)
    adc_map = sitk.Cast(adc_map, sitk.sitkFloat32)

    width = b_low_image.GetWidth()
    print(b_low_image.GetPixelIDTypeAsString())
    height = b_low_image.GetHeight()
    depth = b_low_image.GetDepth()

    for z in range(depth):

        for x in range(width):

            for y in range(height):

                Sb0 = b_low_image[x, y, z] # intensity of b0_image(x,y,z)
                Sb1 = b_high_image[x, y, z] # intensity of b1000_image(x,y,z)
                #print(Sb0, Sb1)

                if Sb0 <= 0 or Sb1 <= 0:
                    adc_map[x, y, z] = 0
                else:
                    b1 = 1000 # DEFINE b-high-value:
                    adc_map[x, y, z] = - (1/ b1) * (math.log(Sb1) - math.log(Sb0)) * 1e6

    print("z = ", z)
    case = case_i.replace("_b0", "")
    sitk.WriteImage(adc_map, os.path.join(dst_adcs, 'ADC_map_' + case))
    print("Finished ADC map of ", case)
