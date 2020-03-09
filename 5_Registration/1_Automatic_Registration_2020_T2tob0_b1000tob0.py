import SimpleITK as sitk
from natsort import natsorted
import os
os.chdir("/Users/monicasilva/Documents/1_Original_data/Final_Cohort_2020_99p/Registration/folder")

# 1. Registers T2 to b0:
    # T2: moving
    # b0: fixed
    # Transform: X

# Get files
b0_b1000_directory = "/Users/monicasilva/Documents/1_Original_data/Final_Cohort_2020_99p/DWI_b0_b1000_w_followup_99p/"
t2ax_directory = "/Users/monicasilva/Documents/1_Original_data/Final_Cohort_2020_99p/T2_w_followup_99p/"
t2axlabel_directory = "/Users/monicasilva/Documents/1_Original_data/Final_Cohort_2020_99p/T2-labels_w_followup_99p/"

DWI_files = os.listdir(b0_b1000_directory)
DWI_files = natsorted(DWI_files)
b0_files = [s for s in DWI_files if s.endswith("b0.nii")] #selects only b0 files
b1000_files = [t for t in DWI_files if t.endswith("b1000.nii")] #selects only b1000 files

for i in range(len(b0_files)):
    idx = b0_files[i]
    b0_files[i] = os.path.join(b0_b1000_directory, idx) # constructs the full path of b0 files

for j in range(len(b1000_files)):
    idx = b1000_files[j]
    b1000_files[j] = os.path.join(b0_b1000_directory, idx)  # constructs the full path of b1000 files

# 0.1 Cropping the Patient Case ID string
total_DWI_cases = len(b0_files)
cases = b0_files[:] #creates a 2nd object, a copy of images' list

for i in range(total_DWI_cases):
    cases[i] = cases[i][95:99]   #substring
print("cases = ",cases)

#movingImage_path = "/Users/monicasilva/Documents/1_Original_data/Final_Cohort_2020_99p/T2_w_followup_99p"
#label_path = "/Users/monicasilva/Documents/1_Original_data/Final_Cohort_2020_99p/T2-labels_w_followup_99p"

#y = 30 # P071
#y = 54 # P101
while y < total_DWI_cases:

    movingImage_path = t2ax_directory
    label_path = t2axlabel_directory

    fixedImage_path = b0_files[y]
    movingImage_path = movingImage_path + "NB_" + cases[y] + "_T2ax.nii"
    label_path = label_path + "NB_" + cases[y] + "_T2ax-label.nii"

    if os.path.isfile(movingImage_path) == False:
        print("break due to ", movingImage_path, "not being a file")
        break

    if os.path.isfile(label_path) == False:
        print("break due to ", label_path, "not being a file")
        break

    fixedImage = sitk.ReadImage(fixedImage_path)
    movingImage = sitk.ReadImage(movingImage_path)
    labelImage = sitk.ReadImage(label_path)
    fixedMask = sitk.BinaryThreshold(fixedImage, 50, 50000, 1, 0)  #
    fixedMask = sitk.Cast(fixedMask, sitk.sitkUInt8)              #
    #movingMask = sitk.BinaryThreshold(movingImage, 50, 50000, 1, 0) #
    #movingMask = sitk.Cast(movingImage, sitk.sitkUInt8)             #


    elastixImageFilter = sitk.ElastixImageFilter()
    elastixImageFilter.SetFixedImage(fixedImage)                    #b0
    elastixImageFilter.SetMovingImage(movingImage)                  #T2
    #elastixImageFilter.SetMovingMask(movingMask)                   #
    elastixImageFilter.SetFixedMask(fixedMask)                     #


    #Defines the parameter map for the registration (Type: rigid)
    rigidparametermap = sitk.GetDefaultParameterMap("rigid")
    rigidparametermap["AutomaticTransformInitialization"] = ["true"]
    rigidparametermap["AutomaticTransformInitializationMethod"] = ["CenterOfGravity"]
    rigidparametermap["OptimizerBase"] = ["Powell"]
    rigidparametermap["ImageSampler"] = ["RandomSparseMask"]        #
    elastixImageFilter.PrintParameterMap(rigidparametermap)
    elastixImageFilter.SetParameterMap(rigidparametermap)

    elastixImageFilter.LogToConsoleOn()
    elastixImageFilter.Execute() #T2 to b0, with X

    # 2. Apply X to T2 LABEL registration to b0
    # Gets the previously computed transform, X,
    # to register T2-Label to b0 (original):
        # T2-label: moving
        # using the Nearest Neighbor interpolator
        # so it is in the same space as b0
    print("step 2")

    #Initialize and define the transform, as the parameter map from the previous ("elastixImageFilter")
    resultTransform = sitk.TransformixImageFilter()
    transformParameterMap = elastixImageFilter.GetTransformParameterMap()
    transformParameterMap[0]["ResampleInterpolator"] = ["FinalNearestNeighborInterpolator"]
    sitk.PrintParameterMap(transformParameterMap)

    resultTransform.SetMovingImage(labelImage) #T2-label
    resultTransform.SetTransformParameterMap(transformParameterMap)

    resultTransform.LogToConsoleOn()
    resultTransform.Execute()

    transformedlabel = resultTransform.GetResultImage()
    #sitk.WriteImage(transformedlabel, "P001_b0-label.nii") #the label that I will use for b0 feature extraction

    # 3. Resampling label-b0 (original)
        # accordingly to the b0 characteristics
        # (reference volume)

    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(fixedImage) #b0
    resampler.SetInterpolator(sitk.sitkNearestNeighbor)
    final_label_b0 = resampler.Execute(transformedlabel)
    filename_final_label_b0 = cases[y] + "_DWI-label.nii"
    sitk.WriteImage(final_label_b0, filename_final_label_b0) # save final, and resampled, DWI label


    # 4. CHECKING IF registration of b0 ----> b1000,
        # through transform X2
        # Improves the ROI targetting in b1000
        # (instead of using the same DWI-label for both b0 and b1000)

    # Get files
    fixedImage_path2 = fixedImage_path #b0
    movingImage_path2 = b1000_files[y] #b1000

    # Read images
    fixedImage2 = sitk.ReadImage(fixedImage_path2) #b0
    movingImage2 = sitk.ReadImage(movingImage_path2) #b1000

    #Initialize new transform X2
    elastixImageFilter2 = sitk.ElastixImageFilter()
    elastixImageFilter2.SetFixedImage(fixedImage2)
    elastixImageFilter2.SetMovingImage(movingImage2)

    #Defines the parameter map of X2 for the registration (Type: translation)
    rigidparametermap2 = sitk.GetDefaultParameterMap("translation")
    rigidparametermap2["AutomaticTransformInitialization"] = ["true"]
    rigidparametermap2["AutomaticTransformInitializationMethod"] = ["CenterOfGravity"]
    rigidparametermap2["OptimizerBase"] = ["Powell"]
    elastixImageFilter2.PrintParameterMap(rigidparametermap2)
    elastixImageFilter2.SetParameterMap(rigidparametermap2)

    elastixImageFilter2.LogToConsoleOn()
    elastixImageFilter2.Execute() #registers b1000 ---> b0, by translation

    b1000_transformedImage = elastixImageFilter2.GetResultImage()
    filename_b1000_transformedImage = cases[y] + "_b1000_f.nii"
    sitk.WriteImage(b1000_transformedImage, filename_b1000_transformedImage)

    print("case = ", cases[y])
    y += 1

    # This way, b0, b0-label (DWI-label) and b1000 will all be in the same space as b0.