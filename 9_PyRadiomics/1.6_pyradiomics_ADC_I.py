# ADC radiomics, for multiparametric MRI, for PCa * MMS, 2019
# ANALYSIS I : empty background (only prostate, no surrounding tissue)

import os
import radiomics
from radiomics import featureextractor
import SimpleITK as sitk
import glob
import csv
import statistics

def _writeResults(featureVector):
    global HEADERS, OUTPUTCSV

    with open(OUTPUTCSV, 'a') as outputFile:
        writer = csv.writer(outputFile, lineterminator='\n')
        if HEADERS is None:
            HEADERS = list(featureVector.keys())
            writer.writerow(HEADERS)

        row = []
        for h in HEADERS:
            row.append(featureVector.get(h, "N/A"))
        writer.writerow(row)

def _createParameterFile(h_final):
    import yaml

    data = dict(
        imageType=dict(
            Original={},
            LoG=dict(
                sigma=[1.0, 2.0, 3.0, 4.0, 5.0]
            ),
            Wavelet={},
            Gradient={},
            # LBP2D={},
        ),
        featureClass=dict(
            firstorder=[],
            glcm=[],
            gldm=[],
            glrlm=[],
            glszm=[],
            ngtdm=[],
            shape=[],
        ),
        setting=dict(
            normalize=False,  # !!!!
            binWidth=h_final,  ##
            resampledPixelSpacing=[1.95, 1.95, 7],  ##
            padDistance=5,
            distances=[1, 2, 3],
            correctMask=True,
            force2D=True,
        ),
    )

    with open(paramPath, 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)


# 0_Setting up data
src = "/Users/monicasilva/Documents/1_Original_data/Final_Cohort_2020_99p/pre_PyRadiomics/6.3_ADC/Analysis_I(nobackground)/5_pyRadiomics_input"
dst = "/Users/monicasilva/Documents/1_Original_data/Final_Cohort_2020_99p/pre_PyRadiomics/6.3_ADC/Analysis_I(nobackground)/6_pyRadiomics_output"
OUTPUTCSV = os.path.join(dst,"results_2020_adc_I_bin_method_99p.csv")
HEADERS = None
os.chdir(src)

ls_max = []
ls_min = []
ls_int_range = []
ls_h = []
calculated_bins = []

paramPath = "/Users/monicasilva/PycharmProjects/ProBCR/radiomics_kit/data_ADC_I.yml"

list_images = sorted(glob.glob("*BB.nii"))
cases = list_images
for i in range(len(list_images)):
    cases[i] = list_images[i][0:4]

list_images = sorted(glob.glob("*BB.nii"))

list_masks = sorted(glob.glob("*mask.nii"))

#For each image (and mask), calculates the (max,min) and saves them in a list.
for (i, j) in zip(list_images, list_masks):

    case = i[0:4]
    imagePath = os.path.join(src, i)
    labelPath = os.path.join(src, j)

    image = sitk.ReadImage(imagePath)
    label = sitk.ReadImage(labelPath)

    # ADDED
    filter = sitk.MinimumMaximumImageFilter()
    filter.Execute(image)
    ImageHighestIntensity = filter.GetMaximum()
    ImageLowestIntensity = filter.GetMinimum()

    int_range = (ImageHighestIntensity - ImageLowestIntensity)

    ls_max.append(ImageHighestIntensity)
    ls_min.append(ImageLowestIntensity)
    ls_int_range.append(int_range)
    h = (int_range / 80)
    ls_h.append(h)

    #out_bin = radiomics.imageoperations.getBinEdges(image,binWidth=h)
    #calculated_bins.append(len(out_bin) + 1)

h_mean = statistics.mean(ls_h)

print("h mean = ", h_mean)
print("round h = ", round(h_mean, 7))


h_final = round(h_mean, 7)

"""
#Auxiliary: Loop to check how many bins is this h_final is producing
for i in list_images:

    case = i[0:4]
    imagePath = os.path.join(src,i)

    image = sitk.ReadImage(imagePath)

    out_bin = radiomics.imageoperations.getBinEdges(image,binWidth=h_final)
    calculated_bins.append(len(out_bin) + 1)

print("calculated_bins = ", calculated_bins)
print("mean calculated_bins = ", statistics.mean(calculated_bins))
"""

#FINALLY: give the value of the calculated h_final to the parameter file
_createParameterFile(h_final)

#Perform pyRadiomics feature extraction with such parameter
for (i, j) in zip(list_images, list_masks):

    case = i[0:4]
    imagePath = os.path.join(src, i)
    labelPath = os.path.join(src, j)

    image = sitk.ReadImage(imagePath)
    label = sitk.ReadImage(labelPath)

    #print(len(out_bin)) #number of bins
    radiomics.setVerbosity(level=10)

    extractor = radiomics.featureextractor.RadiomicsFeatureExtractor(paramPath)

    print("Enabled filters:\n\t", extractor.enabledImagetypes, "")
    print("Extraction parameters:\n\t", extractor.settings)
    print("Enabled features:\n\t", extractor.enabledFeatures)

    # Extract features
    print("Calculating features")
    featureVector = extractor.execute(image, label)

    d1 = featureVector
    d1.update({'Case': case})
    d1.move_to_end('Case', last=False)

    _writeResults(d1)
    print("ty,next", case)




    #for featureName in featureVector.keys():
    #    print("Computed %s: %s" % (featureName, featureVector[featureName]))

