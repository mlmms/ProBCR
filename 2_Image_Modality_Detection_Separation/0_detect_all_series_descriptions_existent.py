# Detects DWI sequences
# Copies them into another folder
# Creates .csv file with the meta-data fields of interest

import os
import csv
from deid.dicom import get_identifiers
import shutil

def createFolder(folder):
    try:
        if not os.path.exists(folder):
            os.makedirs(folder)
    except OSError:
        print ('Error: Creating folder. ' + folder)

#  source and destin directories
#src ="/Users/monicasilva/Documents/Deidentification/Test_folder"

#Internal
src="/Users/monicasilva/Documents/1_Original_data/AllImageTypes_ProBCR_Deidentified_1_157_cohort2_139p"
dst="/Users/monicasilva/Documents/1_Original_data/All_possible_sequences_cohort2"

csv_deid_coding = []

listOfDirectories = []
listOfFiles = []
listOfIDs4Files = []
pString = "P"

createFolder(dst)

for root, directories, fileNames in os.walk(src):
    if root == "/Users/monicasilva/Documents/1_Original_data/AllImageTypes_ProBCR_Deidentified_1_157_cohort2_139p/P001":
        break
    listOfDirectories = directories
    listOfDirectories.sort()

    for directory in listOfDirectories:
        patientDir = directory

        for rootPatient, directoriesPatient, filenamesPatient in os.walk(os.path.join(src, directory)):
            directoriesPatient.sort()
            count_patientDir = 1
            for directorySequence in directoriesPatient:
                for rootImage, directoriesImage, filenamesImage in os.walk(os.path.join(src, directory, directorySequence)):
                    if len(filenamesImage) > 2:
                        for filename in filenamesImage:
                            possibleFilename = os.path.join(rootImage, filename)

                            if filename.startswith('I'):
                                ids = get_identifiers(possibleFilename)   #  deid function: gets identifiers from a dicom file
                                print(possibleFilename)

                                for image, fields in ids.items():

                                    series_description = fields['SeriesDescription']
                                    series_description = series_description.lower()

                                    #count_patientDir = 1

                                    #if series_description.startswith("dif") or series_description.startswith("ddif") or "dwi" in series_description or series_description.startswith("b1"):
                                    aux_str = patientDir + "_" + str(count_patientDir)
                                    #createFolder(os.path.join(dst, aux_str))


                                    #  save these items into .csv
                                    # patient_id = fields['PatientID']
                                    # slice_thickness = fields['SliceThickness']
                                    # repetition_time = fields['RepetitionTime']
                                    # echo_time = fields['EchoTime']
                                    # field_strength = fields['MagneticFieldStrength']
                                    # pixel_spacing = fields['PixelSpacing']
                                    # b_value = fields['DiffusionBValue']


                                    csv_deid_coding.append([aux_str, series_description])
                                    print(csv_deid_coding)
                                    shutil.copytree(os.path.join(src,patientDir,directorySequence), os.path.join(dst,aux_str))

                                    count_patientDir += 1

                                    with open(os.path.join(dst, "All_series_existent" + ".csv"), "a") as f:
                                        writer = csv.writer(f, delimiter='/')
                                        writer.writerow([aux_str, series_description])

                                break
                            break
                        break
                    break
