# Detects T2 W Axial sequences
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
#src="/Users/monicasilva/Documents/Deidentification/ProBCR_Deidentified_1_157"
src = "/Users/monicasilva/Documents/1_Original_data/ProBCR_Deidentified_1_157"
dst="/Users/monicasilva/Desktop/"

csv_deid_coding = []

listOfDirectories = []
listOfFiles = []
listOfIDs4Files = []
pString = "P"

createFolder(dst)

for root, directories, fileNames in os.walk(src):
    if root == "Users/monicasilva/Documents/Deidentification/ProBCR_Deidentified_1_157/P001":
        break
    listOfDirectories = directories
    listOfDirectories.sort()

    for directory in listOfDirectories:
        patientDir = directory

        #  walks inside Patient Folder ('rootPatient')
        for rootPatient, directoriesPatient, filenamesPatient in os.walk(os.path.join(src, directory)):
            directoriesPatient.sort()
            for directorySequence in directoriesPatient:   #  for each sequence folder
              #  enters each Sequence folder
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

                                    if "t2" in series_description and "ax" in series_description:

                                        #createFolder(os.path.join(dst, patientDir))

                                        #  save these items into .csv
                                        patient_id = fields['PatientID']
                                        slice_thickness = fields['SliceThickness']
                                        repetition_time = fields['RepetitionTime']
                                        echo_time = fields['EchoTime']
                                        field_strength = fields['MagneticFieldStrength']
                                        pixel_spacing = fields['PixelSpacing']


                                        csv_deid_coding.append([patientDir, series_description, slice_thickness, repetition_time, echo_time, field_strength])
                                        print(csv_deid_coding)

                                        shutil.copy(os.path.join(src,patientDir,directorySequence,filename), os.path.join(dst,patientDir,filename))

                                        with open(os.path.join(dst, "T2W_parameters" + ".csv"), "a") as f:
                                            writer = csv.writer(f, delimiter='/')
                                            writer.writerow([patientDir, series_description, slice_thickness, repetition_time, echo_time, field_strength])
                                    break
                                break
                            break
                        break
                    break