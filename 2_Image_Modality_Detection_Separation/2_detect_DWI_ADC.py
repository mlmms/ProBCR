# Detects DWI sequences (not sb000 or sbXXXX)
# and ADC maps
# Creates .csv file with the meta-data fields of interest

import os
import csv
from deid.dicom import get_identifiers

#  source and destin directories
#src ="/Users/monicasilva/Documents/Deidentification/Test_folder"

#Internal
src="/Users/monicasilva/Documents/Deidentification/ProBCR_Deidentified_1_157"
dst="/Users/monicasilva/Documents/Deidentification/ProBCR_DWI_Data"

csv_deid_coding = []

listOfDirectories = []
listOfFiles = []
listOfIDs4Files = []
pString = "P"

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

                                    if series_description.startswith("DIF") or series_description.startswith("Dif") or "DWI" in series_description or series_description.startswith("dD") or series_description == "dadc" or series_description.startswith("dADC") or series_description == "b1400 4mm" or series_description == "b1500" or series_description.startswith("db") or series_description == "sadc" or series_description == "b1400": #  if is a diffusion

                                        #  save these items into .csv
                                        patient_id = fields['PatientID']
                                        slice_thickness = fields['SliceThickness']
                                        repetition_time = fields['RepetitionTime']
                                        echo_time = fields['EchoTime']
                                        field_strength = fields['MagneticFieldStrength']


                                        csv_deid_coding.append([patientDir, series_description, slice_thickness, repetition_time, echo_time, field_strength])
                                        print(csv_deid_coding)

                                        with open(os.path.join(dst, "DWI_ADC_parameters_field" + ".csv"), "a") as f:
                                            writer = csv.writer(f, delimiter='/')
                                            writer.writerow([patientDir, series_description, slice_thickness, repetition_time, echo_time, field_strength])
                                    break
                                break
                            break
                        break
                    break
