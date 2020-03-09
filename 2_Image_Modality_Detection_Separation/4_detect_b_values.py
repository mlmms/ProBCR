# Detects DWI sequences (as many as the different b values acquired)
# and sb post processing sequences
# and ADC maps
# Creates .csv file with the meta-data fields of interest

import os
import csv
from deid.dicom import get_identifiers

#  source and destin directories
src="/Users/monicasilva/Documents/1_Original_data/DWI_sequences_cohort2_139p"
dst="/Users/monicasilva/Documents/1_Original_data/DWI_sequences_cohort2_139p"

csv_deid_coding = []

listOfDirectories = []
listOfFiles = []
listOfIDs4Files = []
pString = "P"

for root, directories, fileNames in os.walk(src):
    if root == "/Users/monicasilva/Documents/1_Original_data/DWI_sequences_cohort2_139p/P001":
        break
    listOfDirectories = directories
    listOfDirectories.sort()
    #caseNumber = 1

    for directory in listOfDirectories:
        patientDir = directory
        #caseID = pString + str(caseNumber).zfill(3) #  zero padding to the left
        #  csv_deid_coding.append([patientDir, caseID])

        #  walks inside Patient Folder ('rootPatient')
        for rootPatient, directoriesPatient, filenamesImage in os.walk(os.path.join(src, directory)):
            directoriesPatient.sort()
            #for directorySequence in directoriesPatient:   #  for each sequence folder
              #  enters each Sequence folder
                #for rootImage, directoriesImage, filenamesImage in os.walk(os.path.join(src, directory, directorySequence)):

            if len(filenamesImage) > 2:
                #filenamesImage.sort()
                for filename in filenamesImage:

                    possibleFilename = os.path.join(rootPatient, filename)

                    if filename.startswith('I') and filename != "I10":
                        ids = get_identifiers(possibleFilename)   #  deid function: gets identifiers from a dicom file
                        print(possibleFilename)
                        #count = 0
                        for image, fields in ids.items():

                            series_description = fields['SeriesDescription']

                            #if series_description.startswith("sb") or series_description.startswith("b") or series_description == "dadc": #  if is a diffusion
                            #  save these items into .csv
                            patient_id = fields['PatientID']
                            slice_thickness = fields['SliceThickness']
                            repetition_time = fields['RepetitionTime']
                            echo_time = fields['EchoTime']
                            b_value = fields['DiffusionBValue']


                            csv_deid_coding.append([patientDir, series_description, b_value, slice_thickness, repetition_time, echo_time])
                            print(csv_deid_coding)

                            with open(os.path.join(dst, "DWI_b_values_and_parameters2020" + ".csv"), "a") as f:
                                writer = csv.writer(f, delimiter='/')
                                writer.writerow([patientDir, series_description, b_value, slice_thickness, repetition_time, echo_time])
                        #break
                    #break
                #break
            break
