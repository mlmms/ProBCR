# Detects DWI sequences of a specific value
# copies entire tree
# Creates .csv file with the meta-data fields of interest

import os
import shutil
import csv
from deid.dicom import get_identifiers

#  source and destin directories
#src="/Users/monicasilva/Documents/2_Preprocessing/DWI"
src ="/Users/monicasilva/Documents/1_Original_data/DWI_sequences_cohort2_139p"
dst="/Users/monicasilva/Documents/1_Original_data/DWI_b1500_folders_separatedPycharm_script1_2020"

csv_deid_coding = []

listOfDirectories = []
listOfFiles = []
listOfIDs4Files = []
pString = "P"

stop = False

for root, directories, fileNames in os.walk(src):
    if root == "/Users/monicasilva/Documents/1_Original_data/DWI_sequences_cohort2_139p/P001":
        break
    listOfDirectories = directories
    listOfDirectories.sort()
    #caseNumber = 1

    for directory in listOfDirectories:
        patientDir = directory
        stop = False

        for rootPatient, directoriesPatient, filenamesImage in os.walk(os.path.join(src, directory)):
            filenamesImage.sort()
            #for directorySequence in directoriesPatient:   #  for each sequence folder
              #  enters each Sequence folder
            for rootImage, directoriesImage, filenamesImage in os.walk(os.path.join(src, directory)):
                if stop is True:
                    break

                if len(filenamesImage) > 2:
                    #filenamesImage.sort()
                    for filename in filenamesImage:
                        if stop is True:
                            break

                        possibleFilename = os.path.join(rootImage, filename)

                        if filename.startswith('I') and filename != "I10":
                            ids = get_identifiers(possibleFilename)   #  deid function: gets identifiers from a dicom file
                            print(possibleFilename)
                            #count = 0
                            for image, fields in ids.items():

                                b_value = fields['DiffusionBValue']

                                if b_value == "1000.0":
                                    shutil.copytree(os.path.join(src, directory), os.path.join(dst,patientDir))
                                    stop = True
                                    break

                            #break
                        #break
                    #break
                #break