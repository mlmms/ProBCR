# Detects DWI sequences (as many as the different b values acquired)
# and sb post processing sequences
# and ADC maps
# Creates .csv file with the meta-data fields of interest

import os
import shutil
from deid.dicom import get_identifiers

#  source and destin directories
src="/Users/monicasilva/Documents/1_Original_data/DWI_b1000_folders_separatedPycharm_script1_2020"
dst="/Users/monicasilva/Documents/1_Original_data/DWI_b0_b1000_stack_separation_2020"

listOfDirectories = []
listOfFiles = []
listOfIDs4Files = []
pString = "P"

stop = False

for root, directories, fileNames in os.walk(src):
    if root == "/Users/monicasilva/Documents/1_Original_data/DWI_b1000_folders_separatedPycharm_script1_2020/P001_2":
        break
    listOfDirectories = directories
    listOfDirectories.sort()

    for directory in listOfDirectories:
        directorySequence = directory
        stop = False

        #for rootPatient, directoriesPatient, filenamesPatient in os.walk(os.path.join(src, directory)):
        #   directoriesPatient.sort()
        #for directorySequence in directoriesPatient:   #  for each sequence folder
          #  enters each Sequence folder
        string1000 = directorySequence + "_b1000"
        string0 = directorySequence + "_b0"
        folder_directory1000 = os.path.join(dst, directorySequence, string1000)
        folder_directory0 = os.path.join(dst, directorySequence, string0)
        os.makedirs(folder_directory1000, exist_ok=True)
        os.makedirs(folder_directory0, exist_ok=True)


        for rootImage, directoriesImage, filenamesImage in os.walk(os.path.join(src, directorySequence)):
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

                        for image, fields in ids.items():

                            b_value = fields['DiffusionBValue']

                            if b_value == "1000.0":
                                shutil.copy(os.path.join(src, directorySequence, filename), folder_directory1000)

                            if b_value == "0.0":
                                shutil.copy(os.path.join(src, directorySequence, filename), folder_directory0)
                            else:
                                break