b0b1000_files = "/Users/monicasilva/Documents/1_Original_data/DWI_b0_b1000_stack_separation_2020"

# Detects DWI sequences (as many as the different b values acquired)
# and sb post processing sequences
# and ADC maps
# Creates .csv file with the meta-data fields of interest

import os
from natsort import natsorted
import csv
from deid.dicom import get_identifiers
import pandas as pd

#  source and destin directories
src=b0b1000_files
dst="/Users/monicasilva/Documents/1_Original_data/Parameter_gathering"

####################################################
#Auxiliary: gather UID of final cohort
final_cohort = "/Users/monicasilva/Documents/1_Original_data/Final_Cohort_2020_99p/Registration"
uids_final_cohort = os.listdir(final_cohort)
uids_final_cohort = natsorted(uids_final_cohort)
df = pd.DataFrame(uids_final_cohort)
# saving the dataframe
df.to_csv(os.path.join(dst, "uids_final_cohort.csv"), header=False, index=False)
####################################################


csv_deid_coding = []

listOfDirectories = []
listOfFiles = []
listOfIDs4Files = []
pString = "P"

for root, directories, fileNames in os.walk(src):
    if root == os.path.join(src, "P001_2"):
        break
    listOfDirectories = directories
    listOfDirectories.sort()
    #caseNumber = 1

    for directory in listOfDirectories:
        patientDir = directory

        for rootPatient, directoriesPatient, filenamesImage in os.walk(os.path.join(src, directory)):
            directoriesPatient.sort()
            for directorySequence in directoriesPatient:   # for each sequence folder
                # enters each Sequence folder
                for rootImage, directoriesImage, filenamesImage in os.walk(os.path.join(src, directory, directorySequence)):

                    if len(filenamesImage) > 2:
                        #filenamesImage.sort()
                        for filename in filenamesImage:

                            possibleFilename = os.path.join(rootPatient, directorySequence, filename)

                            if filename.startswith('I') and filename != "I10":
                                ids = get_identifiers(possibleFilename)   #  deid function: gets identifiers from a dicom file
                                print(possibleFilename)
                                #count = 0
                                for image, fields in ids.items():

                                    #  save these items into .csv
                                    patient_id = fields['PatientID']
                                    series_description = fields['SeriesDescription']

                                    magnetic_field = fields['MagneticFieldStrength']
                                    b_value = fields['DiffusionBValue']
                                    repetition_time = fields['RepetitionTime']
                                    echo_time = fields['EchoTime']
                                    flip_angle = fields['FlipAngle']
                                    manufacturer = fields['Manufacturer']
                                    manufacturer_model_name = fields['ManufacturerModelName']
                                    protocol_name = fields['ProtocolName']
                                    receive_coil_name = fields['ReceiveCoilName']
                                    scanning_sequence = fields['ScanningSequence']
                                    sequence_variant = fields['SequenceVariant']

                                    rows = fields['Rows']
                                    columns = fields['Columns']
                                    pixel_spacing = fields['PixelSpacing']
                                    slice_thickness = fields['SliceThickness']
                                    spacing_between_slices = fields['SpacingBetweenSlices']


                                    csv_deid_coding.append([patientDir, series_description, magnetic_field, b_value,
                                                            repetition_time, echo_time, flip_angle, manufacturer, manufacturer_model_name,
                                                            protocol_name, receive_coil_name, scanning_sequence, sequence_variant,
                                                            rows, columns, pixel_spacing, slice_thickness, spacing_between_slices])

                                    print(csv_deid_coding)

                                    with open(os.path.join(dst, "DWI_parameters_2020" + ".csv"), "a") as f:
                                        writer = csv.writer(f, delimiter='/')
                                        writer.writerow([patientDir, series_description, magnetic_field, b_value,
                                                        repetition_time, echo_time, flip_angle, manufacturer, manufacturer_model_name,
                                                        protocol_name, receive_coil_name, scanning_sequence, sequence_variant,
                                                        rows, columns, pixel_spacing, slice_thickness, spacing_between_slices])
                                break
                            #break
                        #break
                    break