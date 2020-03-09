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
src = "/Users/monicasilva/Documents/1_Original_data/AllImageTypes_ProBCR_Deidentified_1_157_cohort2_139p"
dst="/Users/monicasilva/Documents/1_Original_data/Parameter_gathering"

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
                                                                repetition_time, echo_time, flip_angle, manufacturer,
                                                                manufacturer_model_name,
                                                                protocol_name, receive_coil_name, scanning_sequence,
                                                                sequence_variant,
                                                                rows, columns, pixel_spacing, slice_thickness,
                                                                spacing_between_slices])
                                        print(csv_deid_coding)

                                        #shutil.copy(os.path.join(src,patientDir,directorySequence,filename), os.path.join(dst,patientDir,filename))

                                        with open(os.path.join(dst, "T2_parameters_2020" + ".csv"), "a") as f:
                                            writer = csv.writer(f, delimiter='/')
                                            writer.writerow([patientDir, series_description, magnetic_field, b_value,
                                                             repetition_time, echo_time, flip_angle, manufacturer,
                                                             manufacturer_model_name,
                                                             protocol_name, receive_coil_name, scanning_sequence,
                                                             sequence_variant,
                                                             rows, columns, pixel_spacing, slice_thickness,
                                                             spacing_between_slices])
                                    #patientDir, series_description, magnetic_field, b_value, repetition_time, echo_time, flip_angle, manufacturer,manufacturer_model_name,protocol_name, receive_coil_name, scanning_sequence,sequence_variant,rows, columns, pixel_spacing, slice_thickness,spacing_between_slices
                                    break
                                break
                            break
                        break
                    break