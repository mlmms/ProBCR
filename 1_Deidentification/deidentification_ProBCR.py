import os
import csv
from deid.dicom import get_identifiers, replace_identifiers
from deid.config import DeidRecipe

recipe = DeidRecipe(deid='./deidProBCR.dicom')

def createFolder(folder):
    try:
        if not os.path.exists(folder):
            os.makedirs(folder)
    except OSError:
        print ('Error: Creating folder. ' + folder)

#source and destin directories
src="/Volumes/TOSHIBA_EXT/ProBCR"
dst="/Volumes/TOSHIBA_EXT/ProBCR_Deid"


createFolder(dst)

csv_deid_coding = []

listOfDirectories = []
listOfFiles = []
listOfIDs4Files = []
pString = "P"

for root, directories, fileNames in os.walk(src):
    if root == "/Volumes/TOSHIBA_EXT/ProBCR/S50710":
        break
    listOfDirectories = directories
    caseNumber = 000

    for directory in listOfDirectories:
        patientDir = directory
        caseID = pString + str(caseNumber).zfill(3) #zero padding to the left
        #csv_deid_coding.append([patientDir, caseID])
        createFolder(os.path.join(dst, caseID))
        #walks inside Patient Folder ('rootPatient')

        for rootPatient, directoriesPatient, filenamesPatient in os.walk(os.path.join(src, directory)):

            for directorySequence in directoriesPatient: #for each sequence folder
                createFolder(os.path.join(dst, caseID, directorySequence))
                #enters each Sequence folder
                for rootImage, directoriesImage, filenamesImage in os.walk(os.path.join(src, directory, directorySequence)):

                    for filename in filenamesImage:
                        possibleFilename = os.path.join(rootImage, filename)
                        if filename.startswith('I'):
                        #if 'DIRFILE' not in possibleFilename and 'dirty' not in possibleFilename and '.DS_Store' not in possibleFilename and '.bmp' not in possibleFilename and '._I10' not in possibleFilename and '._I11' not in possibleFilename and '._I00' not in possibleFilename and '._I200' not in possibleFilename and '._I460' not in possibleFilename and '._I880' not in possibleFilename:
                            #print(possibleFilename)
                            ids = get_identifiers(possibleFilename) #deid function to get identifiers from a dicom file

                            recipe.deid #changing header values

                            #print(recipe.deid)
                            #print(recipe.get_actions()) # check the actions that are defined

                            updated_ids = dict()
                            count = 0
                            for image, fields in ids.items():
                                #save these items to put into .csv
                                patientName = fields['PatientName']
                                study_date = fields['StudyDate']
                                #institution_name = fields['InstitutionName']
                                #patient_age = fields['PatientAge']
                                patient_birthdate = fields['PatientBirthDate']
                                patient_id = fields['PatientID']
                                #patient_weight = fields['PatientWeight']
                                #performing_physician_name = fields['PerformingPhysicianName']
                                #requesting_physician = fields['RequestingPhysician']
                                #scheduled_performingphysicianname = fields['ScheduledPerformingPhysicianName']

                                #replace the items
                                fields['patient_name'] = caseID
                                fields['patient_id'] = caseID
                                fields['std_number'] = "000000"
                                fields['std_date'] = "20010101"
                                fields['std_name'] = "AaaaaBbbbbb"


                                updated_ids[image] = fields
                                count += 1


                            # use the deid recipe and updated to create new files
                            cleaned_files = replace_identifiers(possibleFilename,
                                                                deid=recipe,
                                                                ids=updated_ids,
                                                                output_folder=os.path.join(dst,caseID,directorySequence))

                            ids_cleaned_files = get_identifiers(cleaned_files)
                            #print(ids_cleaned_files)



                    #os.rename(rootPatient, os.path.join(root, caseID))
        caseNumber = caseNumber + 1
        #print(caseNumber)

        csv_deid_coding.append([patientName, patient_id, patientDir, caseID, patient_birthdate, study_date])


        #print(csv_deid_coding)

        with open(os.path.join(dst, "DeID_Coding_BCR_Prostate" + ".csv"), "a") as f:
            writer = csv.writer(f)
            writer.writerows([patientName, patient_id, patientDir, caseID, patient_birthdate, study_date])
