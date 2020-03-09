import os
import dicom2nifti
from natsort import natsorted

#  source and destin directories
src="/Users/monicasilva/Documents/2_Preprocessing/2.1_Preprocessing_T2ax/T2ax_dicom"
dst="/Users/monicasilva/Documents/2_Preprocessing/2.1_Preprocessing_T2ax/T2ax_nifty/T2ax_nii_files"

listOfDirectories = []
listOfFiles = []
listOfIDs4Files = []
pString = "P"


for root, directories, fileNames in os.walk(src):
    if root == "/Users/monicasilva/Documents/2_Preprocessing/2.1_Preprocessing_T2ax/T2ax_dicom/P001":
        break
    listOfDirectories = directories
    listOfDirectories.sort()

    for directory in listOfDirectories:
        patientDir = directory
        a = os.walk(os.path.join(src, directory))
        for rootPatient, directoriesPatient, filenamesPatient in os.walk(os.path.join(src, directory)):
            directoriesPatient.sort()
            aux = 0
            for directorySequence in directoriesPatient:   #  for each sequence folder

                string = str(directorySequence) + ".nii"
                dicom2nifti.dicom_series_to_nifti(os.path.join(src, directory, directorySequence), os.path.join(dst, string), reorient_nifti=True)
                aux += 1
            if aux == 2:
                break