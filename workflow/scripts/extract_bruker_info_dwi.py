import sys
import pprint
import json
import pydicom
from pybruker import jcamp
import numpy as np
from glob import glob


# adapted from conversion scripts written by Naila Rahman (Corey Baron Lab)

#glob niftis to figure out series number
nii = sorted(glob( snakemake.input.nii_folder + f'/*cfmmDtiEpi_10B0_30B1k_60B2k_Marm_ISO400.nii'))[-1] #get last one only



series_num = nii.split('/')[-1].split('_')[0]


dcm = glob(snakemake.input.dcm_folder + f'/*/*/*/*/*/{series_num}/*.dcm')[-1]



# read dicom header
H = pydicom.read_file(dcm, stop_before_pixels=True)


# read bruker headers
method = jcamp.jcamp_parse(
    H[0x0177, 0x1100].value.decode('utf-8').splitlines()
)
#visu_pars = jcamp.jcamp_parse(
#    H[0x0177, 0x1101].value.decode('utf-8').splitlines()
#)

with open(snakemake.output.method_json, 'w') as fp:
    json.dump(method,  fp, indent=4)

# Bvalue information
bval = method["$PVM_DwEffBval"]["value"]
bvec = method["$PVM_DwDir"]["value"]
bvec = np.reshape(bvec, (int(len(bvec)/3+0.5), 3)).T
bmat = method["$PVM_DwBMat"]["value"]
bmat = np.reshape(bmat, (int(len(bmat)/9+0.5), 9)).T

with open(snakemake.output.bval, 'w') as fp:
    for item in bval:
        fp.write("%s " % item)
    fp.write("\n")

with open(snakemake.output.bmat, 'w') as fp:
    for row in bmat:
        for item in row:
            fp.write("%s " % item)
        fp.write("\n")

# Determine bvec from bmat
bmat = np.transpose(bmat)
bvecFromMat = np.zeros((len(bval), 3))
for idx, row in enumerate(bmat):
    row = np.reshape(row, (3, 3))
    u, s, vh = np.linalg.svd(row, full_matrices=True)
    bvecFromMat[idx] = vh[0]

# Get polarity of bvec correct based on input vectors (since polarity is arbitrary after svd)
bvec = np.transpose(bvec)
if len(bval) > len(bvec):
    # Fill in b0 acquisitions that were not in input dir vector
    for n in range(len(bval) - len(bvec)):
        idx = np.argsort(bval)
        if len(bval) != 55:                                       #if not OGSE -- not needed for OGSE since b0 is included in dir vector
            bvec = np.insert(bvec, idx[0], 0, axis=0)
for idx, dir_n in enumerate(bvec):
    mind1 = np.argmax(np.abs(dir_n))
    mind2 = np.argmax(np.abs(bvecFromMat[idx]))
    fact = np.sign(dir_n[mind1]*bvecFromMat[idx][mind2])
    if fact < 0:
        bvecFromMat[idx] = fact*bvecFromMat[idx]

bvecFromMat = np.transpose(bvecFromMat)
with open(snakemake.output.bvec, 'w') as fp:
    for row in bvecFromMat:
        for item in row:
            fp.write("%s " % item)
        fp.write("\n")
