# conmat
Connectivity Matrices

+ Camino
+ FSL
+ Matlab
+ Paraview

## Pre-Processing
- fsl2scheme
- flirt
- convert_xfm
- fslmaths

## Deterministic Tractography
+ RK4
+ track
- input models
- dt, diffusion tensor data produced by modelfit
- tracking algorithms
- rk4 vs. euler

- interpolation algorithms 
- linear vs tend vs nn

Parameters determined from: Reproducibility of graph metrics of human brain structural networks

wdtfit - wrapper to modelfit
track
fa


## Probabilistic Tractography
+ probability look-up table
+ PICo

- modelfit
- estimatesnr
- dtlutgen
- picopdfs
- track 


## Connectivity Matrix
- conmat

## Visualisation Tools
- vtkstreamlines
- pdview 
1. starting in archive/data
2. input: projectName
3. go into ./${projectName}/data/nii/
4. produce a list from the folder names
FOR subject IN subject_list
5. go into cd ../../
6. go into ./pipelines/dtifit/${subject}
## this is the directory you want to be in to get the files

copy the entire folder to /tmp
- enter the folder
- make a second copy of the 
bval, bvec - rename with no extension
*_eddy_correct_b0_bet_mask.nii.gz - rename to nodif_brain_mask.nii.gz
*_eddy_correct.nii.gz - rename to data.nii.gz
run bedpostx in folder - bedpostx ./
copy the output *.bedpostx to targetPath = '/archive/data-2.0/SPINS/pipelines/bedpostX/'
delete the entire subject folder in tmp


** starting from directory you want to get the files
copy the entire folder to /tmp
- enter the folder
- run registrations - make a registration .py

- run streamlining using bedpostx data, and save as .Bfloat
- run conmat in folder
copy the output *_sc.csv to targetPath = '/archive/data-2.0/SPINS/pipelines/bedpostX/conmat/', renamed as subjectName_connectivity.csv
delete the entire subject folder in tmp





location of MNI and shen
/projects/lliu/ImportantFiles/

deterministic or probabilistic?
if deterministic, run deterministic.sh



##SETUP
. $MODULESHOME/init/bash
module load FSL
module load camino

##INPUTS
bvecFile='*.bvec'
bvalFile='*.bval'
multiVol='*_correct.nii.gz'
b0Mask='*_b0_bet_mask.nii.gz'
b0Brain='*_b0_bet.nii.gz'
T1wBrain='/archive/data/' + projectName + '/pipelines/hcp/' + subject + '/T1w/T1w_brain.nii.gz'
wmParc='/archive/data/' + projectName + '/pipelines/hcp/' + subject + '/T1w/wmparc.nii.gz'
mniBrain='/projects/lliu/ImportantFiles/MNI152_T1_2mm_brain.nii.gz'
shenParc='/projects/lliu/ImportantFiles/shen_2mm_268_parcellation.nii.gz'

#get T1w_brain, wmparc, shen, MNI

echo "making scheme file"
fsl2scheme \
-bvecfile ${bvecFile} \
-bvalfile ${bvalFile} > bVectorScheme.scheme
echo

echo "registering b0-T1"
flirt \
-in ${b0Brain} \
-ref ${T1wBrain} \
-omat tfm.mat

convert_xfm \
-omat tfm_invert.mat \
-inverse tfm.mat
echo

echo "registering white matter mask"
flirt \
-in ${wmParc} \
-ref ${b0Brain} \
-applyxfm \
-init tfm_invert.mat \
-o wmparc_invert.nii.gz
echo
fslmaths \
wmparc_invert.nii.gz \
-thr 2500 \
-bin wmparc_invert_bin.nii.gz
echo

echo "registering MNI"
flirt \
-in ${mniBrain} \
-ref ${b0Brain} \
-interp nearestneighbour \
-omat mni.mat
echo

echo "registering atlas"
flirt \
-in ${shenParc} \
-ref ${b0Brain} \
-interp nearestneighbour \
-applyxfm \
-init mni.mat \
-o atlas.nii.gz
echo

#DETERMINISTIC TRACTOGRAPHY
echo "streamlining"
# track \
# -bedpostxdir '../' subject + '.bedpostX/' \
# -inputmodel bedpostx_dyad \
# -seedfile wmparc_invert_bin.nii.gz \
# -curvethresh 90 \
# -curveinterval 2.5 \
# -anisthresh 0.2 \
# -tracker rk4 \
# -interpolator linear \
# -iterations 100 \
# -brainmask nodif_brain_mask.nii.gz | procstreamlines \
# -endpointfile atlas.nii.gz \
# -outputfile bedDetTracts.Bfloat

if probabilistic, run probabilistic.sh