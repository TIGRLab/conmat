#!/bin/bash

##SETUP
. $MODULESHOME/init/bash
module load camino
module load FSL
module load matlab
module load Paraview

mkdir pico

##INPUTS
ptID=SPN01_CMH_P001_02_01_DTI60-1000_08_Ax-DTI-60plus5-20iso_eddy_correct

bvecFile=${ptID}.bvec
bvalFile=${ptID}.bval
multiVol=${ptID}.nii.gz
b0Mask=${ptID}_b0_bet_mask.nii.gz
b0Brain=${ptID}_b0_bet.nii.gz
T1wBrain=T1w_brain.nii.gz
wmParc=wmparc.nii.gz
mniBrain=MNI152_T1_2mm_brain.nii.gz
shenParc=shen_2mm_268_parcellation.nii.gz

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

echo "converting DWI to Camino format"
image2voxel \
-4dimage ${multiVol} \
-outputfile dwi.Bfloat
echo

#***************************************************************************
#PROBABILISTIC TRACTOGRAPHY
echo "fitting tensors"
modelfit \
-inputfile dwi.Bfloat \
-schemefile bVectorScheme.scheme \
-model ldt \
-brainmask ${b0Mask} \
-outputfile dt.Bdouble
echo

for PROG in fa md 
do cat dt.Bdouble | ${PROG} | voxel2image -outputroot ${PROG} -header ${multiVol}
done
cat dt.Bdouble | dteig > dteig.Bdouble
# pdview -inputfile dteig.Bdouble -scalarfile fa.nii.gz

echo "calculating snr"
signalNoiseRatio=$(estimatesnr \
	-inputfile dwi.Bfloat \
	-schemefile bVectorScheme.scheme \
	-bgmask ${b0Mask} | grep "SNR mult" | tr "\t" " " | tr -s " " | cut -d " " -f3)
echo

echo "generating look-up table: snr=${signalNoiseRatio}"
dtlutgen \
-schemefile bVectorScheme.scheme \
-snr ${signalNoiseRatio} > picoTable.dat
echo

echo "generating Probability Density Functions (PDFs) in each voxel"
picopdfs \
-inputmodel dt \
-luts picoTable.dat < dt.Bdouble > pdfs.Bdouble
echo

echo "streamlining"
track \
-inputmodel pico \
-seedfile wmparc_invert_bin.nii.gz \
-iterations 100 < pdfs.Bdouble > picoTracts.Bfloat
echo

# ## BedpostX
# track \
# -bedpostxdir ../SPN01_CMH_P001_02.bedpostX/ \
# -inputmodel bedpostx \
# -seedfile wmparc_invert_bin.nii.gz \
# -curvethresh 90 \
# -curveinterval 2.5 \
# -anisthresh 0.2 \
# -tracker rk4 \
# -interpolator linear \
# -iterations 100 \
# -brainmask nodif_brain_mask.nii.gz | procstreamlines \
# -endpointfile atlas.nii.gz \
# -outputfile bedProbTracts.Bfloat

#***************************************************************************
#CONNECTIVITY MATRIX
md \
-inputfile wdt.nii.gz \
-outputfile md.nii.gz

fa \
-inputfile wdt.nii.gz \
-outputfile fa.nii.gz

echo "calculating connectivity matrix"
conmat \
-inputfile picoTracts.Bfloat \
-targetfile atlas.nii.gz \
-tractstat length \
-outputroot conmat_prob_
echo

# get matrices of fa and md

#***************************************************************************
#visuals
# for visualizing picoTracts and ROI probability

# track -inputmodel pico -seedfile wmparc_invert_bin.nii.gz < pdfs.Bdouble | procstreamlines -seedfile wmparc_invert_bin.nii.gz -outputacm -outputroot pico/
# track -inputmodel pico -seedfile wmparc_invert_bin.nii.gz < pdfs.Bdouble | procstreamlines -seedfile wmparc_invert_bin.nii.gz -outputcp -outputroot pico/
# imagestats -stat max -outputroot roiprobs -images pico/*

# fslview pico/acm_sc.nii.gz
# fslview roiprobs.nii.gz

# matlab
# conmat_path = 'conmat_prob_ts.csv';
# myconmat = csvread(conmat_path, 1, 0);
# figure, imagesc(myconmat)
