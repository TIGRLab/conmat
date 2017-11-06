#!/bin/bash

mkdir pico

subject="SPN01_CMH_0001_01_01"
tag="DTI60-1000"

# raw data
bvecfile=$(ls ${subject}_${tag}_*.bvec)
bvalfile=$(ls ${subject}_${tag}_*.bval)
dwifile="SPN01_CMH_0001_01_01_DTI60-1000_20_Ax-DTI-60plus5.nii.gz" # need to be clever here

# from dtifit
mask=$(ls *_b0_bet_mask.nii.gz)
b0=$(ls *_b0_bet.nii.gz)

# from HCP pipeline
wm=wmparc.nii.gz
t1=T1w_brain.nii.gz

# included in datman
mni=MNI152_T1_2mm_brain.nii.gz
atlas=shen_2mm_268_parcellation.nii.gz

## NB: runs after DTIFIT, requires some of the outputs

## preprocessing so all analysis can be done in single subject space
# make scheme file
fsl2scheme -bvecfile ${bvecFile} -bvalfile ${bvalFile} > "${subject}_${tag}_bvec.scheme"

# register b0 to T1 and vice-versa
flirt -in ${b0} -ref ${t1} -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -cost mutualinfo -omat reg_b0_to_t1.mat -o reg_b0_to_t1.nii.gz
convert_xfm -omat reg_t1_to_b0.mat -inverse reg_b0_to_t1.mat

# register wm mask to b0
flirt -in ${wm} -ref ${b0} -interp nearestneighbour -applyxfm -init reg_t1_to_b0.mat -o reg_wm_to_b0_tmp.nii.gz
fslmaths reg_wm_to_b0_tmp.nii.gz -thr 2500 -bin reg_wm_to_b0.nii.gz
rm reg_wm_to_b0_tmp.nii.gz

# register atlas to b0
flirt -in ${t1} -ref ${mni} -omat reg_t1_to_mni.mat -o reg_t1_to_mni.nii.gz
convert_xfm -omat reg_mni_to_t1.mat -inverse reg_t1_to_mni.mat
convert_xfm -omat reg_mni_to_b0.mat -concat reg_t1_to_b0.mat reg_mni_to_t1.mat
flirt -in ${atlas} -ref ${b0} -interp nearestneighbour -applyxfm -init reg_mni_to_b0.mat -o reg_atlas_to_b0.nii.gz


## deterministic tractographt camino
# fitting tensors
wdtfit ${multiVol} bVectorScheme.scheme -brainmask ${b0Mask} -outputfile wdt.nii.gz
# calculate streamlines
track -inputfile wdt.nii.gz -inputmodel dt -seedfile wmparc_invert_bin.nii.gz \
    -curvethresh 90 -curveinterval 2.5 -anisthresh 0.2 -tracker rk4
    -interpolator linear -stepsize 0.5 -iterations 100 -brainmask ${b0Mask} | procstreamlines \
    -endpointfile atlas.nii.gz -outputfile camDetTracts.Bfloat

# voxel-wise fa and md
md -inputfile wdt.nii.gz -outputfile md.nii.gz
fa -inputfile wdt.nii.gz -outputfile fa.nii.gz

# calculate connectivity matrix
conmat -inputfile camDetTracts.Bfloat -targetfile atlas.nii.gz -tractstat length -outputroot conmat_det_


## probabilistic tractography camino
# convert DWI to camino format
image2voxel -4dimage ${multiVol} -outputfile dwi.Bfloat

# fit tensors
modelfit -inputfile dwi.Bfloat -schemefile bVectorScheme.scheme -model ldt -brainmask ${b0Mask} -outputfile dt.Bdouble

for metric in fa md; do
    cat dt.Bdouble | ${metric} | voxel2image -outputroot ${PROG} -header ${multiVol}
done

cat dt.Bdouble | dteig > dteig.Bdouble

# calculate snr
signalNoiseRatio=$(estimatesnr -inputfile dwi.Bfloat -schemefile bVectorScheme.scheme -bgmask ${b0Mask} | grep "SNR mult" | tr "\t" " " | tr -s " " | cut -d " " -f3)
echo "snr=${signalNoiseRatio}"

# generate look up table (lut)
dtlutgen -schemefile bVectorScheme.scheme -snr ${signalNoiseRatio} > picoTable.dat

# generate Probability Density Functions (PDFs) in each voxel
picopdfs -inputmodel dt -luts picoTable.dat < dt.Bdouble > pdfs.Bdouble

# generate streamlines
track -inputmodel pico -seedfile wmparc_invert_bin.nii.gz -iterations 100 < pdfs.Bdouble > picoTracts.Bfloat

# voxel-wise fa + md
md -inputfile wdt.nii.gz -outputfile md.nii.gz
fa -inputfile wdt.nii.gz -outputfile fa.nii.gz

# calculating connectivity matrix
conmat -inputfile picoTracts.Bfloat -targetfile atlas.nii.gz -tractstat length -outputroot conmat_prob_




## BedpostX probabilistic
track -bedpostxdir ../${subject}.bedpostX/ -inputmodel bedpostx -seedfile wmparc_invert_bin.nii.gz \
    -curvethresh 90 -curveinterval 2.5 -anisthresh 0.2 -tracker rk4 -interpolator linear \
    -iterations 100 -brainmask nodif_brain_mask.nii.gz | procstreamlines \
    -endpointfile atlas.nii.gz -outputfile bedProbTracts.Bfloat

## BedpostX deterministic
track -bedpostxdir ../SPN01_CMH_P001_02.bedpostX/ -inputmodel bedpostx_dyad \
    -seedfile wmparc_invert_bin.nii.gz -curvethresh 90 -curveinterval 2.5 \
    -anisthresh 0.2 -tracker rk4 -interpolator linear -iterations 100 \
    -brainmask nodif_brain_mask.nii.gz | procstreamlines \
    -endpointfile atlas.nii.gz -outputfile bedDetTracts.Bfloat



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
