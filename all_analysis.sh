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
fsl2scheme -bvecfile ${bvecfile} -bvalfile ${bvalfile} > "${subject}_${tag}_bvec.scheme"

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


#################################################################
## deterministic tractographt camino
# fitting tensors
wdtfit ${dwifile} "${subject}_${tag}_bvec.scheme" -bgmask ${mask} -outputfile ${subject}_${tag}_wdfit_determinisitc.nii.gz

# calculate streamlines
track -inputfile ${subject}_${tag}_wdfit_determinisitc.nii.gz -inputmodel dt -seedfile reg_wm_to_b0.nii.gz -curvethresh 90 -curveinterval 2.5 -anisthresh 0.2 -tracker rk4 -interpolator linear -stepsize 0.5 -iterations 100 -brainmask ${mask} | procstreamlines -endpointfile reg_atlas_to_b0.nii.gz -outputfile ${subject}_${tag}_wdfit_determinisitc_tracts.Bfloat

# voxel-wise fa and md
md -inputfile ${subject}_${tag}_wdfit_determinisitc.nii.gz -outputfile ${subject}_${tag}_wdfit_determinisitc_md.nii.gz
fa -inputfile ${subject}_${tag}_wdfit_determinisitc.nii.gz -outputfile ${subject}_${tag}_wdfit_determinisitc_fa.nii.gz

# calculate connectivity matrix
conmat -inputfile ${subject}_${tag}_wdfit_determinisitc_tracts.Bfloat -targetfile reg_atlas_to_b0.nii.gz -tractstat length -outputroot ${subject}_${tag}_wdfit_determinisitc_conn_


#####################################
## probabilistic tractography camino
# convert DWI to camino format
image2voxel -4dimage ${dwifile} -outputfile "${subject}_${tag}_dwi.Bfloat"

# fit tensors
modelfit -inputfile "${subject}_${tag}_dwi.Bfloat" -schemefile "${subject}_${tag}_bvec.scheme" -model ldt -brainmask ${mask} -outputfile ${subject}_${tag}_modelfit_probabilistic.Bdouble

cat ${subject}_${tag}_modelfit_probabilistic.Bdouble | dteig > ${subject}_${tag}_modelfit_probabilistic_dteig.Bdouble

# generate Probability Density Functions (PDFs) in each voxel
snr=$(estimatesnr -inputfile "${subject}_${tag}_dwi.Bfloat" -schemefile "${subject}_${tag}_bvec.scheme" -bgmask ${mask} | grep "SNR mult" | tr "\t" " " | tr -s " " | cut -d " " -f3)
dtlutgen -schemefile "${subject}_${tag}_bvec.scheme" -snr ${snr} > ${subject}_${tag}_dtlut.dat
picopdfs -inputmodel dt -luts ${subject}_${tag}_dtlut.dat < ${subject}_${tag}_modelfit_probabilistic.Bdouble > ${subject}_${tag}_modelfit_probabilistic_pdfs.Bdouble

# generate streamlines and connectivity matrix
track -inputmodel pico -seedfile reg_wm_to_b0.nii.gz -iterations 1000 < ${subject}_${tag}_modelfit_probabilistic_pdfs.Bdouble > ${subject}_${tag}_modelfit_probabilistic_tracts.Bfloat
conmat -inputfile picoTracts.Bfloat -targetfile atlas.nii.gz -tractstat length -outputroot conmat_prob_

# voxel-wise fa + md
for metric in fa md; do
    cat ${subject}_${tag}_modelfit_probabilistic.Bdouble | ${metric} | voxel2image -outputroot ${subject}_${tag}_modelfit_probabilistic_${metric} -header ${dwifile}
done

# run FSL's bedpostx thing
mkdir bedpostx
cp ${dwifile} bedpostx/data.nii.gz
cp ${mask} bedpostx/nodif_brain_mask.nii.gz
cp ${bvecfile} bedpostx/bvecs
cp ${bvalfile} bedpostx/bvals
cd bedpostx
bedpostx ./

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
