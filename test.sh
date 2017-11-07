#!/bin/bash
PATH=${PATH}:${PWD}

dti_conn.py \
    ~/dti_test/T1w_brain.nii.gz \
    ~/dti_test/wmparc.nii.gz \
    ~/dti_test/shen_2mm_268_parcellation.nii.gz \
    ~/dti_test/MNI152_T1_2mm_brain.nii.gz \
    ~/dti_test/SPN01_CMH_0001_01_01_DTI60-1000_20_Ax-DTI-60plus5_eddy_correct_b0_bet.nii.gz \
    ~/dti_test/SPN01_CMH_0001_01_01_DTI60-1000_20_Ax-DTI-60plus5_eddy_correct_b0_bet_mask.nii.gz \
    ~/dti_test/SPN01_CMH_0001_01_01_DTI60-1000_20_Ax-DTI-60plus5.nii.gz \
    ~/dti_test/dti_conn_test
