#!/bin/bash
PATH=${PATH}:${PWD}

dti_conn.py -v \
    -b ~/dti_conn/bedpostx \
    ~/dti_conn/T1w_brain.nii.gz \
    ~/dti_conn/wmparc.nii.gz \
    ~/dti_conn/shen_2mm_268_parcellation.nii.gz \
    ~/dti_conn/MNI152_T1_2mm_brain.nii.gz \
    ~/dti_conn/SPN01_CMH_0001_01_01_DTI60-1000_20_Ax-DTI-60plus5_eddy_correct_b0_bet.nii.gz \
    ~/dti_conn/SPN01_CMH_0001_01_01_DTI60-1000_20_Ax-DTI-60plus5_eddy_correct_b0_bet_mask.nii.gz \
    ~/dti_conn/SPN01_CMH_0001_01_01_DTI60-1000_20_Ax-DTI-60plus5.nii.gz \
    ~/dti_conn/dti_conn_test
