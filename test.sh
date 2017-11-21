#!/bin/bash
PATH=${PATH}:${PWD}

dti_conn.py -v \
    -b dti_conn_test/bedpostx \
    dti_conn_test/T1w_brain.nii.gz \
    dti_conn_test/wmparc_ero.nii.gz \
    dti_conn_test/shen_2mm_268_parcellation.nii.gz \
    dti_conn_test/MNI152_T1_2mm_brain.nii.gz \
    dti_conn_test/subject_b0_bet.nii.gz \
    dti_conn_test/subject_b0_bet_mask.nii.gz \
    dti_conn_test/subject.nii.gz \
    dti_conn_test/dti_conn_test
