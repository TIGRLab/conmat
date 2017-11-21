#!/bin/bash
#
# an example method for producing a white matter mask for tractography seeding
# aparc+aseg.nii.gz comes from freesurfer
#

if [ ! -f wmparc_ero.nii.gz ]; then
    3dcalc \
        -a aparc+aseg.nii.gz \
        -expr "equals(a,2)  + \
               equals(a,7)  + \
               equals(a,41) + \
               equals(a,46) + \
               equals(a,251)+ \
               equals(a,252)+ \
               equals(a,253)+ \
               equals(a,254)+ \
               equals(a,255)" \
        -prefix wmparc.nii.gz

    3dcalc \
        -a wmparc.nii.gz \
        -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
        -expr 'a*(1-amongst(0,b,c,d,e,f,g))' \
        -prefix wmparc_ero.nii.gz
fi

