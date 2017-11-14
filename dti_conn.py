#!/usr/bin/env python
"""
calculates structural connectivity between ROIs in the input atlas using
multiple diffusion weighted imaging (DTI) methods implemented in CAMINO, both
probabilistic and deterministic. Can optionally be used to estimate connectivity
using FSL's BEDPOSTX outputs as well.

Main references:
+ J. T. Duda, P. A. Cook, and J. C. Gee, “Reproducibility of graph metrics of
  human brain structural networks.,” Front Neuroinform, vol. 8, p. 46, 2014.
+ Parker GJM et al, A Framework for a Streamline-Based Probabilistic Index of
  Connectivity (PICo) using a Structural Interpretation of MRI Diffusion
  Measurements, Journal of Magnetic Resonance Imaging, 18, 242-254, 2003.
+ P. A. Cook, Y. Bai, S. Nedjati-Gilani, K. K. Seunarine, M. G. Hall, G. J.
  Parker, D. C. Alexander, Camino: Open-Source Diffusion-MRI Reconstruction and
  Processing, 14th Scientific Meeting of the International Society for Magnetic
  Resonance in Medicine, Seattle, WA, USA, p. 2759, May 2006.
+ http://camino.cs.ucl.ac.uk/index.php?n=Main.Citations
"""

import os, sys
from glob import glob
import shutil
import sys
import argparse
import logging
import subprocess as proc

logging.basicConfig(level=logging.WARN, format="[%(name)s] %(levelname)s: %(message)s")
logger = logging.getLogger(os.path.basename(__file__))
#tmpdir = tempfile.mkdtemp(prefix='dti-conn-')

# eventually could be used to a yaml configuration file or similar
OPTIONS = {
    'iterations': 1000,
    'curvethresh': 90,
    'curveinterval': 2.5,
    'anisthresh': 0.2,
    'tracker': 'rk4',
    'interpolator': 'linear',
    'stepsize': 0.5
}


def camino_probabilisitc(dwi, bvec, bval, b0_mask, output_prefix):
    """
    calculated connectivity matrix using single-tensor dti
    Parker GJM et al, A Framework for a Streamline-Based Probabilistic Index of Connectivity (PICo) using a Structural Interpretation of MRI Diffusion Measurements, Journal of Magnetic Resonance Imaging, 18, 242-254, 2003.
    Cook PA et al, Modelling noise-induced fibre-orientation error in diffusion-tensor MRI, IEEE International Symposium on Biomedical Imaging, 332-335, 2004.
    Cook PA, Alexander DC,Modelling uncertainty in two fibre-orientation estimates within a voxel, Proc. Intl. Soc. Mag. Reson. Med. 14 (2006) pp 1629.
    """

    if not os.path.isfile('{op}_modelfit_probabilistic.Bdouble'.format(op=output_prefix)):
        # make scheme file
        run('fsl2scheme -bvecfile {bvec} -bvalfile {bval} > {op}_bvec.scheme'.format(
            bvec=bvec, bval=bval, op=output_prefix))

        # convert DWI to camino format
        run('image2voxel -4dimage {dwi} -outputfile {op}_dwi.Bfloat'.format(
            dwi=dwi, op=output_prefix))

        # fit tensors
        run('modelfit -inputfile {op}_dwi.Bfloat -schemefile {op}_bvec.scheme -model ldt -brainmask {b0_mask} -outputfile {op}_modelfit_probabilistic.Bdouble'.format(
            b0_mask=b0_mask, op=output_prefix))

    if not os.path.isfile('{op}_modelfit_probabilistic_dteig.Bdouble'.format(op=output_prefix)):
        # calculate eigenvectors + eigenvalues from tensors
        run('cat {op}_modelfit_probabilistic.Bdouble | dteig > {op}_modelfit_probabilistic_dteig.Bdouble'.format(
            op=output_prefix))

    if not os.path.isfile('{op}_modelfit_probabilistic_tracts.Bfloat'.format(op=output_prefix)):
        # generate Proability Density Functions (PDFs) in each voxel
        run('dtlutgen -schemefile {op}_bvec.scheme -snr $(estimatesnr -inputfile {op}_dwi.Bfloat -schemefile {op}_bvec.scheme -bgmask {b0_mask} | grep "SNR mult" | tr "\t" " " | tr -s " " | cut -d " " -f3) > {op}_dtlut.dat'.format(
            b0_mask=b0_mask, op=output_prefix))
        run('picopdfs -inputmodel dt -luts {op}_dtlut.dat < {op}_modelfit_probabilistic.Bdouble > {op}_modelfit_probabilistic_pdfs.Bdouble'.format(
            op=output_prefix))

        # generate streamlines
        run('track -inputmodel pico -seedfile {op}_reg_wm_to_b0.nii.gz -iterations {iterations} < {op}_modelfit_probabilistic_pdfs.Bdouble > {op}_modelfit_probabilistic_tracts.Bfloat'.format(
            op=output_prefix, **OPTIONS))

    # generate voxel-wise fa + md metrics
    if not os.path.isfile('{op}_modelfit_probabilistic_md.nii.gz'.format(op=output_prefix)):
        run('cat {op}_modelfit_probabilistic.Bdouble | fa | voxel2image -outputroot {op}_modelfit_probabilistic_fa -header {dwi}'.format(
            dwi=dwi, op=output_prefix))
        run('cat {op}_modelfit_probabilistic.Bdouble | md | voxel2image -outputroot {op}_modelfit_probabilistic_md -header {dwi}'.format(
            dwi=dwi, op=output_prefix))

    # calculate structural connectivity
    if not os.path.isfile('{op}_modelfit_probabilistic_sc.csv'.format(op=output_prefix)):
        run('conmat -inputfile {op}_modelfit_probabilistic_tracts.Bfloat -targetfile {op}_reg_atlas_to_b0.nii.gz -tractstat length -outputroot {op}_modelfit_probabilistic_'.format(
            op=output_prefix))


def camino_deterministic(dwi, bvec, bval, b0_mask, output_prefix):
    """
    calculates connectivity matrix using single-tensor dti
    J. T. Duda et al, “Reproducibility of graph metrics of human brain structural networks.,” Front Neuroinform, vol. 8, p. 46, 2014.
    Peter J. Basser, et al. In vivo fiber tractography using DT-MRI data. Magnetic Resonance in Medicine, 2000.
    Mariana Lazar, et al.. White matter tractography using diffusion tensor deflection. Human Brain Mapping, 2003.
    """
    # calculate deterministic streamlines
    if not os.path.isfile('{op}_wdfit_determinisitc_tracts.Bfloat'.format(op=output_prefix)):
        run('fsl2scheme -bvecfile {bvec} -bvalfile {bval} > {op}_bvec.scheme'.format(
            bvec=bvec, bval=bval, op=output_prefix))
        run('wdtfit {dwi} {op}_bvec.scheme -bgmask {b0_mask} -outputfile {op}_wdfit_determinisitc.nii.gz'.format(
            dwi=dwi, b0_mask=b0_mask, op=output_prefix))
        run('track -inputfile {op}_wdfit_determinisitc.nii.gz -inputmodel dt -seedfile {op}_reg_wm_to_b0.nii.gz -curvethresh {curvethresh} -curveinterval {curveinterval} -anisthresh {anisthresh} -tracker {tracker} -interpolator {interpolator} -stepsize {stepsize} -iterations {iterations} -brainmask {b0_mask} | procstreamlines -endpointfile {op}_reg_atlas_to_b0.nii.gz -outputfile {op}_wdfit_determinisitc_tracts.Bfloat'.format(
            b0_mask=b0_mask, op=output_prefix, **OPTIONS))

    # calculate MD and FA per voxel (single tensor)
    if not os.path.isfile('{op}_wdfit_determinisitc_md.nii.gz'.format(op=output_prefix)):
        run('md -inputfile {op}_wdfit_determinisitc.nii.gz -outputfile {op}_wdfit_determinisitc_md.nii.gz'.format(
            op=output_prefix))
    if not os.path.isfile('{op}_wdfit_determinisitc_fa.nii.gz'.format(op=output_prefix)):
        run('fa -inputfile {op}_wdfit_determinisitc.nii.gz -outputfile {op}_wdfit_determinisitc_fa.nii.gz'.format(
            op=output_prefix))

    # max, mean, min FA and MD, and fiber count per connection
    if not os.path.isfile('{op}_wdfit_determinisitc_conn_sc.csv'.format(op=output_prefix)):
        run('conmat -inputfile {op}_wdfit_determinisitc_tracts.Bfloat -targetfile {op}_reg_atlas_to_b0.nii.gz -scalarfile {op}_wdfit_determinisitc_fa.nii.gz -tractstat max -outputroot {op}_wdfit_deterministic_fa_max_'.format(
            op=output_prefix))
        run('conmat -inputfile {op}_wdfit_determinisitc_tracts.Bfloat -targetfile {op}_reg_atlas_to_b0.nii.gz -scalarfile {op}_wdfit_determinisitc_fa.nii.gz -tractstat mean -outputroot {op}_wdfit_deterministic_fa_mean_'.format(
            op=output_prefix))
        run('conmat -inputfile {op}_wdfit_determinisitc_tracts.Bfloat -targetfile {op}_reg_atlas_to_b0.nii.gz -scalarfile {op}_wdfit_determinisitc_fa.nii.gz -tractstat min -outputroot {op}_wdfit_deterministic_fa_min_'.format(
            op=output_prefix))
        run('conmat -inputfile {op}_wdfit_determinisitc_tracts.Bfloat -targetfile {op}_reg_atlas_to_b0.nii.gz -scalarfile {op}_wdfit_determinisitc_md.nii.gz -tractstat max -outputroot {op}_wdfit_deterministic_md_max_'.format(
            op=output_prefix))
        run('conmat -inputfile {op}_wdfit_determinisitc_tracts.Bfloat -targetfile {op}_reg_atlas_to_b0.nii.gz -scalarfile {op}_wdfit_determinisitc_md.nii.gz -tractstat mean -outputroot {op}_wdfit_deterministic_md_mean_'.format(
            op=output_prefix))
        run('conmat -inputfile {op}_wdfit_determinisitc_tracts.Bfloat -targetfile {op}_reg_atlas_to_b0.nii.gz -scalarfile {op}_wdfit_determinisitc_md.nii.gz -tractstat min -outputroot {op}_wdfit_deterministic_md_min_'.format(
            op=output_prefix))
        # fiber count
        run('conmat -inputfile {op}_wdfit_determinisitc_tracts.Bfloat -targetfile {op}_reg_atlas_to_b0.nii.gz -tractstat length -outputroot {op}_wdfit_determinisitc_conn_'.format(
            op=output_prefix))


def bedpostx(bpost_dir, b0_mask, output_prefix):
    """computes probabilisitc and deterministic bedpostx outputs"""

    if not os.path.isfile('{op}_bedpost_probabilistic_sc.csv'.format(op=output_prefix)):
        # bedpostx probabilistic tractography / structural connectivity
        run('track -bedpostxdir {bpost_dir} -inputmodel bedpostx -seedfile {op}_reg_wm_to_b0.nii.gz -curvethresh {curvethresh} -curveinterval {curveinterval} -anisthresh {anisthresh} -tracker {tracker} -interpolator {interpolator} -iterations {iterations} -brainmask {b0_mask} | procstreamlines -endpointfile {op}_reg_atlas_to_b0.nii.gz -outputfile {op}_bedpost_probabilistic_tracts.Bfloat'.format(
            bpost_dir=bpost_dir, b0_mask=b0_mask, op=output_prefix, **OPTIONS))
        run('conmat -inputfile {op}_bedpost_probabilistic_tracts.Bfloat -targetfile {op}_reg_atlas_to_b0.nii.gz -tractstat length -outputroot {op}_bedpost_probabilistic_'.format(
            op=output_prefix))

    if not os.path.isfile('{op}_bedpost_deterministic_sc.csv'.format(op=output_prefix)):
        # bedpostx deterministic tractography / structural connectivity
        run('track -bedpostxdir {bpost_dir} -inputmodel bedpostx_dyad -seedfile {op}_reg_wm_to_b0.nii.gz -curvethresh {curvethresh} -curveinterval {curveinterval} -anisthresh {anisthresh} -tracker {tracker} -interpolator {interpolator} -iterations {iterations} -brainmask {b0_mask} | procstreamlines -endpointfile {op}_reg_atlas_to_b0.nii.gz -outputfile {op}_bedpost_deterministic_tracts.Bfloat'.format(
            bpost_dir=bpost_dir, b0_mask=b0_mask, op=output_prefix, **OPTIONS))
        run('conmat -inputfile {op}_bedpost_deterministic_tracts.Bfloat -targetfile {op}_reg_atlas_to_b0.nii.gz -tractstat length -outputroot {op}_bedpost_deterministic_'.format(
            op=output_prefix))


def register(b0, t1, t1_wm, atlas, atlas_anat, output_prefix):
    """moves all data to dwi-native single subject space"""

    # register b0 to T1 and vice-versa
    if not os.path.isfile('{op}_reg_t1_to_b0.mat'.format(op=output_prefix)):
        run('flirt -in {b0} -ref {t1} -cost mutualinfo -omat {op}_reg_b0_to_t1.mat -o {op}_reg_b0_to_t1.nii.gz'.format(
            b0=b0, t1=t1, op=output_prefix))
        run('convert_xfm -omat {op}_reg_t1_to_b0.mat -inverse {op}_reg_b0_to_t1.mat'.format(
            op=output_prefix))

    # register wm mask to b0
    if not os.path.isfile('{op}_reg_wm_to_b0.nii.gz'.format(op=output_prefix)):
        run('flirt -in {wm} -ref {b0} -interp nearestneighbour -applyxfm -init {op}_reg_t1_to_b0.mat -o {op}_reg_wm_to_b0_tmp.nii.gz'.format(
            wm=t1_wm, b0=b0, op=output_prefix))
        run('fslmaths {op}_reg_wm_to_b0_tmp.nii.gz -thr 2500 -bin {op}_reg_wm_to_b0.nii.gz'.format(
            op=output_prefix))
        os.remove('{op}_reg_wm_to_b0_tmp.nii.gz'.format(op=output_prefix))

    # register atlas to b0
    if not os.path.isfile('{op}_reg_atlas_to_b0.nii.gz'.format(op=output_prefix)):
        run('flirt -in {t1} -ref {mni} -omat {op}_reg_t1_to_mni.mat -o {op}_reg_t1_to_mni.nii.gz'.format(
            t1=t1, mni=atlas_anat, op=output_prefix))
        run('convert_xfm -omat {op}_reg_mni_to_t1.mat -inverse {op}_reg_t1_to_mni.mat'.format(
            op=output_prefix))
        run('convert_xfm -omat {op}_reg_mni_to_b0.mat -concat {op}_reg_t1_to_b0.mat {op}_reg_mni_to_t1.mat'.format(
            op=output_prefix))
        run('flirt -in {atlas} -ref {b0} -interp nearestneighbour -applyxfm -init {op}_reg_mni_to_b0.mat -o {op}_reg_atlas_to_b0.nii.gz'.format(
            atlas=atlas, b0=b0, op=output_prefix))


def run(cmd):
    """runs command in default shell, returning STDOUT and a return code."""

    logger.debug("executing command: {}".format(cmd))
    p = proc.Popen(cmd, shell=True, stdout=proc.PIPE, stderr=proc.PIPE)
    out, err = p.communicate()

    if p.returncode:
        logger.error('command {} failed with code {}.\n STDERR: {}'.format(
            cmd, p.returncode, err))
        raise Exception('{} failed'.format(cmd))

    return p.returncode, out


def check_inputs(args, arglist, directory=False):
    """
    checks for the existance of files submitted in the arglist. if directory is
    true, then checks for directories instead
    """
    for arg in arglist:
        if directory:
            if not os.path.isdir(args.__dict__[arg]):
                logger.error('input directory {}={} not found'.format(
                    arg, args.__dict__[arg]))
                sys.exit(1)
        else:
            if not os.path.isfile(args.__dict__[arg]):
                logger.error('input file {}={} not found'.format(
                    arg, args.__dict__[arg]))
                sys.exit(1)


if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    argparser.add_argument('t1', help='T1 weighted image (deskulled)')
    argparser.add_argument('t1_wm', help='White matter mask in register wiht T1 (voxel dimensions must match)')
    argparser.add_argument('atlas', help='atlas file with integer ROIs')
    argparser.add_argument('atlas_anat', help='T1 weighted image associated with atlas')
    argparser.add_argument('b0', help='represenative b0 volume from dwi')
    argparser.add_argument('b0_mask', help='brain mask of b0')
    argparser.add_argument('dwi', help='diffusion weighted imaging data')
    argparser.add_argument('output_prefix', help='full path to output prefix')
    argparser.add_argument('-v', '--verbose', action="count", help='turns on debug messages')
    argparser.add_argument('-b', '--bedpostx', help='optional fsl bedpostx input directory')
    args = argparser.parse_args()

    # set debugging
    logger.info('starting')
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # check input directories
    check_inputs(args, ['t1', 't1_wm', 'atlas', 'atlas_anat', 'b0', 'b0_mask', 'dwi'])
    if args.bedpostx:
        check_inputs(args, ['bedpostx'], directory=True)

    # collect bvec and bvals
    if args.dwi.endswith('nii.gz'):
        dwistem = os.path.splitext(os.path.splitext(args.dwi)[0])[0]
    else:
        dwistem = os.path.splitext(args.dwi)[0]

    bvec = '{}.bvec'.format(dwistem)
    bval = '{}.bval'.format(dwistem)

    try:
        register(args.b0, args.t1, args.t1_wm, args.atlas, args.atlas_anat, args.output_prefix)
    except:
        logger.error('failed to coregister all inputs to subject native space')
        sys.exit(1)

    try:
        camino_deterministic(args.dwi, bvec, bval, args.b0_mask, args.output_prefix)
    except:
        logger.error('failed to calculate deterministic structural connectivity using camino')
        sys.exit(1)

    try:
        camino_probabilisitc(args.dwi, bvec, bval, args.b0_mask, args.output_prefix)
    except:
        logger.error('failed to calculate probabilistic structural connectivity using camino')
        sys.exit(1)

    if args.bedpostx:
        try:
            bedpostx(args.bedpostx, args.b0_mask, args.output_prefix)
        except:
            logger.error('failed you calculate structural connectivity using bedpostx')
            sys.exit(1)

