#!/usr/bin/env python

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


def connectivity():
    os.chdir(tempSubjDir)
    if not os.path.exists(conmatPipeDir):
        os.makedirs(conmatPipeDir)

    if os.listdir(conmatPipeDir)==[]:
        register()
        os.system('echo "streamlining"')
        os.system('track \
            -bedpostxdir ' + tempSubjName + '.bedpostX' + ' \
            -inputmodel bedpostx_dyad \
            -seedfile wmparc_invert_bin.nii.gz \
            -curvethresh 90 \
            -curveinterval 2.5 \
            -anisthresh 0.2 \
            -tracker rk4 \
            -interpolator linear \
            -iterations 1 \
            -brainmask nodif_brain_mask.nii.gz | procstreamlines \
            -endpointfile atlas.nii.gz \
            -outputfile bedDetTracts.Bfloat')

        os.system('echo "calculating connectivity matrix"')
        os.system('conmat \
            -inputfile bedDetTracts.Bfloat \
            -targetfile atlas.nii.gz \
            -tractstat length \
            -outputroot ' + subjectID + '_bedpostX_det_')

        conExt = '*_bedpostX_det_sc.csv'
        lenExt = '*_bedpostX_det_ts.csv'

        shutil.copyfile(tempSubjDir + glob('bedDetTracts.Bfloat')[0], conmatPipeDir + subjectID + '_detTracts.Bfloat')
        shutil.copyfile(tempSubjDir + glob('*.scheme')[0], conmatPipeDir + subjectID + '.scheme')
        shutil.copyfile(tempSubjDir + glob('atlas.nii.gz')[0], conmatPipeDir + subjectID + '_registered_shen.nii.gz')
        shutil.copyfile(tempSubjDir + glob(conExt)[0], conmatPipeDir + subjectID + '_bedpostX_det_connectivity.csv')
        shutil.copyfile(tempSubjDir + glob(lenExt)[0], conmatPipeDir + subjectID + '_bedpostX_det_length.csv')

def get_fa_matrix(stat):
    os.chdir(tempSubjDir)
    faCSV = conmatPipeDir + subjectID + '_bedpostX_det_fa_' + stat + '.csv'
    if not os.path.exists(faCSV):
        try:
            faFile = tempSubjDir + glob('*_FA.nii.gz')[0]

            logger.info('calculating connectivity matrix')
            os.system('conmat \
                -inputfile bedDetTracts.Bfloat \
                -targetfile atlas.nii.gz \
                -scalarfile {} -tractstat {} -outputroot {}_fa_{}_'.format(
                    faFile, stat, subjectID, stat))
        except:
            register()
            logger.info('calculating FA map')
            os.system('fa -inputfile wdt.nii.gz -outputfile fa.nii.gz')

            logger.info('calculating connectivity matrix')
            os.system('conmat \
                -inputfile bedDetTracts.Bfloat \
                -targetfile atlas.nii.gz \
                -scalarfile fa.nii.gz \
                -tractstat {} -outputroot {}_fa_{}_'.format(stat, subjectID, stat))

    faExt = '*_fa_' + stat + '_ts.csv'
    shutil.copyfile(tempSubjDir + glob(faExt)[0], faCSV)


def get_md_matrix(stat):
    os.chdir(tempSubjDir)
    mdCSV = conmatPipeDir + subjectID + '_bedpostX_det_md_' + stat + '.csv'
    mdFile = tempSubjDir + glob('*_MD.nii.gz')[0]
    mdExt = '*_md_' + stat + '_ts.csv'

    if not os.path.exists(mdCSV):
        try:
            logger.info('calculating connectivity matrix')
            os.system('conmat \
                -inputfile bedDetTracts.Bfloat \
                -targetfile atlas.nii.gz \
                -scalarfile {} -tractstat {} -outputroot {}_md_{}_ '.format(
                    mdFile, stat, subjectID, stat))

            shutil.copyfile(tempSubjDir + glob(mdExt)[0], conmatPipeDir + subjectID + '_bedpostX_det_md_' + stat + '.csv')
        except:
            register()
            os.system('md -inputfile wdt.nii.gz -outputfile md.nii.gz')

            os.system('echo "calculating connectivity matrix"')
            os.system('conmat \
                -inputfile bedDetTracts.Bfloat \
                -targetfile atlas.nii.gz \
                -scalarfile md.nii.gz \
                -tractstat {} -outputroot {}_md_{}_'.format(stat, subjectID, stat))
            shutil.copyfile(tempSubjDir + glob(mdExt)[0], conmatPipeDir + subjectID + '_bedpostX_det_md_' + stat + '.csv')


def camino_probabilisitc(dwi, bvec, bval, b0_mask, output_prefix):
    """calculated connectivity matrix using single-tensor dti"""
    run('fsl2scheme -bvecfile {bvec} -bvalfile {bval} > {op}_bvec.scheme'.format(
        bvec=bvec, bval=bval, op=output_prefix))
    run('modelfit -inputfile {op}_dwi.Bfloat -schemefile {op}_bvec.scheme -model ldt -brainmask {b0_mask} -outputfile {op}_modelfit_probabilistic.Bdouble'.format(
        b0_mask=b0_mask, op=output_prefix))
    run('cat {op}_modelfit_probabilistic.Bdouble | dteig > {op}_modelfit_probabilistic_dteig.Bdouble'.format(
        op=output_prefix))
    run('dtlutgen -schemefile {op}_bvec.scheme -snr $(estimatesnr -inputfile {op}_dwi.Bfloat -schemefile {op}_bvec.scheme -bgmask {b0_mask} | grep "SNR mult" | tr "\t" " " | tr -s " " | cut -d " " -f3) > {op}_dtlut.dat'.format(
        b0_mask=b0_mask, op=output_prefix))
    run('picopdfs -inputmodel dt -luts {op}_dtlut.dat < {op}_modelfit_probabilistic.Bdouble > {op}_modelfit_probabilistic_pdfs.Bdouble'.format(
        op=output_prefix))
    run('track -inputmodel pico -seedfile {op}_reg_wm_to_b0.nii.gz -iterations 1000 < {op}_modelfit_probabilistic_pdfs.Bdouble > {op}_modelfit_probabilistic_tracts.Bfloat'.format(
        op=output_prefix))
    run('conmat -inputfile {op}_modelfit_probabilistic_tracts.Bfloat -targetfile {op}_reg_atlas_to_b0.nii.gz -tractstat length -outputroot {op}_conmat_prob_'.format(
        op=output_prefix))
    run('cat {op}_modelfit_probabilistic.Bdouble | fa | voxel2image -outputroot {op}_modelfit_probabilistic_fa -header {dwi}'.format(
        dwi=dwi, op=output_prefix))
    run('cat {op}_modelfit_probabilistic.Bdouble | md | voxel2image -outputroot {op}_modelfit_probabilistic_md -header {dwi}'.format(
        dwi=dwi, op=output_prefix))


def camino_deterministic(dwi, bvec, bval, b0_mask, output_prefix):
    """calculates connectivity matrix using single-tensor dti"""
    run('fsl2scheme -bvecfile {bvec} -bvalfile {bval} > {op}_bvec.scheme'.format(
        bvec=bvec, bval=bval, op=output_prefix))
    run('wdtfit {dwi} {op}_bvec.scheme -bgmask {b0_mask} -outputfile {op}_wdfit_determinisitc.nii.gz'.format(
        dwi=dwi, b0_mask=b0_mask, op=output_prefix))
    run('track -inputfile {op}_wdfit_determinisitc.nii.gz -inputmodel dt -seedfile {op}_reg_wm_to_b0.nii.gz -curvethresh 90 -curveinterval 2.5 -anisthresh 0.2 -tracker rk4 -interpolator linear -stepsize 0.5 -iterations 100 -brainmask {b0_mask} | procstreamlines -endpointfile {op}_reg_atlas_to_b0.nii.gz -outputfile {op}_wdfit_determinisitc_tracts.Bfloat'.format(
        b0_mask=b0_mask, op=output_prefix))
    run('md -inputfile {op}_wdfit_determinisitc.nii.gz -outputfile {op}_wdfit_determinisitc_md.nii.gz'.format(
        op=output_prefix))
    run('fa -inputfile {op}_wdfit_determinisitc.nii.gz -outputfile {op}_wdfit_determinisitc_fa.nii.gz'.format(
        op=output_prefix))
    run('conmat -inputfile {op}_wdfit_determinisitc_tracts.Bfloat -targetfile {op}_reg_atlas_to_b0.nii.gz -tractstat length -outputroot {op}_wdfit_determinisitc_conn_'.format(
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
    if not os.path.isfile('{op}_reg_atlas_to_b0.nii.gz'.format(op=output_prefix):
        run('flirt -in {t1} -ref {mni} -omat {op}_reg_t1_to_mni.mat -o {op}_reg_t1_to_mni.nii.gz'.format(
            t1=t1, mni=atlas_anat, op=output_prefix))
        run('convert_xfm -omat {op}_reg_mni_to_t1.mat -inverse {op}_reg_t1_to_mni.mat'.format(
            op=output_prefix))
        run('convert_xfm -omat {op}_reg_mni_to_b0.mat -concat {op}_reg_t1_to_b0.mat {op}_reg_mni_to_t1.mat'.format(
            op=output_prefix))
        run('flirt -in {atlas} -ref {b0} -interp nearestneighbour -applyxfm -init {op}_reg_mni_to_b0.mat -o {op}_reg_atlas_to_b0.nii.gz'.format(
            atlas=atlas, b0=b0, op=output_prefix))


def run(cmd, dryrun=False, specialquote=True, verbose=True):
    """runs command in default shell, returning STDOUT and a return code."""
    # perform shell quoting for special characters in filenames
    if specialquote:
        cmd = _escape_shell_chars(cmd)

    logger.debug("Executing command: {}".format(cmd))

    p = proc.Popen(cmd, shell=True, stdout=proc.PIPE, stderr=proc.PIPE)
    out, err = p.communicate()

    if p.returncode and verbose:
        logger.error('command {} failed with code {}.\n STDERR: {}'.format(
            cmd, p.returncode, err))

    return p.returncode, out


def check_inputs(args, arglist, directory=False):
    """
    checks for the existance of files submitted in the arglist. if directory is
    true, then checks for directories instead
    """
    for arg in arglist:
        if not os.path.isfile(args.__dict__[arg]):
            logger.error('input {}={} not found'.format(arg, args.__dict__[arg]))
            sys.exit(1)


if __name__ == '__main__':

    argparser = argparse.ArgumentParser()

    # collect input files
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

    import IPython; IPython.embed()

    register(args.b0, args.t1, args.t1_wm, args.atlas, args.atlas_anat, args.output_prefix)
    camino_deterministic(args.dwi, bvec, bval, args.b0_mask, args.output_prefix)
    camino_probabilistic(args.dwi, bvec, bval, args.b0_mask, args.output_prefix)

    # put these inside each camino fxn?
    get_fa_matrix('max')
    get_fa_matrix('mean')
    get_fa_matrix('min')
    get_md_matrix('max')
    get_md_matrix('mean')
    get_md_matrix('min')


