#!/usr/bin/env python

import os, sys
from glob import glob
import shutil
import sys
import argparse
import logging

logging.basicConfig(level=logging.WARN, format="[%(name)s] %(levelname)s: %(message)s")
logger = logging.getLogger(os.path.basename(__file__))
#tmpdir = tempfile.mkdtemp(prefix='qa-dti-')


def remove_temp_files():
    """unused"""
	shutil.rmtree(tempSubjDir)
	shutil.rmtree(tempSubjName + '.bedpostX')


def copy_bedpost_to_temp():
	shutil.copytree(bedpostXPipeDir + subjectID + '.bedpostX', tempSubjName + '.bedpostX')


def copy_bedpost_from_temp():
	shutil.copytree(tempSubjName + '.bedpostX', bedpostXPipeDir + subjectID + '.bedpostX')


def run_bedpost():
	os.chdir(tempSubjDir)
	os.system('bedpostx ./')


def get_files():
    """finds all input files and copies them to a working directory"""
	os.chdir(dtiPipeDir)
	if glob('*_eddy_correct.nii.gz') != []:
		shutil.copytree(dtiPipeDir, tempSubjDir)
		shutil.copy(hcpPipeDir + 'wmparc.nii.gz', tempSubjDir)
		shutil.copy(hcpPipeDir + 'T1w_brain.nii.gz', tempSubjDir)
		shutil.copy(mniAtlasDir + 'shen_2mm_268_parcellation.nii.gz', tempSubjDir)
		shutil.copy(mniAtlasDir + 'MNI152_T1_2mm_brain.nii.gz', tempSubjDir)

		os.chdir(tempSubjDir)

		# renaming necessary files to proper input names
		shutil.copyfile(tempSubjDir + glob('*.bval')[0], tempSubjDir + 'bvals')
		shutil.copyfile(tempSubjDir + glob('*.bvec')[0], tempSubjDir + 'bvecs')
		shutil.copyfile(tempSubjDir + glob('*_eddy_correct_b0_bet_mask.nii.gz')[0], tempSubjDir + 'nodif_brain_mask.nii.gz')
		shutil.copyfile(tempSubjDir + glob('*_eddy_correct.nii.gz')[0], tempSubjDir + 'data.nii.gz')
	else:
		return


def register():
	"""
    co-registers all necessary B0s from the DWI to the subject's T1, and the
    subject's T1 to the MNI space associated with the atlas used. Outputs all
    transforms, images, and masks in the same voxel space.
    """
	os.chdir(tempSubjDir)

	b0Mask = glob('*_b0_bet_mask.nii.gz')[0]
	b0Brain = glob('*_b0_bet.nii.gz')[0]
	multiVol = glob('*_eddy_correct.nii.gz')[0]
	T1wBrain = 'T1w_brain.nii.gz'
	wmParc = 'wmparc.nii.gz'
	mniBrain = 'MNI152_T1_2mm_brain.nii.gz'
	shenParc = 'shen_2mm_268_parcellation.nii.gz'
	bvecFile = '*.bvec'
	bvalFile = '*.bval'

    logger.info('making scheme file')
	os.system('fsl2scheme -bvecfile {} -bvalfile {} > bVectorScheme.scheme'.format(
        bvecFile, bvalFile))

	logger.info('registering B0 to T1')
	os.system('flirt -in {} -ref {} -omat tfm.mat'.format(b0Brain, T1wBrain))
	os.system('convert_xfm -omat tfm_invert.mat -inverse tfm.mat')

    logger.info('registering white matter mask')
	os.system('flirt -in {} -ref {} -applyxfm -init tfm_invert.mat -o wmparc_invert.nii.gz'.format(
        wmParc, b0Brain))
	os.system('fslmaths wmparc_invert.nii.gz -thr 2500 -bin wmparc_invert_bin.nii.gz')

	logger.info('registering MNI')
	os.system('flirt -in {} -ref {} -interp nearestneighbour -omat mni.mat'.format(
        mniBrain, b0Brain))

	logger.info('registering atlas')
	os.system('flirt -in {} -ref {} -interp nearestneighbour -applyxfm -init mni.mat -o atlas.nii'.format(
        shenParc, b0Brain))

	logger.info('fitting tensors')
	os.system('wdtfit {} bVectorScheme.scheme -brainmask {} -outputfile wdt.nii.gz'.format(
        multiVol, b0Mask))


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



def main():
    """
    1) runs bedpostX
    2) runs streamline tractography
    3) creates connectivity matrices
    """
	print(subjectID)
	if os.path.exists(dtiPipeDir):
		if os.path.exists(tempSubjDir):
			flag = 0
			while flag != 1:
				if os.path.exists(tempSubjName + '.bedpostX'):
					if os.path.exists(bedpostXPipeDir):
						flag = 1
					elif os.path.exists(eyeFile):
						copyBedpostX()
						flag = 1
					else:
						flag = 0
				else:
					if os.path.exists(bedpostXPipeDir):
						copy_bedpost_to_temp()
						flag = 1
					else:
						run_bedpost()
						flag = 0
			flag = 0
			while flag != 1:
				connectivity()
				get_fa_matrix('max')
				get_fa_matrix('mean')
				get_fa_matrix('min')
				get_md_matrix('max')
				get_md_matrix('mean')
				get_md_matrix('min')
				flag = 1
		else:
			copy_bedpost_from_temp()


if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    argparser.add_argument('t1', help='T1 weighted image (deskulled)')
    argparser.add_argument('t1_wm', help='White matter mask in register wiht T1 (voxel dimensions must match)')
    argparser.add_argument('atlas' help='atlas file with integer ROIs')
    argparser.add_argument('atlas_anat' help='T1 weighted image associated with atlas')



    argparser.add_argument('bval', help='bval file from spherical phantom')
    argparser.add_argument('output_prefix', help='full path to output prefix')
    argparser.add_argument('-a', '--accel', help='y/n, y=nyquist accelerated data [default=n]')
    argparser.add_argument('-v', '--verbose', action="count", help='turns on debug messages')
    args = argparser.parse_args()

    projectName = sys.argv[1]
    subjectID = sys.argv[2]

    # input directories
    dtiPipeDir = '/archive/data-2.0/' + projectName + '/pipelines/dtifit/' + subjectID + '/'
    hcpPipeDir = '/archive/data-2.0/' + projectName + '/pipelines/hcp/' + subjectID + '/T1w/'
    mniAtlasDir = '/projects/lliu/ImportantFiles/'
    eyeFile = tempSubjName + '.bedpostX/xfms/eye.mat'

    # output directories
    bedpostXPipeDir = '/scratch/lliu/' + projectName + '/pipelines/bedpostX/' + subjectID + '/'
    conmatPipeDir = '/scratch/lliu/' + projectName + '/pipelines/conmat/' + subjectID + '/'
    tempSubjDir = '/scratch/lliu/tmp/' + subjectID + '/'
    tempSubjName = '/scratch/lliu/tmp/' + subjectID

	main()

