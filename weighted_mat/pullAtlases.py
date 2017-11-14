import os
import glob
import shutil
import csv

#this program pulls the connectivity matrices and puts them in a folder, and then makes a text file that has a string of all the file names separated by spaces

# enter project directory
projectName = raw_input('Project Name: ')
projectPath = '/archive/data-2.0/' + projectName + '/data/nii/'

os.chdir(projectPath)

projectSubjects = glob.glob('*')
projectSubjects.sort()

fileNames = ''
# groupedConmatDir = '/scratch/lliu/' + projectName + '/pipelines/grouped/' + parameter + '/'
binningDir = '/home/lliu/Desktop/BINNING/'
demographicsFile = '/archive/data-2.0/SPINS/data/clinical/demographics.csv'

#------------------------------------------------------------------------------------------------
csv_file = csv.reader(open(demographicsFile, "rb"), delimiter = ",")
for row in csv_file:
	os.chdir(projectPath)
	extract = glob.glob(row[0] + '*')
	if extract != []:
		subjectID = extract[0]
		# if "PHA" not in subjectID and "CMH" in subjectID and "test" not in subjectID:
		if "PHA" not in subjectID and "test" not in subjectID and "9999" not in subjectID:
			subjConmatDir = '/scratch/lliu/' + projectName + '/pipelines/conmat/' + subjectID +'/'
			if os.path.exists(subjConmatDir):
				os.chdir(subjConmatDir)
				subjFile = glob.glob('*registered_shen.nii.gz')
				if subjFile != []:
					shutil.copy(subjConmatDir + subjFile[0], binningDir)

