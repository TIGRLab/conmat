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

parameter = raw_input('Parameter [fmri, connectivity, fa_stat, md_stat, length]: ')
status = raw_input('Status [control_arm_1, case_arm_2]: ')

fileNames = ''
groupedConmatDir = '/scratch/lliu/' + projectName + '/pipelines/grouped/' + parameter + '/'
demographicsFile = '/archive/data-2.0/SPINS/data/clinical/demographics.csv'

#------------------------------------------------------------------------------------------------
# for subjectID in projectSubjects:
# 	if "PHA" not in subjectID and "CMH" in subjectID:
# 		subjConmatDir = '/scratch/lliu/' + projectName + '/pipelines/conmat/' + subjectID +'/'
# 		if os.path.exists(subjConmatDir):
# 			os.chdir(subjConmatDir)
# 			subjFile = glob.glob('*' + parameter + '.csv')
# 			if subjFile != []:
# 				shutil.copy(subjConmatDir + subjFile[0], groupedConmatDir)
# 				fileNames = fileNames + ' \'' + subjFile[0] + '\''

# os.chdir(groupedConmatDir)
# file = open(parameter + "_fileNames.txt", "w")
# file.write(fileNames)
# file.close()

#------------------------------------------------------------------------------------------------
csv_file = csv.reader(open(demographicsFile, "rb"), delimiter = ",")
for row in csv_file:
	if row[1] == status:
		os.chdir(projectPath)
		extract = glob.glob(row[0] + '*')
		if extract != []:
			subjectID = extract[0]
			# if "PHA" not in subjectID and "CMH" in subjectID and "test" not in subjectID:
			if "PHA" not in subjectID and "test" not in subjectID and "9999" not in subjectID:
				subjConmatDir = '/scratch/lliu/' + projectName + '/pipelines/conmat/' + subjectID +'/'
				if parameter == 'fmri' and os.path.exists(subjConmatDir):
					projectPath = '/archive/data-2.0/' + projectName + '/pipelines/fmri/rest/'
					subjConmatDir = projectPath + subjectID + '/'
					if os.path.exists(subjConmatDir):
						os.chdir(subjConmatDir)
						subjFile = glob.glob('*nonlin_roi-corrs.csv')
						if subjFile != []:
							shutil.copy(subjConmatDir + subjFile[0], groupedConmatDir)
							fileNames = fileNames + ' \'' + subjFile[0] + '\''
				else:
					subjConmatDir = '/scratch/lliu/' + projectName + '/pipelines/conmat/' + subjectID +'/'
					if os.path.exists(subjConmatDir):
						os.chdir(subjConmatDir)
						subjFile = glob.glob('*' + parameter + '.csv')
						print subjFile
						if subjFile != []:
							shutil.copy(subjConmatDir + subjFile[0], groupedConmatDir)
							fileNames = fileNames + ' \'' + subjFile[0] + '\''


os.chdir(groupedConmatDir)
file = open(parameter + "_" + status + "_fileNames.txt", "w")
file.write(fileNames)
file.close()
