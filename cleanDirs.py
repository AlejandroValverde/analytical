import os
import sys
import getopt

def recursiveFunction():
	cwd = os.getcwd() #Get working directory
	print('-> Exploring: '+cwd)
	for file in os.listdir(cwd):

		if os.path.isdir(file) and not 'git' in file:
			os.chdir(cwd + '\\' + file)
			recursiveFunction()
			os.chdir(cwd)

		elif file.startswith('Thumbs.db'):
			print('--> Thumbs removed in: '+cwd)
			os.remove(file)

		elif file.endswith('.eps') and not 'eps' in cwd:

			print('--> Converting: '+file+ ', in: '+cwd)
			os.system('epstopdf '+file)
			if not os.path.isdir('eps'):
				os.mkdir('.\\eps')
				print('--> Making new eps folder in: '+cwd)
			print('--> Moving file: '+file)
			os.rename(cwd+'\\'+file, cwd+'\\eps\\'+file)


recursiveFunction()

print('---> Execution finished')