#!/usr/bin/env python2.7
import os

def main(path = './', logsPath = None):
	if not path[-1] == '/': path += '/'
	if logsPath is None: 
		logsPath = path
	else: 
		if not logsPath[-1] == '/': logsPath += '/'
		
	fileList = [f for f in os.listdir(logsPath) if '_Master.dat' in f]
	incompleteRuns = []
	for file in fileList:
		with open(logsPath + file) as f: lines = f.readlines()
		if 'Program terminates normally' in lines[-1]: continue
		else: incompleteRuns.append(file)
		
	#print 'Runs:', incompleteRuns
	
	runNums = []
	for file in incompleteRuns:
		f = file.replace('_Master.dat','')
		idx = f[::-1].find('_')
		n = int(f[-idx::])
		runNums.append(n)
	
	if len(runNums) == 0:
		print 'All runs are complete. Nothing to do here.'
		return
		
	runNums.sort()
	runList = ','.join([str(item) for item in runNums])
	#print 'Nums:', runNums
	print 'Incomplete runs:', runList
	
	if os.path.isfile(path + 'mafot.sbatch'):
		with open(path + 'mafot.sbatch') as f: lines = f.readlines()
		with open(path + 'mafotRerun.sbatch','w') as f:
			for line in lines:
				if '#SBATCH --array' in line: f.write('#SBATCH --array=' + runList + '\n')
				else: f.write(line)
	else: print 'No job Queue file found'
	
	return
	
	
if __name__ == '__main__':
	import argparse
	import textwrap
	parser = argparse.ArgumentParser(description = 'Find and restart incomplete MAFOT runs', 
                formatter_class = argparse.RawDescriptionHelpFormatter,
                epilog = textwrap.dedent('''\
                Examples: restartMafot.py'''))
	
	parser.add_argument('-d', '--dir', help = "set working dir", type = str, default = './')
	parser.add_argument('-l', '--logs', help = "set log dir, if different from working dir", type = str, default = None)
	args = parser.parse_args()
	
	main(args.dir, args.logs)
