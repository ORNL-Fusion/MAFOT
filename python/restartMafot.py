#!/usr/bin/env python2.7
import os

def findIncomplete(path = './', logsPath = None, jobFile = 'mafot.sbatch'):
	if not path[-1] == '/': path += '/'
	if logsPath is None: 
		logsPath = path
	else: 
		if not logsPath[-1] == '/': logsPath += '/'
	
	allFiles = os.listdir(logsPath)
	fileList = [f for f in allFiles if '_Master.dat' in f]
	inputList = [f for f in allFiles if (('_inner_' in f) | ('_outer_' in f) | ('_shelf_' in f))]
	incompleteRuns = []
	for file in fileList:
		with open(logsPath + file) as f: lines = f.readlines()
		if 'Program terminates normally' in lines[-1]: continue
		else: incompleteRuns.append(file)
		
	#print 'Runs:', incompleteRuns
	runNums = findMissing(fileList, inputList)	
	#runNums = []
	for file in incompleteRuns:
		f = file.replace('_Master.dat','')
		idx = f[::-1].find('_')
		try: n = int(f[-idx::])
		except: continue
		runNums.append(n)
	
	if len(runNums) == 0:
		print 'All runs are complete. Nothing to do here.'
		return []
		
	runNums.sort()
	runList = ','.join([str(item) for item in runNums])
	#print 'Nums:', runNums
	print 'List of runs incomplete:', runList
	
	if os.path.isfile(path + jobFile):
		with open(path + jobFile) as f: lines = f.readlines()
		with open(path + 'rerun_' + jobFile,'w') as f:
			for line in lines:
				if '#SBATCH --array' in line: f.write('#SBATCH --array=' + runList + '\n')
				else: f.write(line)
	else: 
		print 'No job Queue file found'
		return []
	
	return runNums
	
	
def findMissing(fileList, inputList):
	runNums = []
	for file in fileList:
		f = file.replace('_Master.dat','')
		idx = f[::-1].find('_')
		try: n = int(f[-idx::])
		except: continue
		runNums.append(n)
	
	N = 0
	for file in inputList:
		f = file.replace('.dat','')
		idx = f[::-1].find('_')
		try: n = int(f[-idx::])
		except: continue
		if n > N: N = n
				
	N += 1
	#print N, 'runs found'
	
	missing = []
	for i in xrange(N):
		if i not in runNums: missing.append(i)
	
	#if len(missing) > 0: print 'List of runs missing:', ','.join([str(item) for item in missing])
	return missing
	

def keepWorking(path = './', logsPath = None, jobnumber = None, jobFile = 'mafot.sbatch'):
	from subprocess import call,check_output
	import getpass,time
	user = getpass.getuser()
	if not path[-1] == '/': path += '/'	
	
	# wait until no more of your jobs in 'preemptable' queue
	print 'Waiting for jobs to finish...'
	bla = check_output(['squeue','-u',user])
	if jobnumber is None: jobnumber = 'preemptab' # use queue keyword instead
	while(jobnumber in bla):
		time.sleep(600)		# wait 10 minutes; no CPU usage
		bla = check_output(['squeue','-u',user])
	
	runNums = findIncomplete(path, logsPath, jobFile = jobFile)
	while(len(runNums) > 0):
	
		# resubmit jobs
		bla = check_output(['sbatch', path + 'rerun_' + jobFile])
		jobnumber = bla.strip().split()[3]
		print 'Submitted batch job', jobnumber
		time.sleep(60)	# wait 60 sec
	
		# wait until no more of your jobs in 'preemptable' queue
		print 'Waiting for jobs to finish...'
		bla = check_output(['squeue','-u',user])
		while(jobnumber in bla):
			time.sleep(600)		# wait 10 minutes; no CPU usage
			bla = check_output(['squeue','-u',user])
		
		# check new incomplete list
		runNums = findIncomplete(path, logsPath, jobFile = jobFile)
	
	return
	
	
def restartJob(batchfile, jobnumber = None, logsPath = None):
	from subprocess import call,check_output
	import getpass,time
	user = getpass.getuser()

	idx = batchfile[::-1].find('/')	# returns location of last '/' in batchfile or -1 if not found
	if(idx == -1): path = './'
	else: 
		path = batchfile[0:-idx]	# path with a final '/'
		batchfile = batchfile[-idx::]

	#if not path[-1] == '/': path += '/'
	if logsPath is None: 
		logsPath = path
	else: 
		if not logsPath[-1] == '/': logsPath += '/'
		
	# read batchfile and get mpirun call
	if os.path.isfile(path + batchfile):
		with open(path + batchfile) as f: lines = f.readlines()
	else: 
		print 'Job file not found'
		return
		
	for line in lines[::-1]:
		if 'mpirun' in line: job = line.strip().split()
	
	# get mafot tool
	if 'laminar' in job[3]: task = 'lam'
	if 'foot' in job[3]: task = 'foot'
	if 'plot' in job[3]: task = 'plot'
	
	# get file tag
	chk = [True for item in job]
	for i in xrange(4): chk[i] = False
	for i,item in enumerate(job):
		if (item[0] == '_'): chk[i] = False
		if (item[0] == '-'): 
			chk[i] = False
			chk[i+1] = False
	
	tag = [job[i] for i in xrange(len(job)) if chk[i]][0]

	print 'Waiting for job to finish...'
	bla = check_output(['squeue','-u',user])
	if jobnumber is None: jobnumber = 'preemptab' # use queue keyword instead
	while(jobnumber in bla):
		time.sleep(60)		# wait 1 minute; no CPU usage
		bla = check_output(['squeue','-u',user])
	
	while(1):
		fileList = [f for f in os.listdir(logsPath) if(('_Master.dat' in f) & (tag in f) & (task in f))]
		incompleteRuns = []
		for file in fileList:
			with open(logsPath + file) as f: lines = f.readlines()
			if 'Program terminates normally' in lines[-1]: continue
			else: incompleteRuns.append(file)
			
		if len(incompleteRuns) == 0:
			print 'Run complete.'
			break
		else: 
			print 'Restarting: ', incompleteRuns[0].replace('log_','').replace('_Master.dat','')
		
		bla = check_output(['sbatch', path + batchfile])
		jobnumber = bla.strip().split()[3]
		print 'Submitted batch job', jobnumber
		time.sleep(10)	# wait 10 sec for job to start

		print 'Waiting for job to finish...'
		bla = check_output(['squeue','-u',user])
		while(jobnumber in bla):
			time.sleep(60)		# wait 1 minute; no CPU usage
			bla = check_output(['squeue','-u',user])
	
	
if __name__ == '__main__':
	import argparse
	import textwrap
	parser = argparse.ArgumentParser(description = 'Find and restart incomplete MAFOT runs', 
                formatter_class = argparse.RawDescriptionHelpFormatter,
                epilog = textwrap.dedent('''\
                Examples: restartMafot.py
                          restartMafot.py -p'''))
	
	parser.add_argument('-d', '--dir', help = "set working dir", type = str, default = './')
	parser.add_argument('-f', '--find', help = 'Just find the incomplete runs, job arrays only', action = 'store_true', default = False)
	parser.add_argument('-i', '--id', help = 'keep resubmitting job with initial number ID until complete. Single run or job arrays.', type = str, default = None)
	parser.add_argument('-j', '--job', help = 'keep resubmitting jobfile JOB until complete, default is mafot.sbatch', type = str, default = None)
	parser.add_argument('-l', '--logs', help = "set log dir, if different from working dir", type = str, default = None)
	parser.add_argument('-p', '--persist', help = 'keep resubmitting until all jobs complete, job arrays only', action = 'store_true', default = False)
	args = parser.parse_args()
	
	if args.job is None: job = 'mafot.sbatch'
	else: job = args.job
	
	if args.persist:
		keepWorking(args.dir, args.logs, args.id, job)
	elif args.find:
		_ = findIncomplete(args.dir, args.logs, job)
	else:
		restartJob(job, args.id, args.logs)
