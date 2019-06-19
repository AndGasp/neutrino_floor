#!/usr/bin/env python


##########################################################################
##  Based on:
##  goUniversal.py
##
##  A python script to run or submit jobs for the common use cases
##  of the IMSRG++ code. We check whether there is a pbs or slurm
##  scheduler, assign the relevant input parameters, set names
##  for the output files, and run or submit.
##  						-Ragnar Stroberg
##  						TRIUMF Nov 2016
##
##  Now:
##  go_oak.py
##  Changed for customized use by Andrea Gaspert (Apr-Jun 2019)
##  To compute neutrino floors for various experiments 
######################################################################

from os import path,environ,mkdir,remove
from sys import argv
from subprocess import call,PIPE
from time import time,sleep
from datetime import datetime
import argparse
import random
import numpy as np


### Check to see what type of batch submission system we're dealing with
BATCHSYS = 'NONE'
if call('type '+'qsub', shell=True, stdout=PIPE, stderr=PIPE) == 0: BATCHSYS = 'PBS'
elif call('type '+'srun', shell=True, stdout=PIPE, stderr=PIPE) == 0: BATCHSYS = 'SLURM'

### The code uses OpenMP and benefits from up to at least 24 threads
NTHREADS=32
exe = 'python run_floor_multi.py'

mail_address = 'andrea.gaspert@gmail.com' #YOUR EMAIL HERE!!!!! :)

### Make a directory for the log files, if it doesn't already exist
if not path.exists('nu_floor_log_multi'): mkdir('nu_floor_log_multi')


### Create the 'script' that we need for execution
FILECONTENT = """#!/bin/bash
#PBS -N %s
#PBS -q oak
#PBS -d %s
#PBS -l walltime=00:60:00
#PBS -l nodes=1:ppn=%d
#PBS -l vmem=10gb
#PBS -m ae
#PBS -M %s
#PBS -j oe
#PBS -o nu_floor_log_multi/%s.o
#PBS -e nu_floor_log_multi/%s.e
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=%d
%s
"""

### Loop parameters
batch_mode = False

#TRANSMIT ARGS TO PYTHON PROGRAM
### ARGS is a (string => string) dictionary of input variables that are passed to the main program
ARGS  =  {}
#element in detector, if no other detector properties specified, code runs with ideal detector properties for this element
ARGS['typ'] = 'Xe' 
# search for optimal cross-section or optimal exposure "cs" or 'exp'
# if mod == 'exp', cross-section needs to be specified
# if mod == 'cs' and no exposure is specified, exposure will be of 1e4 TY for Xenon or 1e7 TY for Argon (exp. to reach 'hard floor')
ARGS['mod'] = 'cs' 

ARGS['ER'] = True

#uncomment for quick test
ARGS['n_mc'] = 30
ARGS['n_dicho'] = 10


m_tab = np.ones(30)*100 #values of DM masses
exp_tab = np.logspace(-2,7,30)

ARGS['speed_dist'] = 'reg'



m_values = m_tab.tolist()
exp_values = exp_tab.tolist()

ARGS['m'] = m_values
ARGS['exp'] = exp_values


jobname = 'floor_{}_test_exposures_ER'.format(ARGS['typ'])
logname = jobname +"{:.3}".format(13*random.random()+random.random()) + datetime.fromtimestamp(time()).strftime('_%y%m%d%H%M.log')
cmd = ' '.join([exe] + ['%s=%s'%(x,ARGS[x]) for x in ARGS])

print(cmd)


if batch_mode==True:

	print("Submitting to cluster...")
	sfile = open(jobname+'.batch','w')
	sfile.write(FILECONTENT%(jobname,environ['PWD'],NTHREADS,mail_address,logname,logname,NTHREADS,cmd))
	sfile.close()
	call(['qsub', jobname+'.batch'])

else:
	call(cmd.split())

"""
ARGS['ER'] = True

for e in exp_tab:
	ARGS['exp'] = e

	jobname = 'floor_{}_test_pp_ER'.format(ARGS['typ'])
	logname = jobname +"{:.3}".format(13*random.random()+random.random()) + datetime.fromtimestamp(time()).strftime('_%y%m%d%H%M.log')
	cmd = ' '.join([exe] + ['%s=%s'%(x,ARGS[x]) for x in ARGS])

	print(cmd)


	if batch_mode==True:

		print("Submitting to cluster...")
		sfile = open(jobname+'.batch','w')
		sfile.write(FILECONTENT%(jobname,environ['PWD'],NTHREADS,mail_address,logname,logname,NTHREADS,cmd))
		sfile.close()
		call(['qsub', jobname+'.batch'])

	else:
		call(cmd.split())

		ARGS['ER'] = True


ARGS['typ'] = 'Ar' 
ARGS['ER'] = True

for e in exp_tab:
	ARGS['exp'] = e

	jobname = 'floor_{}_test_pp_ER'.format(ARGS['typ'])
	logname = jobname +"{:.3}".format(13*random.random()+random.random()) + datetime.fromtimestamp(time()).strftime('_%y%m%d%H%M.log')
	cmd = ' '.join([exe] + ['%s=%s'%(x,ARGS[x]) for x in ARGS])

	print(cmd)


	if batch_mode==True:

		print("Submitting to cluster...")
		sfile = open(jobname+'.batch','w')
		sfile.write(FILECONTENT%(jobname,environ['PWD'],NTHREADS,mail_address,logname,logname,NTHREADS,cmd))
		sfile.close()
		call(['qsub', jobname+'.batch'])

	else:
		call(cmd.split())
"""
