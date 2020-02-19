#!/usr/bin/env python
# coding=utf-8
#################################
# Andrea Gaspert 04/2019
# Main function to run code on cluster
#code to run up to 32 mass points to find the neutrino floor simultaneously 
#on 1 node (oak), or whichever number of cores are available on another cluster/device)
###################################

#from platform import python_version
#print(python_version())

# TO OD: OPTIMIZE THIS USING BUILT-IN ARG FUNCTION?

import sys,os
import numpy as np
from constants import *
from acceptance import acc_100
from optim import *
import time
from extract_file import text_to_numpy
from figure_plot import plot_floor
from multiprocessing import Pool


#Use command line arguments: possible args are m,cs,exp,E_min,E_max,ER,typ,mod,n_mc,n_dicho,n_bin,speed_dist (pp to profile over SHM++ parameters, else to use other speed distribution), dist to specify other speed distribution (in halo frame, integrated by hand, VERY LONG TO COMPUTE)

fullCmdArguments = sys.argv #all arguments

# - further arguments
argumentList = fullCmdArguments[1:] #arguments after the name

def keep_going(i):
	try:
		next = argumentList[i+1]
		add = 1
	except:
		add = 0
	return add

#print(argumentList)
arg = ' '
i=0
while i<(len(argumentList)-1):
	#print(i)
	if 'm=[' in argumentList[i]:
		m_start = i
		arg_last = ' '

		while ']' not in arg_last:
			i+=1
			arg_last = argumentList[i]
	
		m_end = i
		m_s = argumentList[m_start]

		for j in range(m_start+1,m_end+1):
			m_s += argumentList[j]

		m_string = m_s.split('=')[1]
		#print(m_string)
		m_tab = [float(s) for s in m_string[1:-1].split(',')]

		i+=keep_going(i)

	if 'cs=[' in argumentList[i]:
		cs_start = i
		arg_last = ' '
		i+=1
		while ']' not in arg_last:
			i+=1
			arg_last = argumentList[i]

		cs_end = i

		cs_s = argumentList[cs_start]
		for j in range(cs_start+1,cs_end+1):
			cs_s += argumentList[j]

		cs_string = cs_s.split('=')[1]
		#print(cs_string)
		cs_tab = [float(s) for s in cs_string[1:-1].split(',')]

		i+=keep_going(i)

	if 'exp=[' in argumentList[i]:
		exp_start = i
		arg_last = ' '

		while ']' not in arg_last:
			i+=1
			arg_last = argumentList[i]
	
		exp_end = i
		exp_s = argumentList[exp_start]
		for j in range(exp_start+1,exp_end+1):
			exp_s += argumentList[j]

		exp_string = exp_s.split('=')[1]
		#print(exp_string)
		exp_tab = [float(s) for s in exp_string[1:-1].split(',')]

		i+=keep_going(i)

	if 'E_min=' in argumentList[i]:
		E_min = float(argumentList[i].split('=')[1])
		i+=keep_going(i)

	if 'E_max=' in argumentList[i]:
		E_max = float(argumentList[i].split('=')[1])
		i+=keep_going(i)

	if 'ER=' in argumentList[i]:
		ER = argumentList[i].split('=')[1] == "True"
		i+=keep_going(i)

	if 'er_rej=' in argumentList[i]:
		eff_er = argumentList[i].split('=')[1]
		i+=keep_going(i)

	if 'typ=' in argumentList[i]:
		typ = argumentList[i].split('=')[1]
		i+=keep_going(i)

	if 'mod=' in argumentList[i]:
		mod = argumentList[i].split('=')[1]
		i+=keep_going(i)

	if 'n_mc=' in argumentList[i]:
		n_mc = int(argumentList[i].split('=')[1])
		i+=keep_going(i)

	if 'n_dicho=' in argumentList[i]:
		n_dicho = int(argumentList[i].split('=')[1])
		i+=keep_going(i)

	if 'n_bin=' in argumentList[i]:
		n_bin = int(argumentList[i].split('=')[1])
		i+=keep_going(i)

	if 'speed_dist=' in argumentList[i]:
		speed_dist = argumentList[i].split('=')[1]
		i+=keep_going(i)

	if 'dist=' in argumentList[i]:
		dist = argumentList[i].split('=')[1]
		i+=keep_going(i)

	if 'save=' in argumentList[i]:
		save_name = argumentList[i].split('=')[1]
		i+=keep_going(i)

	if 'plot_floor=' in argumentList[i]:
		plot_floor_truth = argumentList[i].split('=')[1] == "True"
		i+=keep_going(i)


def empty(v):
	return v

try:
	det = typ
except NameError:
	print('Type of detector needs to be specified in command (typ=Xe or typ=Ar)')

try:
	print('mode='+mod)
except:
	print('A mode needs to be specified in command (mod=cs or mod=exp)')

try:
	print('save name='+save_name)
except NameError :
	print('Specify name of file to save info to (save_name=)')

try:
	print('plot={}'.format(plot_floor_truth))
except NameError :
	plot_floor_truth = True


if det == "Ar":
	print('Argon Detector')
	AA = A_argon_s
	ZZ = Z_argon

	neu_arr = np.load('dist_neutrino_argon_long.npy')


	try: 
		print('E_min = {} keV'.format(E_min)) 
	except NameError : 
		E_min = E_min_Ar
		print('E_min = {} keV'.format(E_min)) 
	try: 
		print('E_max = {} keV'.format(E_max)) 
	except NameError :
		E_max = E_max_Ar
		print('E_max = {} keV'.format(E_max)) 
	try:
		print('ER = {}'.format(ER))
	except NameError:
		print('Presence or not of ER in calculations needs to be specified!')
	else:
		if ER:
			try:
				print('eff_er = {}'.format(eff_er))
			except NameError:
				eff_er = rej_er_Ar
				print('eff_er = {}'.format(eff_er))
		if not ER:
			eff_er = 1e15

	if mod == 'cs':
		try:
			print('Exposure = {} Ty'.format(exp_tab))
		except NameError:
			exp_tab = [1e7] #have to check
			print('Exposure = {} Ty'.format(exp_tab))


if det == "Xe":
	print('Xenon Detector')
	AA = A_xenon_s
	ZZ = Z_xenon

	neu_arr = np.load('dist_neutrino_xenon_long.npy')


	try: 
		print('E_min = {} keV'.format(E_min)) 
	except NameError : 
		E_min = E_min_Xe
		print('E_min = {} keV'.format(E_min)) 
	try: 
		print('E_max = {} keV'.format(E_max)) 
	except NameError :
		E_max = E_max_Xe
		print('E_max = {} keV'.format(E_max)) 
	try:
		print('ER = {}'.format(ER))
	except NameError:
		print('Presence or not of ER in calculations needs to be specified!')
	else:
		if ER:
			try:
				print('eff_er = {}'.format(eff_er))
			except NameError:
				eff_er = rej_er_Xe
				print('eff_er = {}'.format(eff_er))
		if not ER:
			eff_er = 1e15

	if mod == 'cs':
		try:
			print('Exposure = {} Ty'.format(exp_tab))
		except NameError:
			exp_tab = 1e4 #have to check
			print('Exposure = {} Ty'.format(exp_tab))


try:
	print('n_mc = {}'.format(n_mc))
except NameError:
	n_mc = 50
	print('n_mc = {}'.format(n_mc))

try:
	print('n_dicho = {}'.format(n_dicho))
except NameError:
	n_dicho = 8
	print('n_dicho = {}'.format(n_dicho))

try:
	print('n_bin = {}'.format(n_bin))
except NameError:
	n_bin = 50
	print('n_bin = {}'.format(n_bin))

try:
	print('Masses : {}'.format(m_tab))
except NameError:
	print('DM mass list to compute not specified!')

try:
	print('Speed distribution ='+speed_dist)
except NameError:
	speed_dist = 'reg'
	print('Regular SHM speed distribution')

try:
	print('Distribution ='+dist)
except NameError:
	dist = empty


if mod == "cs": 

	list_info_cs = [n_dicho,n_mc,E_min,E_max,n_bin,neu_arr,AA,ZZ,eff_er,speed_dist,dist,save_name]

	def multi_task(task_tab,co_args=list_info_cs):
		n_dicho = co_args[0]
		n_mc = co_args[1]
		E_min = co_args[2]
		E_max = co_args[3]
		n_bin = co_args[4]
		neu_arr = co_args[5]
		AA = co_args[6]
		ZZ = co_args[7]
		eff_er = co_args[8]
		speed_dist = co_args[9]
		dist = co_args[10]
		save_name = co_args[11]

		m = task_tab[0]
		exp = task_tab[1]

		if m<20: #TO-DO, change to avoid this knit-picking?
			sigmin = -52
			sigmax = -44

		if m>=20:
			sigmin = -56
			sigmax = -47

		sigma_floor = dicho_search_sig(sigmin,sigmax,n_dicho,n_mc,m*1e6,E_min,E_max,n_bin,neu_arr,AA,ZZ,eff_er,acc_100,exp*1000*365*24*3600,speed_dist,dist)
		cross_section = 10000*10**sigma_floor

		info_line = '###m={} GeV, cs={} cm2, exposure={} Ty \n'.format(m,cross_section,exp)

		print(info_line)

		#CHANGE TO CRASH EXISTING FILE IF ALREADY EXISTS

		file1 = open(save_name+".txt","a+") 
		file1.write(info_line)
		file1.close()

		return m,cross_section,exp

	n_tasks = len(m_tab)
	print('number of tasks:{}'.format(n_tasks))
	m_pairs = np.zeros((n_tasks,2))
	m_pairs[:,0] = m_tab

	print(exp_tab)

	try: 
		if len(exp_tab) != n_tasks:
			print('List of exposures needs to have the same length as list of DM masses (or only contain one value)')
		m_pairs[:,1] = exp_tab

	except:
		m_pairs[:,1] = exp_tab*np.ones(n_tasks)

	tasks = []
	for j in range(n_tasks):
		tasks.append(m_pairs[j,:].tolist())


if mod == 'exp': 

	try:
		print(cs_tab)

	except NameError:
		print('You need to specify cross-section list for the given mass point list(m^2)!')

	list_info_exp = [n_bin,n_mc,E_min,E_max,n_bin,neu_arr,AA,ZZ,eff_er,speed_dist,dist,save_name]

	def multi_task(task_tab,co_args_exp=list_info_exp):

		m = task_tab[0]
		cs= task_tab[1]

		if cs>1e-40 or cs<1e-52:
			print('Cross-section is out of range (1e-40 - 1e-52 cm^2)')
			exp = 0

		else:
			n_dicho = co_args_exp[0]
			n_mc = co_args_exp[1]
			E_min = co_args_exp[2]
			E_max = co_args_exp[3]
			n_bin = co_args_exp[4]
			neu_arr = co_args_exp[5]
			AA = co_args_exp[6]
			ZZ = co_args_exp[7]
			eff_er = co_args_exp[8]
			speed_dist = co_args[9]
			dist = co_args[10]
			save_name = co_args[11]

			exp_min = -4
			exp_max = 10

			exp_floor = dicho_search_exp(exp_min,exp_max,n_dicho,n_mc,m*1e6,E_min,E_max,n_bin,neu_arr,AA,ZZ,eff_er,acc_100,cs/1e4,speed_dist,dist)

			info_line = '###m={} GeV, cs={} cm2, exposure={} Ty \n'.format(m,cs,exp_floor)

			print(info_line)

			file1 = open(save_name+".txt","a+") 
			file1.write(info_line)
			file1.close()


		return m,cs,exp_floor

	n_tasks = len(m_tab)
	print('number of tasks:{}'.format(n_tasks))
	m_pairs = np.zeros((n_tasks,2))
	m_pairs[:,0] = m_tab

	if len(cs_tab)>1:
		if len(cs_tab) != n_tasks:
			print('List of cross-sections needs to have the same length as list of DM masses (or only contain one value)')
		m_pairs[:,1] = cs_tab

	else:
		m_pairs[:,1] = cs_tab[0]*np.ones(n_tasks)
	tasks = []
	for j in range(n_tasks):
		tasks.append(m_pairs[j,:].tolist())


if __name__ == '__main__':

	#mark the start time
	total_tasks = len(m_tab) #number of mass points wanted
	startTime = time.time()
	#count number of cores available on the computer
	n_cpu = os.cpu_count()
	print('number of cpu used:{}'.format(n_cpu))

	#create pool with number of processes equal to the number of available cpu cores
	max_number_processes = n_cpu
	pool = Pool(max_number_processes)

	#map single_task to available Pool processes
	results = pool.map(multi_task, tasks)
	pool.close()
	pool.join()

	#mark the end time
	endTime = time.time()
	#calculate the total time it took to complete the work
	workTime =  endTime - startTime

	#print results
	print("The job took " + str(workTime) + " seconds to complete")

	#read txt file created and put information in numpy array
	
	results = text_to_numpy(save_name+".txt")

	#plot neutrino floor
	if plot_floor_truth: 
		plot_floor(results,save_name,exp_show=True)