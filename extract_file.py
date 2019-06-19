# coding=utf-8
# Andrea Gaspert 05/2019
# Code to find all file from the same 'job', compile information on neutrino floor and save it in numpy array

import numpy as np
import matplotlib.pyplot  as plt
import os

def extract_num(s):
	s1 = s.split('=')[1] #remove variable name
	print(s1)
	s2 = ''.join(reversed(s1)) #reverse order of string
	print(s2)
	s3 = s2[1:] #remove potential number in units (cm2)
	print(s3)

	i=0
	while True:
		try:
			val = int(s3[i])
			print(val)
			ind = i
			break
		except ValueError:
			i+=1

	s4 = ''.join(reversed(s3[ind:]))
	print(s4)
	num = float(s4)

	return num


directory = '.\nu_floor_log'

list_file = os.listdir(directory)

m = []
cs = []
exp = [] 

d_type = [('mass', float), ('cs', float), ('exp', float)] #create data type for neutrino floor points

d_list = []
i = 0
for filename in list_file:
	if "Xe" in filename and "floor" in filename and "noER" in filename:
		f = open(directory+filename)
		for line in f:
			if "###" in line: #finds line with final information
				line_elems = line[3:].split(', ') #slits 3 different elements (mass, cross section, exposure)	
				d_list.append((line_elems[0],line_elems[1].line_elems[2]))
				i+=1

#reorder mass points, since processes don't necessary finish in order

d_arr = np.array(d_list, dtype = d_type) #convert m list to array

print(d_arr)
d_ord =  np.sort(d_arr, order='mass') #arrange array element by mass order

print(d_ord)

name = "floor_xenon_ER.npy"
np.save(name, d_ord)