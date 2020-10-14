# coding=utf-8
# Code to compute th spectrum of DM recoils in detectors for a given DM mass

import sys,os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import scipy.integrate as scint
from neutrino_events import *
from constants import *
import math


def int_new(e,mdm,AA,ZZ,funn,iso=iso_1,v_esc=v_esc,v_0=v_0,fp=1,fn=1,inel='none',delta=0,m_phi=0):
	#mdm mass of DM particle, AA mass number, Z electron number, funn acceptance as a function of energy for detector of interest
	#new improved version of integral with options inel=yes if inelastic DM, delta mass difference (<0 if exothermic)
	#can change  fp and fn for isospin violating DM, m_phi option to have a low mass mediator
	#this analytical solution of integral for maxwellian speed distribution from: arXiv: 1605.05098 (Geng and al. 2016) (also derived by David)
	# e is numpy array comtaining all values of e at which to evaluate speed integral
	dR=0
	for j in range(len(AA)): #loop over different elements or isotopes if needed
		dens_target=6.022e23/(AA[j]*1e-3)
		m_N=AA[j]*m_p
		red_p=mdm*m_p/(mdm+m_p)
		red_n=mdm*m_N/(mdm+m_N)

		f=nuc_F(e,AA[j]) #nuclear form factor

		if inel=='none':
			v_min=((m_N*e)/(2*red_n**2))**(0.5)

		elif inel=='yes': #minimal speed for inelastic DM, from arXiv: 1605.05098 (Geng and al. 2016)
			v_min=(2*m_N*e)**(-0.5)*np.abs(delta+m_N*e/red_n)

		coup=(fp*ZZ[j]+fn*(AA[j]-ZZ[j])) 
		C=coup**2*red_n**2/red_p**2
		x=v_min/v_0
		y=v_earth/v_0
		z=v_esc/v_0

		ind_1 = x<abs(y-z) 
		ind_2 = (x>=abs(y-z)) * (x<(y+z))

		N_esc=scipy.special.erf(z)-2*z/((math.pi)**0.5)*np.exp(-z**2)

	
		eta_1=1/(2*N_esc*v_0*y)*(scipy.special.erf(x+y)-scipy.special.erf(x-y)-4/(math.pi**(0.5))*y*np.exp(-z**2))
	
		eta_2=1/(2*N_esc*v_0*y)*(scipy.special.erf(z)-scipy.special.erf(x-y)-2/(math.pi**(0.5))*(y+z-x)*np.exp(-z**2))
	
		eta = eta_1 * ind_1 + eta_2 * ind_2
	

		dR+=iso[j]*cc*coup**2*funn(e)*eta*dens_target*f**2*m_N/(mdm*2*red_p**2)

		if m_phi != 0: #if mediator mass given (for 2 fermionic majorana DM particles with small mass difference coupled with dark photon of mass m_phi)
			#this comes from arXiv 1605.05098 (Geng and al. 2016) and arXiv: 1412.6220 (Li and al. 2015)
			q2=2*m_N*e
			qref2=4*red_p**2*v_ref**2
			fac=(1+qref2/(m_phi**2))/(1+q2/(m_phi**2))**2
			dR*=fac

		return dR


def int_other_dist(e,mdm,AA,ZZ,fun_int,funn,iso=iso_1,v_esc=v_esc,v_0=v_0,fp=1,fn=1,inel='none',delta=0,m_phi=0):
	#compute distribution of DM events with a speed distribution that is ot Maxwellian
	#fun_int speed distribution in lab frame/v to integrate
	dR=0
	for j in range(len(AA)):
		dens_target=6.022e23/(AA[j]*1e-3)
		m_N=AA[j]*m_p
		red_p=mdm*m_p/(mdm+m_p)
		red_n=mdm*m_N/(mdm+m_N)

		f=nuc_F(e,AA[j])	

		if inel=='none':
			v_min=((m_N*e)/(2*red_n**2))**(0.5)

		elif inel=='yes': #minimal speed for inelastic DM, from arXiv: 1605.05098 (Geng and al. 2016)
			v_min=(2*m_N*e)**(-0.5)*np.abs(delta+m_N*e/red_n)

		coup=(fp*ZZ[j]+fn*(AA[j]-ZZ[j])) 
		C=coup**2*red_n**2/red_p**2
		eta = scint.quad(fun_int,v_min,1000)[0]

		dR+=iso[j]*cc*coup**2*funn(e)*eta*dens_target*f**2*m_N/(mdm*2*red_p**2)

		if m_phi != 0: #if mediator mass given (for 2 fermionic majorana DM particles with small mass difference coupled with dark photon of mass m_phi)
			#this comes from arXiv 1605.05098 (Geng and al. 2016) and arXiv: 1412.6220 (Li and al. 2015)
			q2=2*m_N*e
			qref2=4*red_p**2*v_ref**2
			fac=(1+qref2/(m_phi**2))/(1+q2/(m_phi**2))**2
			dR*=fac

		return dR


def eta_0(v_min,mdm,v_esc=v_esc,v_0=v_0,dens_dm=dens_dm,fp=1,fn=1,inel='none',delta=0):
	#mdm mass of DM particle, AA mass number
	#new improved version of integral with options inel=yes if inelastic DM, delta mass difference (<0 if exothermic)
	#can change  fp and fn for isospin violating DM, m_phi option to have a low mass mediator
	#this analytical solution of integral for maxwellian speed distribution from: arXiv: 1605.05098 (Geng and al. 2016) (also derived by David)
	#THIS GIVES THE VALLUE OF n(vmin) FOR THE SHM TO COMPARE WITH SHM INDEPENDANT EXCLUSION CURVES

	x=v_min/v_0
	y=v_earth/v_0
	z=v_esc/v_0

	N_esc=scipy.special.erf(z)-2*z/((math.pi)**0.5)*np.exp(-z**2)

	if x<abs(y-z):
		eta=1/(2*N_esc*v_0*y)*(scipy.special.erf(x+y)-scipy.special.erf(x-y)-4/(math.pi**(0.5))*y*np.exp(-z**2))

	if (x > abs(y-z)) and (x<(y+z)):
		eta=1/(2*N_esc*v_0*y)*(scipy.special.erf(z)-scipy.special.erf(x-y)-2/(math.pi**(0.5))*(y+z-x)*np.exp(-z**2))

	if x>(y+z):
		eta=0

	dR=cc*eta*dens_dm/(mdm)


	return dR



vmax=1000 #maximal speed in lab frame above which speed distribution 0
#integrate velocity dist/v already in halo frame (SHM ++)

#upload distributions and fit them
arr_v = np.load('speed_beta_eta.npy') #precomputed values of f(v)/v for different values of beta and eta
v_tab = np.linspace(0,1000,200)

vv_tab = np.linspace(0,1000,20000)
dv = 1000/20000

dic = {} 
n_i = len(arr_v[:,0,0])
n_j = len(arr_v[0,:,0])

dic_v={} #dictionary to contain function of f(v)/v for different beta and eta values

for i in range(n_i):
	dic_v[i] = {}
	for j in range(n_j):

		arr_v[i,j] = np.nan_to_num(arr_v[i,j])

		fun_v = scipy.interpolate.interp1d(v_tab,arr_v[i,j],kind='linear',bounds_error = False, fill_value = 0)

		dic_v[i][j] = fun_v(vv_tab)


def int_simpson(v,dist):
	#takes in array with speeds and array with distribution to integrate and returns value of integral using simpson method
	#v need to be equally spaced and ordered!
	if len(v)<4:
		sum_simpson = 0
	else:
		if len(v)/2 == len(v)//2: #simpson method requires even number of subdivisions
			v = v[:-1]
			dist = dist[:-1]

		h = v[1] - v[0] #spacing

		sum_simpson = h/3 * (dist[0] + 4*np.sum(dist[1::2]) + 2*np.sum(dist[:-1][2::2]) + dist[-1])

	return sum_simpson



def int_1(e,red_m,AA,acc,v_fun):

	m_N=AA*m_p
	vmin=(e*m_N/(2*red_m**2))**(0.5)*cc/1000

	if vmax<vmin:
		sum_simp=0
	else:
		ind_min = np.argmin(np.abs(vv_tab-vmin)) #index of minimal speed

		sum_simp = int_simpson(vv_tab[ind_min:],v_fun[ind_min:])

	
	return acc(e)*nuc_F(e,AA)**2*sum_simp

def int_irreg(e,m,AA,acc,beta_ind,eta_ind):
	#e er energy, m DM mass, AA mass number of target material
	#v_fun function of f(v)/v in halo frame
	#beta_ind integer between 0-5 corresponding to a value of beta for SHM++ 0.82-0.97
	#eta_int integer between 0-7 corresponding to value of eta for SHM++ 0.1-0.8

	m_N=AA*m_p
	red_m=m*m_N/(m+m_N)
	red_p=m*m_p/(m+m_p)

	v_fun = dic_v[beta_ind][eta_ind]

	dens_target=6.022e23/(AA*1e-3) 

	int_e=int_1(e,red_m,AA,acc,v_fun) #integral over energies
	R=m_N*AA**2*dens_target/(2*m*red_p**2)*cc**2*int_e/1000

	return R





"""
# For distribution known is sun/earth frame (distributions from papers)
#import np file
v_data =  np.load('vdata.npy')
deltav = v_data[1,0] - v_data[0,0]


oldshm_max = v_data[np.sum(v_data[:,1]!=0)-1,0]
lisanti_max = v_data[np.sum(v_data[:,2]!=0)-1,0]
shmplus_max = v_data[np.sum(v_data[:,3]!=0)-1,0]


def v_int(fun_v,v_min,v_esc):
	#ref = np.where(label_tab==fun_v)
	ref_min = np.argmin(np.abs(v_data[:,0]-v_min)) #index of minimum v value for integral
	ref_max = np.argmin(np.abs(v_data[:,0]-v_esc)) #index of max v value

	int_tab = v_data[fun_v,ref_min:ref_max]

	if len(int_tab)<3:
		int_res = 0
	elif len(int_tab)%2 == 0:
		#print(len(int_tab))
		int_tab = np.concatenate((int_tab,np.array([0])))
		int_res = scint.romb(int_tab,deltav)
	else:
		#print(len(int_tab))
		int_res = scint.romb(int_tab,deltav)

	return int_res


def int_sun(e,red,AA,acc,v_esc,fun_v):
	m_N=AA*con.m_p
	vmin=(e*m_N/(2*red**2))**(0.5)

	if vmin>v_esc:
		sum_v = 0
	else:
		sum_v = v_int(fun_v,vmin,v_esc)
	return acc(e)*nuc_F(e,AA)**2*sum_v



def int_dm_sun(e,m,AA,acc,dens_dm=dens_dm,v_esc=oldshm_max,fun_v=1,const=1):
	m_N=AA*con.m_p
	red_m=m*m_N/(m+m_N)
	red_p=m*con.m_p/(m+con.m_p)

	dens_target=6.022e23/(AA*1e-3) 

	int_e=int_sun(e,red_m,AA,acc,v_esc,fun_v) #integral over energies
	R=dens_dm*m_N*AA**2*dens_target/(2*m*red_p**2)*con.cc*int_e/const

	return R
"""