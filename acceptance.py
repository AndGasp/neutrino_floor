# coding=utf-8
#Acceptances and Quenching factors for various experiments (references included)
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

import constants as con



#ACCEPTANCE FOR DEAP CALCULATION====================================================


#Lindhard-birks Quenching factor as a function of energy for DEAP-3600 (from figure sent by Pietro - not published)
Q=np.array([0.23,0.25,0.26,0.27,0.28,0.29,0.3,0.305,0.31,0.315,0.32])
E_nr=np.array([10,15,20,30,40,60,80,100,120,140,160])

Q_ee=scipy.interpolate.interp1d(Q*E_nr,1/Q,kind='linear') #Quenching factor as a fun of EE energy interpolation

E_ee_deap=np.array([8,10.9,12.2,13.6,14.9,16.3,17.7,19,20.4,21.7,23.1,24.5,25.8,27.2,28.5,30,31.3,32.6,40])
en_k_deap=Q_ee(E_ee_deap)*E_ee_deap
alpha_k_deap=np.array([0,0.04,0.1,0.2,0.3,0.43,0.5,0.55,0.6,0.61,0.62,0.63,0.63,0.63,0.63,0.63,0.63,0.63,0.63]) #from arXiv:1707.08042 FIG.4 (2017)


acc_deap_1=scipy.interpolate.interp1d(en_k_deap,alpha_k_deap,kind='linear') #linear interpolation of acceptance as a function of Enr

def acc_deap(e):
	#acceptance as a function of Enr for DEAP-3600 detector (2017)
	emin=con.E_min_deap #ROI limits from 2017 paper as well
	emax=con.E_max_deap
	if e<emin or e>emax:
		r=0
		print('Recoil energy outside of accepted bounds for DEAP-3600')
	else:
		r=acc_deap_1(e)
	return r 

#ACCEPTANCE FOR XENON-1T CALCULATION====================================================
en_k_xenon=np.array([0,2,3,4.5,5.5,8.2,12,14,16,20,24,28,28,32,35,38,40,45,47,50,55,60,65])
alpha_k_xenon=np.array([0,0,0.1,0.2,0.27,0.6,0.8,0.82,0.85,0.855,0.86,0.855,0.85,0.82,0.8,0.75,0.7,0.45,0.3,0.12,0.05,0.005,0.00005]) #from arXiv: 1805.12562 (2018)


acc_xenon_1=scipy.interpolate.interp1d(en_k_xenon,alpha_k_xenon,kind='linear') #linear interpolation of acceptance as a function of Enr

def acc_xenon(e):
	emin=con.E_min_xenon
	emax=con.E_max_xenon
	if e<emin or e>emax:
		r=0
	else:
		r=acc_xenon_1(e)
	return r


#ACCEPTANCE FOR LUX CALCULATION====================================================
en_k_lux=np.array([0,1.1,2,3,4,5,6,7,10,20,25,30,40,50,60,70,80,90])
alpha_k_lux=np.array([0,5e-3,9e-2,0.3,0.45,0.62,0.72,0.8,0.9,0.98,0.995,1,0.8,0.3,4e-2,2e-3,1e-4,1e-5]) #from arXiv: 1608.07648 (2017)


acc_lux_1=scipy.interpolate.interp1d(en_k_lux,alpha_k_lux,kind='cubic')

def acc_lux(e):
	emin=con.E_min_lux
	emax=con.E_max_lux
	if e<emin or e>emax:
		r=0
	else:
		r=acc_lux_1(e)
	return r

#ACCEPTANCE FOR PANDAX CALCULATION====================================================
en_k_panda=np.array([0,1.1,2,3,4,5,6,7,8,9,10,13,15,18,20,30,40,50,60,70,80,90,100]) #from arXiv: 1607.07400 (2016)
alpha_k_panda=np.array([0,0.001,0.025,0.1,0.236,0.375,0.5,0.625,0.7,0.75,0.8,0.875,0.890,0.9,0.875,0.64,0.275,0.12,0.05,0.03,0.02,0.01,0.001])


acc_panda_1=scipy.interpolate.interp1d(en_k_panda,alpha_k_panda,kind='cubic')

def acc_panda(e):
	emin=con.E_min_panda
	emax=con.E_max_panda
	if e<emin or e>emax:
		r=0
	else:
		r=acc_panda_1(e)
	return r

#ACCEPTANCE FOR SUPERCDMS==============================================================

en_k_cdms=np.array([9,10,10.2,12,14,16,24,44,60,76,100,110]) #from arXiv: 1504.05871 (2015) fig. 10, black curve (classical analysis)
alpha_k_cdms=np.array([0,0,0.2,0.3,0.34,0.36,0.38,0.37,0.36,0.34,0.32,0.31])

acc_cdms_1=scipy.interpolate.interp1d(en_k_cdms,alpha_k_cdms,kind='linear')

def acc_cdms(e):
	emin=con.E_min_cdms
	emax=con.E_max_cdms
	if e<emin or e>emax:
		r=0
	else:
		r=acc_cdms_1(e)
	return r

def acc_100(e): #function for 100% acceptance
	return 1