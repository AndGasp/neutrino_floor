#########################
# Andrea Gaspert
# 11/10/2020
# Plot test and neutrino floor example code
##########################
import numpy as np
import matplotlib.pyplot as plt
from dm_events import *
from optim import *
from constants import *


"""
########## DM dR/DdEr for different targets ###############################
m_tab = np.array([25,50,100])*1e6
c_tab = np.array([1e-45,1e-45,1e-45])*1e-4
AA_tab = np.array([A_xenon_s,A_xenon_s,A_xenon_s])
ZZ_tab = np.array([Z_xenon,Z_xenon,Z_xenon])
col_tab = ['r','g','b']
lab_tab = ['25 GeV', '50 GeV', '100 GeV']

emin = 0
emax = 100
nbin = 15

density = 0.3e12
v0 = 220e3/cc
vesc = 544e3/cc

dm_params = np.array([vesc,v0,0,0])
speed_dist = 'reg'
dist = empty
acc = acc_100 

e_tab = np.linspace(emin,emax,1000)

exposure =  24 * 3600

for i in range(len(m_tab)):
	dR = int_new(e_tab,m_tab[i],AA_tab[i],ZZ_tab[i],acc,iso=iso_1,v_esc=vesc,v_0=v0,fp=1,fn=1,inel='none',delta=0,m_phi=0) * exposure * dens_dm_old * c_tab[i]
	plt.plot(e_tab,dR,col_tab[i],label = lab_tab[i])

plt.legend(loc='best')
plt.yscale('log')
plt.show()
"""
##### Test statistics for various exposures ####################################

m = 30e6
c_tab = np.array([1e-47,1e-48,1e-49,1e-50,1e-51])*1e-4
AA = A_xenon_s
ZZ= Z_xenon
lab_tab = ['1e-47', '1e-48', '1e-49', '1e-50', '1e-51']
color_tab = ['r','g','b','m','c']

emin = 5
emax = 100
nbin = 20

density = 0.3e12
v0 = 220e3/cc
vesc = 544e3/cc

dm_params = np.array([vesc,v0,0,0])
speed_dist = 'reg'
dist = empty
acc = acc_100 

resol = 0.05

exp_tab = np.logspace(3,12,200)*365*24*3600

neu_argon = np.load('dist_neutrino_argon_new.npy')
neu_xenon = np.load('dist_neutrino_xenon_new.npy')

er_rej = 1e12

z_array = np.zeros((len(c_tab),len(exp_tab)))
n_array = np.zeros((len(c_tab),len(exp_tab)))

for i in range(len(c_tab)):
	for j in range(len(exp_tab)):
		z_array[i,j] , n_array[i,j] = stats_run(m,c_tab[i],emin,emax,resol,neu_xenon,AA,ZZ,er_rej,acc,exp_tab[j],speed_dist,dist,Quench_true=0)

	plt.plot(exp_tab/24/3600/365,z_array[i,:]**2,color_tab[i]+'-', label = lab_tab[i])

for i in range(len(c_tab)):
	for j in range(len(exp_tab)):
		z_array[i,j] , n_array[i,j] = stats_run(m,c_tab[i],emin,emax,resol,neu_argon,A_argon_s,Z_argon,er_rej,acc,exp_tab[j],speed_dist,dist,Quench_true=0)

	plt.plot(exp_tab/24/3600/365,z_array[i,:]**2,color_tab[i]+':')

plt.legend(loc='best')
plt.yscale('log')
plt.xscale('log')
plt.show()


"""



bins_tab = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

for i in range(4):
	events = dm_bins(emin,emax,nbin,m_tab[i],acc_100,AA,ZZ,dm_params,speed_dist,dist)*10000*24*365*3600*1000*c_tab[i]*density
	plt.plot(bins_tab,events,'.',label='{:.2e} GeV,{:.2e} cm^2'.format(m_tab[i]/1e6,c_tab[i]*1e4))

plt.legend(loc='best')
plt.xlabel('Bin')
plt.ylabel('events')
plt.yscale('log')
plt.show()


search_min = -55
search_max = -47
niter = 12
n_run = 1
emin = 55
emax = 100
nbin = 50
neu_arr = np.load('dist_neutrino_argon_new.npy')
AA = A_argon_s
ZZ = Z_argon
er_rej = 1e9
acc = acc_100 
speed_dist = 'reg'
dist = empty
exp = 10000

n_try = 1
print('%12 iteration in bi, 200 mc runs, Xe Ideal (1-90 keV) with 1E3 and 1000TY')

for m in m_tab:
	av=0
	for n in range(n_try):
		cs, n_dm = dicho_search_sig(search_min,search_max,niter,n_run,m*1e6,emin,emax,nbin,neu_arr,AA,ZZ,er_rej,acc,exp*24*3600*365*1000,speed_dist,dist, mc=0)
		cross_section = 10**(cs+4)
		av+=cross_section
	cs_av = av/n_try
	print('###m={} GeV, cs={} cm2, exposure={} Ty\n dm events = {}'.format(m,cs_av,exp,n_dm))

"""