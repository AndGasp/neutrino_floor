# coding=utf-8
# Code to produce "normalized" recoil energy distributions for different neutrino sources in different detector materials
# Only need to run once
import sys,os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import scipy.integrate as scint
from neutrino_events import *
from constants import *

absFilePath = os.path.abspath(__file__)                # Absolute Path of the module
fileDir = os.path.dirname(os.path.abspath(__file__))   # Directory of the Module
newPath = os.path.join(fileDir, 'data')   # Get the directory with neutrino flux information
sys.path.append(newPath)


# TYPICAL SOLAR NEUTRINOS==================================================
#8B NEUTRINOS
n_8b=np.genfromtxt('data/8benergytab.txt') #txt file containing 8B neutrino energy spectrum (energy in MeV, spectrum normalized to 1)
en_ns_8b0=n_8b[:,0]*1e3 #txt files for all solar neutrino fluxes taken from J. Bahcall's website: http://www.sns.ias.edu/~jnb/
f_ns_8b0=n_8b[:,1]/1e3 #neutrino fluxes taken from article: arXiv: 1601.00972 (2016)

flux_8b=scipy.interpolate.interp1d(en_ns_8b0,f_ns_8b0,kind='linear',bounds_error = False, fill_value = 0)

#hep NEUTRINOS
n_hep=np.genfromtxt('data/hepenergytab.txt')
en_ns_hep0=n_hep[:,0]*1e3
f_ns_hep0=n_hep[:,1]/1e3

flux_hep=scipy.interpolate.interp1d(en_ns_hep0,f_ns_hep0,kind='linear',bounds_error = False, fill_value = 0)

#pp NEUTRINOS
n_pp=np.genfromtxt('data/ppenergytab.txt')
en_ns_pp0=np.concatenate((n_pp[:,0],n_pp[:,2],n_pp[:,4],n_pp[:,6]))*1e3
f_ns_pp0=np.concatenate((n_pp[:,1],n_pp[:,3],n_pp[:,5],n_pp[:,7]))/1e3

flux_pp=scipy.interpolate.interp1d(en_ns_pp0,f_ns_pp0,kind='linear',bounds_error = False, fill_value = 0)

#CNO CYCLE SOLAR NEUTRINOS===================================================================

#15O NEUTRINOS
n_o15=np.genfromtxt('data/o15energytab.txt')
en_ns_o150=n_o15[:,0]*1e3
f_ns_o150=n_o15[:,1]/1e3

#f17 NEUTRINOS
n_f17=np.genfromtxt('data/f17energytab.txt')
en_ns_f170=n_f17[:,0]*1e3
f_ns_f170=n_f17[:,1]/1e3

#n13 NEUTRINOS
n_n13=np.genfromtxt('data/n13energytab.txt')
en_ns_n130=n_n13[:,0]*1e3
f_ns_n130=n_n13[:,1]/1e3

flux_ns_o15=scipy.interpolate.interp1d(en_ns_o150,f_ns_o150,kind='linear',bounds_error = False, fill_value = 0)
flux_ns_n13=scipy.interpolate.interp1d(en_ns_n130,f_ns_n130,kind='linear',bounds_error = False, fill_value = 0)
flux_ns_f17=scipy.interpolate.interp1d(en_ns_f170,f_ns_f170,kind='linear',bounds_error = False, fill_value = 0)

"""
#CNO fluctuate independantly, had to remove this
int_o15 = scint.quad(flux_ns_o15,0,1e5)[0]
int_n13 = scint.quad(flux_ns_n13,0,1e5)[0]
int_f17 = scint.quad(flux_ns_f17,0,1e5)[0]

int_cno_total = int_o15 + int_f17 + int_n13  #total flux CNO cycle neutrinos

#merge all three sources of CNO neutrinos into 1 spectrum (same uncertainty)
e_cno = np.linspace(0,1e4,10000)
f_o15_cno = flux_ns_o15(e_cno)/int_cno_total
f_n13_cno = flux_ns_n13(e_cno)/int_cno_total
f_f17_cno = flux_ns_f17(e_cno)/int_cno_total


flux_cno = scipy.interpolate.interp1d(e_cno,f_o15_cno+f_n13_cno+f_f17_cno,kind='linear',bounds_error = False, fill_value = 0)
"""

#ATMOSPHERIC NEUTRINOS=====================================================
#data taken from : DOI: 0.1016/j.astropartphys.2005.03.006 (Battistoni, 2005)  !!!data for Grand Sasso, not Snolab!!!
atm_e=np.genfromtxt('data/e_atm_neutrino.txt') #data for atmospheric electron neutrinos , energy in GeV, flux in neutrino/GeV/me/s
en_na_e=atm_e[:,0]*1e6
f_na_e=atm_e[:,1]/1e6 #neutrinos
f_na_antie=atm_e[:,4]/1e6 #anti-neutrinos

atm_mu=np.genfromtxt('data/mu_atm_neutrino.txt') #data for atm. muon neutrinos
en_na_mu=atm_mu[:,0]*1e6
f_na_mu=atm_mu[:,1]/1e6 #neutrinos
f_na_antimu=atm_mu[:,4]/1e6 #anti-neutrinos

#compute total flux of atmospheric neutrinos

#interpolations from data
flux_na_e=scipy.interpolate.interp1d(en_na_e,f_na_e,kind='linear',bounds_error = False, fill_value = 0) 

flux_na_antie=scipy.interpolate.interp1d(en_na_e,f_na_antie,kind='linear',bounds_error = False, fill_value = 0)

flux_na_mu=scipy.interpolate.interp1d(en_na_mu,f_na_mu,kind='linear',bounds_error = False, fill_value = 0)

flux_na_antimu=scipy.interpolate.interp1d(en_na_mu,f_na_antimu,kind='linear',bounds_error = False, fill_value = 0)

int_atm_e = scint.quad(flux_na_e,0,1e6)[0]
int_atm_antie = scint.quad(flux_na_antie,0,1e6)[0]
int_atm_mu = scint.quad(flux_na_mu,0,1e6)[0]
int_atm_antimu = scint.quad(flux_na_antimu,0,1e6)[0]

int_atm_total = int_atm_e + int_atm_antie + int_atm_mu + int_atm_antimu #total flux atmospheric neutrinos

e_atm = np.linspace(1e3,1e6,10000)
f_e_atm = flux_na_e(e_atm)/int_atm_total
f_antie_atm = flux_na_antie(e_atm)/int_atm_total
f_mu_atm = flux_na_mu(e_atm)/int_atm_total
f_antimu_atm = flux_na_antimu(e_atm)/int_atm_total


flux_atm_e = scipy.interpolate.interp1d(e_atm,f_e_atm+f_antie_atm,kind='linear',bounds_error = False, fill_value = 0)
flux_atm_mu = scipy.interpolate.interp1d(e_atm,f_mu_atm+f_antimu_atm,kind='linear',bounds_error = False, fill_value = 0)

# DSNB ==============================================================================
#data taken from: arXiv: 0812.3157 (2009) Fermi-Dirac distribution integrated over expected spectrum of expected supernovas at different 
# z with T=3 for electron neutrinos, T=5 for anti-electron neutrinos and T=8 for other neutrinos

dsnb_nu=np.genfromtxt('data/dsnb.txt') #energy in MeV , flux in neutrino/cm2/MeV/s

#data for electron neutrinos
en_ndsnb_e=dsnb_nu[:,0]*1e3
f_ndsnb_e=dsnb_nu[:,1]/1e3
flux_ndsnb_e=scipy.interpolate.interp1d(en_ndsnb_e,f_ndsnb_e,kind='linear',bounds_error = False, fill_value = 0)

#data for DSNB electron anti-neutrinos
f_ndsnb_antie=dsnb_nu[:,2]/1e3
flux_ndsnb_antie=scipy.interpolate.interp1d(en_ndsnb_e,f_ndsnb_antie,kind='linear',bounds_error = False, fill_value = 0)

#data for other neutrinos and anti-neutrinos (muon and tau)
f_ndsnb_other=dsnb_nu[:,3]/1e3
flux_ndsnb_other=scipy.interpolate.interp1d(en_ndsnb_e,f_ndsnb_other,kind='linear',bounds_error = False, fill_value = 0)


int_dsnb_e = scint.quad(flux_ndsnb_e,0,1e6)[0]
int_dsnb_antie = scint.quad(flux_ndsnb_antie,0,1e6)[0]
int_dsnb_other = scint.quad(flux_ndsnb_other,0,1e6)[0]

int_dsnb_total = int_dsnb_e + int_dsnb_antie +int_dsnb_other #total flux DSNB electron neutrinos

e_dsnb = np.linspace(1e2,1e5,10000)
f_e_dsnb = flux_ndsnb_e(e_dsnb)/int_dsnb_total
f_antie_dsnb = flux_ndsnb_antie(e_dsnb)/int_dsnb_total
f_other_dsnb = flux_ndsnb_other(e_dsnb)/int_dsnb_total


flux_dsnb_e = scipy.interpolate.interp1d(e_dsnb,f_e_dsnb+f_antie_dsnb,kind='linear',bounds_error = False, fill_value = 0)
flux_dsnb_other = scipy.interpolate.interp1d(e_dsnb,f_other_dsnb,kind='linear',bounds_error = False, fill_value = 0)

#uncomment to check validity of spectrum obtained

#make sure all fluxes properly normalized to 1

emin = 0
emax = 1e6

norm_8b = scint.quad(flux_8b,emin,emax)[0]
norm_hep = scint.quad(flux_hep,emin,emax)[0]
norm_pp = scint.quad(flux_pp,emin,1e3)[0]
norm_15o = scint.quad(flux_ns_o15,emin,2e4)[0]
norm_13n = scint.quad(flux_ns_n13,emin,2e4)[0]
norm_17f = scint.quad(flux_ns_f17,emin,2e4)[0]
norm_atm = scint.quad(flux_atm_e,emin,emax)[0] + scint.quad(flux_atm_mu,emin,emax)[0]
norm_dsnb = scint.quad(flux_dsnb_e,emin,emax)[0] + scint.quad(flux_dsnb_other,emin,emax)[0]

print(norm_8b,norm_hep,norm_pp,norm_15o,norm_13n,norm_17f,norm_atm,norm_dsnb)

#show all neutrino fluxes on a plot

e_neutrino = np.logspace(1,6,1000)
dist_f_dsnb = flux_dsnb_e(e_neutrino)*dsnb_flux + flux_dsnb_other(e_neutrino)*dsnb_flux
dist_f_dsnb = np.ma.masked_where(dist_f_dsnb<1e-3,dist_f_dsnb)

dist_f_atm = flux_atm_e(e_neutrino)*atm_flux + flux_atm_mu(e_neutrino)*atm_flux
dist_f_atm = np.ma.masked_where(dist_f_atm<1e-3,dist_f_atm)

dist_f_pp = flux_pp(e_neutrino)*pp_flux
dist_f_8b = flux_8b(e_neutrino)*b8_flux
dist_f_hep = flux_hep(e_neutrino)*hep_flux
dist_f_o15 = flux_ns_o15(e_neutrino)*o15_flux
dist_f_n13 = flux_ns_n13(e_neutrino)*n13_flux
dist_f_f17 = flux_ns_f17(e_neutrino)*f17_flux
#positions for uncertainty labels on plot
def pos(dist,extra=0):
	ind = np.argmax(dist)
	y = dist[ind]
	x = e_neutrino[ind]
	if extra == 1:
		x=x*1.3
	return 1.5*x,y
color_tab = ['black','red','blue','darkorange','darkgreen','limegreen','purple','magenta','saddlebrown','darkgray']
plt.figure(figsize=(15,9))

plt.plot(e_neutrino,dist_f_pp,color = color_tab[0],linewidth=4,label='pp')
plt.plot(e_neutrino,dist_f_8b,color = color_tab[1],linewidth=4,label=r'$^{8}B$')
plt.plot(e_neutrino,dist_f_hep,color = color_tab[2],linewidth=4,label='hep')
plt.plot(e_neutrino,dist_f_o15,color = color_tab[3],linewidth=4,label=r'$^{15}O$')
plt.plot(e_neutrino,dist_f_n13,color = color_tab[4],linewidth=4,label=r'$^{13}N$')
plt.plot(e_neutrino,dist_f_f17,color = color_tab[5],linewidth=4,label=r'$^{17}F$')
plt.plot([e_7be_2,e_7be_2],[0,be7_2_flux],color = color_tab[8],linewidth=4,label = r'$^{7}Be$')
plt.plot([e_7be_1,e_7be_1],[0,be7_1_flux],color = color_tab[8],linewidth=4)
plt.plot([e_pep,e_pep],[0,pep_flux],color = color_tab[9],linewidth=4,label = r'$pep$')
plt.plot(e_neutrino,dist_f_atm,color = color_tab[6],linewidth=4,label='atmospheric')
plt.plot(e_neutrino,dist_f_dsnb,color = color_tab[7],linewidth=4,label='DSNB')


plt.text(*pos(dist_f_pp),r'$\pm 1\%$',color = color_tab[0],fontsize=18)
plt.text(*pos(dist_f_8b),r'$\pm 16\%$',color = color_tab[1],fontsize=18)
plt.text(*pos(dist_f_hep,extra=1),r'$\pm 16\%$',color = color_tab[2],fontsize=18)
plt.text(*pos(dist_f_o15,extra=1),r'$\pm 5\%$',color = color_tab[3],fontsize=18)
plt.text(*pos(dist_f_n13),r'$\pm 5\%$',color = color_tab[4],fontsize=18)
plt.text(*pos(dist_f_f17,extra=1),r'$\pm 5\%$',color = color_tab[5],fontsize=18)
plt.text(e_7be_2,3*be7_2_flux,r'$\pm 10.5\%$',color = color_tab[8],fontsize=18)
plt.text(e_7be_1,3*be7_1_flux,r'$\pm 10.5\%$',color = color_tab[8],fontsize=18)
plt.text(e_pep,3*pep_flux,r'$\pm 2\%$',color = color_tab[9],fontsize=18)
plt.text(*pos(dist_f_atm,extra=1),r'$\pm 20\%$',color = color_tab[6],fontsize=18)
plt.text(*pos(dist_f_dsnb),r'$\pm 50\%$',color = color_tab[7],fontsize=18)

plt.xlabel('Neutrino energy (keV)',fontsize=24)
plt.ylabel(r'flux ($s^{-1} m^{-2} keV^{-1}$)',fontsize=24)
plt.xscale('log')
plt.yscale('log')
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlim([20,9*10**5])
plt.ylim([10**(-3),10**15])
plt.legend(loc='best',fontsize=20)
plt.savefig('fluxes_neu.png')
plt.show()


#PRODUCE A NUMPY ARRAY WITH "NORMALIZED" RECOIL DISTRIBUTIONS FOR ALL 10 SOURCES OF NEUTRINO, FOR NR AND ER

def numpy_energy(emax,num,AA,iso,ZZ,name):
	#creates a .npy file  with differential recoil dist from 0 to emax with num points for element of mass number AA
	#No detector properties included here (acceptance, ER rejection, size, etc...)
	e_tab = np.logspace(-4,2,1000)
	dist = np.zeros([21,1000])
	dist[0,:] = e_tab



	for j,e in enumerate(e_tab):
		#NR
		dist[1,j] = eval_neu(e,AA,iso,ZZ,flux_pp,emax=1e3)
		dist[2,j] = eval_neu(e,AA,iso,ZZ,flux_hep,emax=1e5)
		dist[3,j] = eval_neu(e,AA,iso,ZZ,flux_8b,emax=1e5)
		dist[4,j] = eval_neu(e,AA,iso,ZZ,flux_ns_o15,emax=1e4)
		dist[5,j] = eval_neu(e,AA,iso,ZZ,flux_ns_n13,emax=1e4)
		dist[6,j] = eval_neu(e,AA,iso,ZZ,flux_ns_f17,emax=1e4)
		dist[7,j] = eval_neu(e,AA,iso,ZZ,flux_pp,spectral_line='pep')
		dist[8,j] = eval_neu(e,AA,iso,ZZ,flux_pp,spectral_line='7be')
		dist[9,j] = eval_neu(e,AA,iso,ZZ,flux_atm_e) + eval_neu(e,AA,iso,ZZ,flux_atm_mu)
		dist[10,j] = eval_neu(e,AA,iso,ZZ,flux_dsnb_e) + eval_neu(e,AA,iso,ZZ,flux_dsnb_other)

		#ER
		dist[11,j] = eval_e(e,AA,iso,ZZ,flux_pp,flux_pp,solar=True,emax=1e3)
		dist[12,j] = eval_e(e,AA,iso,ZZ,flux_hep,flux_hep,solar=True,emax=1e5)
		dist[13,j] = eval_e(e,AA,iso,ZZ,flux_8b,flux_8b,solar=True,emax=1e5)
		dist[14,j] = eval_e(e,AA,iso,ZZ,flux_ns_o15,flux_ns_o15,solar=True,emax=1e4)
		dist[15,j] = eval_e(e,AA,iso,ZZ,flux_ns_n13,flux_ns_n13,solar=True,emax=1e4)
		dist[16,j] = eval_e(e,AA,iso,ZZ,flux_ns_f17,flux_ns_f17,solar=True,emax=1e4)
		dist[17,j] = eval_e(e,AA,iso,ZZ,flux_pp,flux_pp,solar=True, spectral_line='pep')
		dist[18,j] = eval_e(e,AA,iso,ZZ,flux_pp,flux_pp,solar=True,spectral_line='7be')
		dist[19,j] = eval_e(e,AA,iso,ZZ,flux_atm_e,flux_atm_mu)
		dist[20,j] = eval_e(e,AA,iso,ZZ,flux_dsnb_e,flux_dsnb_other)

		print(j)
		#check recoil spectrum

	dist[1:,:] = np.ma.masked_where(dist[1:,:]<1e-17,dist[1:,:])

	# NR distribution
	plt.plot(dist[0,:],dist[1,:]*pp_flux,label='pp')
	plt.plot(dist[0,:],dist[2,:]*hep_flux,label='hep')
	plt.plot(dist[0,:],dist[3,:]*b8_flux,label='8b')
	plt.plot(dist[0,:],dist[4,:]*o15_flux,label='15O')
	plt.plot(dist[0,:],dist[5,:]*n13_flux,label='13N')
	plt.plot(dist[0,:],dist[6,:]*f17_flux,label='17F')
	plt.plot(dist[0,:],dist[7,:]*pep_flux,label='pep')
	plt.plot(dist[0,:],dist[8,:]*be7_flux,label='be7')
	plt.plot(dist[0,:],dist[9,:]*atm_flux,label='atm')
	plt.plot(dist[0,:],dist[10,:]*dsnb_flux,label='DSNB')
	plt.xlabel('recoil energy (keV)')
	plt.ylabel(r'flux ($s^{-1} m^{-2} keV^{-1}$)')
	plt.xscale('log')
	plt.yscale('log')
	plt.legend(loc='best')
	plt.show()


	#ER distribuition
	plt.plot(dist[0,:],dist[11,:]*pp_flux,label='pp')
	plt.plot(dist[0,:],dist[12,:]*hep_flux,label='hep')
	plt.plot(dist[0,:],dist[13,:]*b8_flux,label='8b')
	plt.plot(dist[0,:],dist[14,:]*o15_flux,label='15O')
	plt.plot(dist[0,:],dist[15,:]*n13_flux,label='13N')
	plt.plot(dist[0,:],dist[16,:]*f17_flux,label='17F')
	plt.plot(dist[0,:],dist[17,:]*pep_flux,label='pep')
	plt.plot(dist[0,:],dist[18,:]*be7_flux,label='be7')
	plt.plot(dist[0,:],dist[19,:]*atm_flux,label='atm')
	plt.plot(dist[0,:],dist[20,:]*dsnb_flux,label='DSNB')
	plt.xlabel('recoil energy (keV)')
	plt.ylabel(r'flux ($s^{-1} m^{-2} keV^{-1}$)')
	plt.xscale('log')
	plt.yscale('log')
	plt.legend(loc='best')
	plt.show()

	#np.save(name, dist)

	return dist


dist_argon = numpy_energy(200,1000,A_argon_s,iso_1,Z_argon,'dist_neutrino_argon_long.npy')
dist_xenon = numpy_energy(100,1000,A_xenon_s,iso_1,Z_xenon,'dist_neutrino_xenon_long.npy')