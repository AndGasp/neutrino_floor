# coding=utf-8
#Constants=============================================================================================
import numpy as np

#basic constants========================================================================================
cc=3e8 #m/s speed of light
m_p=0.939e6 #proton mass keV
m_e=510.998 #electron mass in keV
sin2_wma=0.2223 #weak mixing angle
g_f=1.166364e-17 #kev-2 (fermi coupling constant)


#Dark matter distribution constants=====================================================================
dens_dm=0.5e12 #dark matter density (keV/m3)
err_dens=0.2e12 #unncertainty on DM density

dens_dm_old = 0.3e12
err_dens_dm_old = 0.1e12

#SHM usual parameters
v_earth=232e3/cc #speed of solar system around galaxy 

v_0=233e3/cc #central speed of DM distribution 
err_v0=5e3/cc

v_0_old = 220e3/cc
err_v0_old = 20e3/cc

v_esc=528e3/cc#escape velocity of DM (m/s)
err_vesc=25e3/cc

v_esc_old = 544e3/cc
err_vsec_old = 70e3/cc
#inelastic DM
v_ref=200e3/cc

#Detector-related constants==============================================================================

A_argon_s=np.array([39.948]) 
A_xenon_s=np.array([131.293])
A_ger_s=np.array([72.64])

Z_argon=np.array([18])
Z_xenon=np.array([54])
Z_ger=np.array([32])

A_c3f8 = np.array([12.0107,18.9984])
iso_c3f8 = np.array([3/11,8/11])
Z_c3f8 = np.array([6,9])

iso_1=np.array([1]) #when only one element, isotope fraction = 1

A_argon_l=np.array([38,36,40]) #use complete version on A arrays only where distinguishing individual isotopes is important (isospin violating DM)
iso_argon=np.array([0.063,0.334,99.604])/100 #isotope distribution in %
A_xenon_l=np.array([126,128,129,130,131,132,134,136])
iso_xenon=np.array([0.089,1.191,26.401,4.071,21.232,26.909,10.436,8.857])/100
A_ger_l=np.array([70,72,73,74,76])
iso_ger=np.array([20.52,27.45,7.76,36.52,7.75])/100


#Actual experimental thresholds

#DEAP-3600
E_min_deap=55#minimum energy ROI zone keV
E_max_deap=111#max energy ROI zone
rej_e_deap=1e9 #electronic recoil rejection efficiency
lim_argon=15 #energy (keV) below which solar neutrinos dominate

#LUX
E_min_lux=1.1
E_max_lux=80
rej_e_xenon=1/(1-0.999)
lim_xenon=6 #energy (keV) below which solar neutrinos dominate 

#PandaX-II
E_min_panda=1
E_max_panda=81

#Xenon1T
E_min_xenon=4.9
E_max_xenon=40.9

#SuperCDMS
E_min_cdms=10
E_max_cdms=100


#actual fiducial masses of experiments, taken from same papers as acceptances
M_deap_new=247*1000*24*3600 #more recent, unpublished data
M_deap=14226*24*3600 #2017 paper
M_lux=3.35e4*24*3600
M_panda=5.4e4*24*3600 #more recent arXiv: 1708.06917 (2017) paper
M_xenon=1300*279*24*3600
M_cdms=220*24*3600


#IDEAL PARAMETERS FOR FUTURE EXPERIMENT
E_min_Ar = 15
E_max_Ar = 100
rej_er_Ar = 1e8

E_min_Xe = 1
E_max_Xe = 90
rej_er_Xe = 1e3

#total background measured in different experiments (supposing that predicted background is the same, since no experiment observed significant #events above predicted background)
back_deap=1e-30 #zero
back_lux=591
back_panda=74
back_xenon=7.36
back_cdms=1.43


#since experiments used Monte-Carlo fits to compute likelihoods, lux and pandax did not publish number of expected events in
#ROI, but know that expected ER background for Xenon1T 0.2mDRU, for 0.8 mDRU for PANDAX and 2.6 mDRU for LUX (mDRU = 10e-3 events /kg /day/KeV)
#XENON1T published that they expected 
back_approx_lux=7.36/M_xenon*M_lux*2.6/0.2
back_approx_panda=7.36/M_xenon*M_panda*0.8/0.2

#Neutrino-related constants==============================================================================

#error on different neutrino fluxes from arXiv: astro-ph/0412096 (Bahcall, 2005) for solar neutrinos and arXiv: 0903.3630 (Strigari, 2009) for others
err_1=16 #8B
err_2=16 #hep
err_3=1 #pp
err_4=5 #CNO
err_5=5 
err_6=5
err_7=20 #DSNB
err_8=20
err_9=50 #atm
err_10=50
err_11=2 #pep
err_12=10.5 # 2 7Be lines
err_13=10.5

#fluxes of different sources of neutrinos
pp_flux = 5.98e14 #neutrinos /m^2/s
n13_flux = 2.96e12
o15_flux = 2.23e12
f17_flux = 5.52e10
b8_flux = 5.58e10
hep_flux = 8.04e7

cno_flux = o15_flux+n13_flux+f17_flux

be7_1_flux = 4.84e12
be7_2_flux = 4.35e13
be7_flux = be7_1_flux + be7_2_flux
pep_flux = 1.44e12

atm_flux = 10.5e4
dsnb_flux = 85.5e4


#energies of Be and pep neutrino fluxes
e_7be_2=861.3 #keV
e_7be_1=384.3
e_pep=1.44*1e3


err_nu=np.array([16,16,1,5,5,5,20,20,50,50])/100 #array with all errors on neutrino flux values


av_values = np.array([v_esc,v_0,dens_dm,pp_flux,hep_flux,b8_flux,o15_flux,n13_flux,f17_flux,pep_flux,be7_flux,atm_flux,dsnb_flux])
sig_values = np.array([err_vesc,err_v0,err_dens,0.01*pp_flux,0.16*hep_flux,0.18*b8_flux,0.05*o15_flux,0.05*n13_flux,0.05*f17_flux,0.02*pep_flux,0.105*be7_flux,0.2*atm_flux,0.5*dsnb_flux])


#electron neutrino survival rate
surv=0.57
err_surv=0.03
