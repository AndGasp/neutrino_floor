# coding=utf-8
#neutrino coherent scattering and electron scattering relevant functions====================================================================
import numpy as np
import scipy.integrate as scint
from scipy.interpolate import interp1d
from constants import *

def nuc_F(e,A):
	#Helm nuclear form factor
	#this form taken from DOI: 10.1142/S0218301392000023 (Hengel and al. 1992), values of constants taken from the Lewin paper (1996)
	s=0.9
	a=0.52
	c=(1.23*A**(1/3)-0.6)
	r_n=(c**2+(7*np.pi**2*a**2)/3-5*s**2)**(0.5)
	q_rn = 6.92e-3*A**(0.5)*(e)**(0.5)*r_n
	q_s = 6.92e-3*A**(0.5)*(e)**(0.5)*s

	f=3*(np.sin(q_rn)-q_rn*np.cos(q_rn))*np.exp((-(q_s)**2)/2) / ((q_rn)**3)
	
	return f

def int_n2_sl(e_r,AA,e_neu,f_neu):
	emin=(AA*m_p*e_r/2)**(0.5)
	con=1/cc**2*e_r*nuc_F(e_r,AA)**2
	rr=0
	for i in range(len(e_neu)):
		if emin<e_neu[i]:
			rr+=f_neu[i]/e_neu[i]**2

	return rr*con

def int_n1_sl(e_r,AA,e_neu,f_neu):
	emin=(AA*m_p*e_r/2)**(0.5)
	con=1/cc**2*nuc_F(e_r,AA)**2
	rr=0
	for i in range(len(e_neu)):
		if emin<e_neu[i]:
			rr+=f_neu[i]
	return	rr*con

def int_n4(e_nu,flux):
	return flux(e_nu)/(e_nu**2)

def int_n3(e_nu,flux):
	return flux(e_nu)

def int_n2(e_r,AA,flux,emax):
	emin=(AA*m_p*e_r/2)**(0.5)
	con=1/cc**2*e_r*nuc_F(e_r,AA)**2
	if emin>emax:
		rr=0
	else:
		rr=scint.quad(int_n4,emin,emax,args=(flux,))[0]
	
	return rr*con

def int_n1(e_r,AA,flux,emax):
	emin=(AA*m_p*e_r/2)**(0.5)
	con=1/cc**2*nuc_F(e_r,AA)**2

	if emin>emax:
		rr=0
	else:
		rr=scint.quad(int_n3,emin,emax,args=(flux,))[0]

	return	rr*con

def eval_neu(e,AA,iso,ZZ,flux,spectral_line='no',emax=1e6): 
	#evaluate differential event rate at specific energy
	#AA table with mass numbers, iso table with isotope fractions, ZZ number of electrons, aaaccc acceptance functions, flux neutrino flux
	R_n=0
	for i in range(len(AA)): #if different mass numbers considered (if isotopes important, i.e. isospin violating DM)
		conss=0.197e-9 #convert unit from keV-2 to m^2 for cross section
		NN=round(AA[i]-ZZ[i]) #number of neutrons in nuclei
		Q=NN-(1-4*sin2_wma)*ZZ[i] #weak charge
		c_n=conss**2*AA[i]*m_p*g_f**2*Q**2/(4*np.pi)*(6.022e23/(AA[i]*1e-3))
		if 'no' in spectral_line:
			R_n+=iso[i]*c_n*cc**2*(int_n1(e,AA,flux,emax)-m_p*AA[i]/2*int_n2(e,AA,flux,emax))

		else:
			if 'pep' in spectral_line:
				e_neu = np.array([e_pep])
				f_neu = np.array([1])
			if '7be' in spectral_line:
				e_neu = np.array([e_7be_1,e_7be_2])
				f_neu = np.array([be7_1_flux/be7_flux,be7_2_flux/be7_flux])

			R_n += iso[i]*c_n*cc**2*(int_n1_sl(e,AA,e_neu,f_neu)-m_p*AA[i]/2*int_n2_sl(e,AA,e_neu,f_neu))

	return R_n


# neutrino electron scattering===========================================================================================================

def int_e1_sl(e_r,AA,typ,e_neu,f_neu):
	#special function to compute e_r distribution if neutrino flux is in the shape of a single spectral line
	emin=0.5*(e_r+(e_r*(e_r+2*m_e))**0.5)

	g_a=0.5
	g_v=2*0.23-0.5
	if 'e' in typ:
		g_a+=1
		g_v+=1

	con=((g_v+g_a)**2+(g_v-g_a)**2)/cc**2

	rr=0

	for i in range(len(e_neu)):
		if emin<e_neu[i]:
			rr+=f_neu[i]

	return	rr*con


def int_e2_sl(e_r,AA,typ,e_neu,f_neu):

	emin=0.5*(e_r+(e_r*(e_r+2*m_e))**0.5)
	g_a=0.5
	g_v=2*0.23-0.5
	if 'e' in typ:
		g_a+=1
		g_v+=1

	con=e_r/cc**2

	rr1=0
	rr2=0
	for i in range(len(e_neu)):

		if emin<e_neu[i]:
			rr1+=f_neu[i]/e_neu[i]
			rr2+=f_neu[i]/e_neu[i]**2

	return con*(-2*(g_v-g_a)**2*rr1+m_e*(g_a**2-g_v**2)*rr2)

def int_e3_sl(e_r,AA,typ,e_neu,f_neu):
	emin=0.5*(e_r+(e_r*(e_r+2*m_e))**0.5)
	g_a=0.5
	g_v=2*0.23-0.5
	if 'e' in typ:
		g_a+=1
		g_v+=1

	con=1/cc**2*e_r**2*(g_v-g_a)**2
	rr=0

	for i in range(len(e_neu)):
		if emin<e_neu[i]:
			rr+=f_neu[i]/e_neu[i]**2

	return con*rr

def int_e6(e_nu,flux):
	return flux(e_nu)/(e_nu**2)

def int_e5(e_nu,flux):
	return flux(e_nu)/(e_nu)

def int_e4(e_nu,flux):
	return flux(e_nu)


def int_e3(e_r,AA,flux,typ,solar,emax):
	emin=0.5*(e_r+(e_r*(e_r+2*m_e))**0.5)
	g_a=0.5
	g_v=2*0.23-0.5
	if 'e' in typ:
		g_a+=1
		g_v+=1
	
	frac = 1

	if solar == True:
		if 'e' in typ:
			frac = surv
		else:
			frac = 1-surv

	con=1/cc**2*e_r**2*(g_v-g_a)**2

	if emin<emax:
		rr=0
	else:
		rr=scint.quad(int_e6,emin,emax,args=(flux,))[0]

	return frac * con*rr

def int_e2(e_r,AA,flux,typ,solar,emax):

	emin=0.5*(e_r+(e_r*(e_r+2*m_e))**0.5)
	g_a=0.5
	g_v=2*0.23-0.5
	if 'e' in typ:
		g_a+=1
		g_v+=1

	frac = 1

	if solar == True:
		if 'e' in typ:
			frac = surv
		else:
			frac = 1-surv

	con=e_r/cc**2

	if emin<emax:
		rr1=0
		rr2=0
	else:
		rr1=scint.quad(int_e5,emin,emax,args=(flux,))[0]
		rr2=scint.quad(int_e6,emin,emax,args=(flux,))[0]

	return frac * con*(-2*(g_v-g_a)**2*rr1+m_e*(g_a**2-g_v**2)*rr2)

def int_e1(e_r,AA,flux,typ,solar,emax):

	emin=0.5*(e_r+(e_r*(e_r+2*m_e))**0.5)

	g_a=0.5
	g_v=2*0.23-0.5

	if 'e' in typ:
		g_a+=1
		g_v+=1

	frac = 1

	if solar == True:
		if 'e' in typ:
			frac = surv
		else:
			frac = 1-surv

	con=((g_v+g_a)**2+(g_v-g_a)**2)/cc**2
	if emin>emax:
		rr=0
	else:
		rr=scint.quad(int_e4,emin,emax,args=(flux,))[0]

	return	frac * rr*con

#list of ionization energies of Argon (in keV)
ion_argon = np.array([ 15.75990912,  27.62907124,   40.74194577,   59.81220276,   75.016587,
   91.00865577,  124.31941988,  143.46222676,  422.44764939,  478.68417896,
  538.96277388,  618.25980444,  686.10431653,  755.74184735,  854.77244818,
  918.02558862, 4120.88561418, 4426.22737567])/1000

ion_xenon = np.array([0]) #not found 


def z_eff(e,Z):
	#effective Z for argon for a given neutrino kinetic energy for step approximation
	if Z == 18:
		arg = np.sum(ion_argon > e)
		z_eff = 18 - arg
	if Z == 54:
		arg = np.sum(ion_xenon > e)
		z_eff = 54 - arg
	return z_eff


#fit data for xenon electron-solar neutrino cross-section
data_cc_xe = np.load('data/er_cc_xe.npy') #array containing extrapolated data for cross section of neutrino-electron int.
ee = data_cc_xe[0,:]
cc_pp = data_cc_xe[1,:] #measured cross section for pp neutrinos
cc_be = data_cc_xe[2,:] #measured corss section for be7 neutrinos (862 and 384 keV lines combined)
#note that for all values above 30 keV, a polynomial fit of order 2 on the tail of the data was used to extrapolate the results

fit_pp = interp1d(ee,cc_pp)
fit_be = interp1d(ee,cc_be)

def eval_e(e,AA,iso,ZZ,fluxe,fluxo,approx,solar=False,spectral_line='no',emax=1e6):
	#evaluate event rate at specific energy point 
	#AA atomic mass (or array of atomic masses of isotopes), iso fractions of each isotopes, ZZ electron number, fluxe elctron neutrino flux, fluxo other neutrinos
	# If solar neutrinos, enter the complete flux for fluxe and fluxo and solar=True so the proper fractions of electron neutrinos and other neutrinos will be considered with respect to electron neutrino survival rate
	# New additon: approx = method to approximate electrons-neutrino scatt used ('SA' for step approx - Argon, or 'data' for exp. results from paper - Xenon)
	R_n_o=0
	R_n_e=0

	for i in range(len(AA)):
		const=g_f**2*m_e/(2*np.pi)*(6.022e23/(AA[i]*1e-3))
		conss=0.197e-9
		if 'no' in spectral_line: 
		#pp
			if approx == 'SA':
				R_n_o+=iso[i]*z_eff(e,ZZ[0])*const*conss**2*cc**2*(int_e1(e,AA[i],fluxo,'o',solar,emax)+int_e2(e,AA[i],fluxo,'o',solar,emax)+int_e3(e,AA[i],fluxo,'o',solar,emax))
				R_n_e+=iso[i]*z_eff(e,ZZ[0])*const*conss**2*cc**2*(int_e1(e,AA[i],fluxe,'e',solar,emax)+int_e2(e,AA[i],fluxe,'e',solar,emax)+int_e3(e,AA[i],fluxe,'e',solar,emax))
			if approx == 'data':
				R_n_e += fit_pp(e)



		else:
			if 'pep' in spectral_line:
				e_neu = np.array([e_pep])
				f_neu = np.array([1])
			if '7be' in spectral_line:
				e_neu = np.array([e_7be_1,e_7be_2])
				f_neu = np.array([be7_1_flux/be7_flux,be7_2_flux/be7_flux])
				if approx == 'SA':
					R_n_o +=iso[i]*z_eff(e,ZZ[0])*const*conss**2*cc**2*(int_e1_sl(e,AA[i],'o',e_neu,f_neu*(1-surv))+int_e2_sl(e,AA[i],'o',e_neu,f_neu*(1-surv))+int_e3_sl(e,AA[i],'o',e_neu,f_neu*(1-surv)))
					R_n_e +=iso[i]*z_eff(e,ZZ[0])*const*conss**2*cc**2*(int_e1_sl(e,AA[i],'e',e_neu,f_neu*surv)+int_e2_sl(e,AA[i],'e',e_neu,f_neu*surv)+int_e3_sl(e,AA[i],'e',e_neu,f_neu*surv))
				if approx == 'data':
					R_n_e += fit_be(e)

	return R_n_o+R_n_e

