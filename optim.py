# coding=utf-8
# Andrea Gaspert 02/2019
# All functions to perform optimisation and monte carlo and compute binned profile likelihood

import sys,os
import numpy as np
import matplotlib.pyplot as plt
import random
import scipy.interpolate
import scipy.integrate as scint
from neutrino_events import *
from constants import *
from acceptance import acc_100
from dm_events import *
import time

start_time = time.clock()



#upload numpy arrays with neutrino events
neu_xenon = np.load('dist_neutrino_xenon_long.npy')
neu_argon = np.load('dist_neutrino_argon_long.npy')

order_fit = np.array([1e14,1e7,1e10,1e13,1e4,1e5])
bnds = ((0.01, 100), (0.01, 100), (0.01, 100), (0.01, 100), (0.01, 100), (0.01, 100)) #limits on floating of fluxes within +/- 2 orders of mignitudes of av. expected value


def isNaN(num): 
	#return true if NaN
    return num != num

def empty(v):
	return v

def int_simpson(e,dist):
	#takes in array with speeds and array with distribution to integrate and returns value of integral using simpson method
	#v need to be equally spaced and ordered!
	if len(e)/2 == len(e)//2: #simpson method requires even number of subdivisions
		e = e[:-1]
		dist = dist[:-1]

	h = e[1] - e[0] #spacing

	sum_simpson = h/3 * (dist[0] + 4*np.sum(dist[1::2]) + 2*np.sum(dist[:-1][2::2]) + dist[-1])

	return sum_simpson


def dm_bins(emin,emax,nbin,m,acc,AA,ZZ,dm_params,speed_dist,dist):
	#Compute number of dark matter events in each bin per kg per second per m^2 (divided by cross-section) for given detector properties
	#emin and emax ROI limits, m dark matter mass, acc acceptance as a function of NR energy, AA mass number of detector medium, ZZ proton number
	#n bin number of energy bins
	#dm_params array containing parameters of the velocity and density distribution of DM (escape velocity, velocity of the sun, density)
	#if speed_dist 'reg' , use analytical computation of dark matter velocity integral (int_new), 'pp' use int_irreg and profile over beta and eta values of SHM++ distribution
	#else use int_other_dist and compute integral by hand 
	#dist to specify distribution if other distribution is used (has to be defined in dm_events)

	vesc = dm_params[0]
	v0 = dm_params[1]
	beta_ind = dm_params[2]
	eta_ind = dm_params[3]

	e_tab = np.logspace(np.log10(emin),np.log10(emax),nbin+1) #definition of energy bins

	dm_events = np.zeros(nbin)
	dm_events_new = np.zeros(nbin)

	if 'reg' in speed_dist: #if want to use SHM

		ee = np.linspace(emin,emax,100000)

		RR = int_new(ee,m,AA,ZZ,acc,iso_1,vesc,v0,1,1,'none',0,0)

		for i in range(nbin-1):
			"""
			#old integration method, time consuming
			start1 = time.time()

			dm_events[i] = scint.quad(int_new, e_tab[i] , e_tab[i+1], args=(m,AA,ZZ,acc,iso_1,vesc,v0,1,1,'none',0,0))[0]
			
			time1 = time.time() - start1
			"""

			#start2 = time.time()

			ind_min = np.argmin(np.abs(ee-e_tab[i]))
			ind_max = np.argmin(np.abs(ee-e_tab[i+1])) 					
	
			dm_events_new[i] = int_simpson(ee[ind_min:ind_max],RR[ind_min:ind_max])

			#time2 = time.time() - start2

			#print('comparison for integrals for DM: {} / {}  vs {} / {}'.format(dm_events[i],time1,dm_events_new[i],time2)) 

	elif 'pp' in speed_dist: #if want to use SHM++

		ee = np.linspace(emin,emax,100000)
		dR = np.zeros(100000)

		#start1 = time.time()
		for k in range(100000):
			dR[k] = int_irreg(ee[k],m,AA,acc,beta_ind,eta_ind) #evaluate dR/dE in SHM++ model

		for i in range(nbin):
			ind_min = np.argmin(np.abs(ee-e_tab[i])) #where speed closest to emin
			ind_max = np.argmin(np.abs(ee-e_tab[i+1])) 	#where speed closest to emax
			dm_events_new[i] = int_simpson(ee[ind_min:ind_max],dR[ind_min:ind_max]) #compute integral using Simpsons approx.
		#time1 = time.time() - start1
		"""
		#old time-consuming method to compute integral
		start2 = time.time()
		for i in range(nbin):

			#print(i)

			dm_events[i] = scint.quad(int_irreg, e_tab[i] , e_tab[i+1], args=(m,AA,acc,beta_ind,eta_ind))[0]
		time2 = time.time() - start2
		"""
		"""
		tabb = np.linspace(0,100,nbin)
		plt.plot(tabb,dm_events)
		plt.yscale('log')
		plt.show()
		"""
		#print('old time {} vs new time {}'.format(time2,time1))
		#print('old values: {}'.format(dm_events))
		#print('new values {}'.format(dm_events_new#))
	else:

		for i in range(nbin): #if other distribution used (dist need to be a function computin speed distribution in lab frame f(v)/v in km/s)

			dm_events[i] = scint.quad(int_other_dist , e_tab[i] , e_tab[i+1], args=(m,AA,ZZ,dist,acc,iso_1,vesc,v0,1,1,'none',0,0))[0]

	return dm_events_new


def neutrino_bin(emin,emax,nbin,neu_array,er_rej,neu_fluxes):
	#Compute number of neutrino events in each bin per kg per second for given detector properties (NR and ER)
	#emin and emax ROI limits, nbin number of energy bins
	#neu_array 17xn array containing normalized distribution of neutrino events as a function of energy for all solar, atm and DSNB neutrinos, NR and ER
	#er_rej electron recoil rejection efficiency of the detector, supposed constant vs energy
	#neu_fluxes array containg neutrino fluxes for all 8 different neutrino sources (pp,hep,8b,cno(15O,17F,13N),pep,7Be,atm.,DSNB)


	e_neu = neu_array[0,:]

	new_array = np.zeros_like(neu_array)


	if type(neu_fluxes) != np.ndarray:
		neu_fluxes = np.array(neu_fluxes)
	#multiple normalized event distribution by appropriate flux
	new_array[1,:]=neu_array[1,:]*neu_fluxes[0] #pp NR
	new_array[11,:]=neu_array[11,:]*neu_fluxes[0] #ER

	new_array[2,:]=neu_array[2,:]*neu_fluxes[1] #hep
	new_array[12,:]=neu_array[12,:]*neu_fluxes[1]

	new_array[3,:]=neu_array[3,:]*neu_fluxes[2] #8b
	new_array[13,:]=neu_array[13,:]*neu_fluxes[2]

	new_array[4,:]=neu_array[4,:]*neu_fluxes[3] #15O
	new_array[14,:]=neu_array[14,:]*neu_fluxes[3]

	new_array[5,:]=neu_array[5,:]*neu_fluxes[4] #13N
	new_array[15,:]=neu_array[15,:]*neu_fluxes[4]

	new_array[6,:]=neu_array[6,:]*neu_fluxes[5] #17F
	new_array[16,:]=neu_array[16,:]*neu_fluxes[5]

	new_array[7,:]=neu_array[7,:]*neu_fluxes[6] #pep
	new_array[17,:]=neu_array[17,:]*neu_fluxes[6]

	new_array[8,:]=neu_array[8,:]*neu_fluxes[7] #7Be
	new_array[18,:]=neu_array[18,:]*neu_fluxes[7]

	new_array[9,:]=neu_array[9,:]*neu_fluxes[8] #atm.
	new_array[19,:]=neu_array[19,:]*neu_fluxes[8]

	new_array[10,:]=neu_array[10,:]*neu_fluxes[9] #DSNB
	new_array[20,:]=neu_array[20,:]*neu_fluxes[9]


	new_array = np.delete(new_array,0,0) #remove column with energy values

	neu_nr, neu_er = np.split(new_array,2,0) #split NR and ER events

	dist_nr = np.sum(neu_nr,axis=0) #sum all different neutrino contributions

	dist_er = np.sum(neu_er, axis=0)/er_rej #divide ER events by ER rejection efficiency

	fit_neu_nr=scipy.interpolate.interp1d(e_neu,dist_nr,kind='linear',bounds_error = False, fill_value = 0) #interpolation fit in all neutrino NR events
	fit_neu_er=scipy.interpolate.interp1d(e_neu,dist_er,kind='linear',bounds_error = False, fill_value = 0) #interpolation fit in all neutrino ER events

	ee = np.linspace(emin,emax,100000)
	RR = fit_neu_nr(ee)+fit_neu_er(ee)

	e_tab = np.logspace(np.log10(emin),np.log10(emax),nbin+1) 

	r_neu = np.zeros(nbin)
	r_neu_new = np.zeros(nbin)

	for i in range(nbin):
		"""
		#old time consuming method

		start1 = time.time()

		r_neu[i] = scint.quad(fit_neu_nr,e_tab[i],e_tab[i+1])[0] + scint.quad(fit_neu_er,e_tab[i],e_tab[i+1])[0]

		time1 = time.time()- start1
		"""
		#start2 = time.time()

		ind_min = np.argmin(np.abs(ee-e_tab[i]))
		ind_max = np.argmin(np.abs(ee-e_tab[i+1]))
			
		r_neu_new[i] = int_simpson(ee[ind_min:ind_max],RR[ind_min:ind_max])

		#time2 = time.time()-start2
		#print('comparison between integral method neutrino:{} / {} vs {} / {}'.format(r_neu[i],time1,r_neu_new[i],time2))

		#print(i)

	return r_neu_new


def likelihood(tab_dm,tab_neu_0,tab_neu_mc,exp,params_0,neu_fluxes_mc,details=0):
	# return z value of DM discovery for given simulated DM and neutrino event distributions
	# tab_dm simulated DM events in each energy bin
	# tab_neu_0 simulated neutrino events in each bin
	# tab_neu_mc neutrino events in each bin that maximize the background only hypothesis
	# exp exposure in kg-seconds
	# params_0 parameters of the neutrino fluxes and DM dist. for simulated events, neu_fluxes_mc neutrino fluxes that maximize background-only likelihood
	# details=1 to have detailed figures at every iteration (for debugging only)


	deno_sig=scipy.special.gammaln(tab_dm*exp+tab_neu_0*exp+1)
	
	p_1_log=-(tab_dm*exp+tab_neu_0*exp)-deno_sig+(exp*tab_neu_0+tab_dm*exp)*np.log(tab_neu_0*exp+tab_dm*exp)
	p_1_log_corrected = np.nan_to_num(p_1_log) #poisson likelihood of events being caused by background + DM signal

	p_0_log=-(tab_neu_mc*exp)-deno_sig+(exp*(tab_dm+tab_neu_0))*np.log(tab_neu_mc*exp) #poisson likelihood of events only beoing caused by background
	p_0_log_corrected = np.nan_to_num(p_0_log)


	p_1_log_tot=np.sum(p_1_log_corrected) #total likelihood after taking product of likelihood in every bin
	p_0_log_tot=np.sum(p_0_log_corrected)

	gauss_1_log = np.log(np.sqrt(2*np.pi*sig_values[3:]**2))-((params_0[3:]-av_values[3:])**2/(2*sig_values[3:]**2)) #gaussian term to penalize neutrino fluxes with values far away from expected value
	gauss_1_log_tot = np.sum(gauss_1_log)

	gauss_0_log = np.log(np.sqrt(2*np.pi*sig_values[3:]**2))-((neu_fluxes_mc-av_values[3:])**2/(2*sig_values[3:]**2))
	gauss_0_log_tot = np.sum(gauss_0_log)

	z_sq = -2*p_0_log_tot-2*gauss_0_log_tot+2*p_1_log_tot+2*gauss_1_log_tot

	if z_sq<0 or z_sq!=z_sq: #error in computation, happens when under neutrino floor
		z_sq = 0

	z= np.sqrt(z_sq)

	

	if details: #show distribution of simulated signal, best fits and likelihood in each bin

		bini = np.arange(len(tab_dm))
		
		#plot number of events in each bin for each category
		plt.hist([tab_dm*exp,tab_neu_0*exp], bini, stacked=True, density=True)
		#plt.plot(bini,tab_dm*exp,'r.',ms=10,label='Simulated DM signal')
		#plt.plot(bini,tab_neu_0*exp,'b.',ms=12,label='Simulated neutrino signal')
		plt.plot(bini,tab_neu_mc*exp,'k.',ms=10,label='Best fit for neutrino only hypothesis')
		plt.yscale('log')
		plt.xlabel('Energy bins')
		plt.ylabel('Number of events')
		plt.text('Expsosure of {:.1f} TY'.format(exp/365/24/3600/1000))
		plt.legend(['Dark Matter','Neutrinos'],loc='best')
		plt.show()
		"""
		#plot likelihoods in every bin
		#plt.plot(bini,p_1_log_corrected,'g.',label='Poisson 1')
		#plt.plot(bini,p_0_log_corrected,'b.',label='Poisson 0')
		plt.plot(bini, 2*(-p_0_log_corrected+p_1_log_corrected),'r.',label='sum')
		plt.ylabel(r'-2log($L_0/L_1$)')
		plt.xlabel('Energy bin #')
		plt.show()
		"""

	return z


def like_0(tab_dm,tab_neu_0,tab_neu_mc,exp,neu_fluxes_mc):
	# background-only likelihood to maximize, return value of the log-likelihood
	# tab_dm simulated DM events in each energy bin
	# tab_neu_0 simulated neutrino events in each bin
	# tab_neu_mc neutrino events in each bin that maximize the background only hypothesis
	# exp exposure in kg-seconds
	# neu_fluxes_mc neutrino fluxes that maximize background-only likelihood
	# details=1 to have detailed figures at every iteration (for debugging only)


	deno_sig=scipy.special.gammaln(tab_dm*exp+tab_neu_0*exp+1)
	
	p_0_log=-(tab_neu_mc*exp)-deno_sig+(exp*(tab_dm+tab_neu_0))*np.log(tab_neu_mc*exp) #poisson likelihood of events only being caused by background

	#total likelihood after taking product of likelihood in every bin
	p_0_log_tot=np.sum(p_0_log)

	gauss_0_log = np.log(np.sqrt(2*np.pi*sig_values[3:]**2))-((neu_fluxes_mc-av_values[3:])**2/(2*sig_values[3:]**2))

	gauss_0_log_tot = np.sum(gauss_0_log)

	#print(p_0_log_tot + gauss_0_log_tot)

	return p_0_log_tot + gauss_0_log_tot


def optimize_this(flux_fit,flux_nofit,tab_dm,tab_neu_0,emin,emax,nbin,neu_arr,er_rej,exp):
	#fonction to MINIMIZE
	# fluxe neutrino fluxes for background only events that maximisee background-only likelihood
	# tab_dm distribution of DM events
	# tab_neu_0 distribution of simulated neutrino events


	fluxes = np.concatenate((flux_fit[:3]*order_fit[:3],flux_nofit,flux_fit[3:]*order_fit[3:])) #multilply value of fluxe (between 0 and 1) by order of magnitude (did this for stability of the minimisation algorithm)

	tab_mc = neutrino_bin(emin,emax,nbin,neu_arr,er_rej,fluxes) #compute neutrino event disttribution 

	like = -like_0(tab_dm,tab_neu_0,tab_mc,exp,fluxes) #compute likelihood

	#print(like)
	return like



def one_mc_run(m,sig,emin,emax,nbin,neu_arr,AA,ZZ,er_rej,acc,exp,speed_dist,dist,dm_events,fac_dm=0.001,fac_nr=1,fac_er=1):
	#run this function for every iteration of the Monte-Carlo

	#randomly pick DM and neu parameters from normal dist.
	#fac factor to include or not DM parameters as nuisance parameters in MC
	dens = np.random.normal(loc=dens_dm,scale=err_dens)

	f_pp = np.random.normal(loc=pp_flux,scale=fac_er*0.01*pp_flux)
	f_hep = np.random.normal(loc=hep_flux,scale=0.16*hep_flux)
	f_8b = np.random.normal(loc=b8_flux,scale=0.16*b8_flux)
	f_15o = np.random.normal(loc= o15_flux,scale=fac_nr*0.05*o15_flux)
	f_13n = np.random.normal(loc= n13_flux,scale=fac_nr*0.05*n13_flux)
	f_17f = np.random.normal(loc= f17_flux,scale=fac_nr*0.05*f17_flux)
	f_pep = np.random.normal(loc=pep_flux,scale=fac_nr*0.02*pep_flux)
	f_be7 = np.random.normal(loc=be7_flux,scale=fac_er*0.105*be7_flux)
	f_atm = np.random.normal(loc=atm_flux,scale=0.2*atm_flux)
	f_dsnb = np.random.normal(loc=dsnb_flux,scale=0.5*dsnb_flux)

	fluxes_neu = np.array([f_pp,f_hep,f_8b,f_15o,f_13n,f_17f,f_pep,f_be7,f_atm,f_dsnb]) #neutrino fluxes
	fluxes_neu_fit = np.array([f_pp,f_hep,f_8b,f_be7,f_atm,f_dsnb]) #fluxes that are fitted (remove CNO and pep because negligeable)
	flux_others = np.array([f_15o,f_13n,f_17f,f_pep]) #other fluxes that are not fitted
	params_0 = np.array([v_esc,v_0,dens,f_pp,f_hep,f_8b,f_15o,f_13n,f_17f,f_pep,f_be7,f_atm,f_dsnb]) #all nuisance parameters

	#compute observed signal distribution from these parameters
	tab_dm = dm_events * dens



	tab_neu_tot = neutrino_bin(emin,emax,nbin,neu_arr,er_rej,fluxes_neu)

	#minimize - likelihood of background hypothesis MAKE THIS FASTER????????
	res = scipy.optimize.minimize(optimize_this,fluxes_neu_fit/order_fit,bounds=bnds,args=(flux_others,tab_dm,tab_neu_tot,emin,emax,nbin,neu_arr,er_rej,exp))
	flux_optim = res.x*order_fit #optimized neutrino fluxes

	flux_good = np.concatenate((flux_optim[:3],flux_others,flux_optim[3:]))


	tab_mc = neutrino_bin(emin,emax,nbin,neu_arr,er_rej,flux_good) #compute event in energy bins for optimized neutrino flux
	
	#compute z value for this MC iteration
	z = likelihood(tab_dm,tab_neu_tot,tab_mc,exp,params_0,flux_good)

	return z


def n_mc_run(n_run,m,sig,emin,emax,nbin,neu_arr,AA,ZZ,er_rej,acc,exp,speed_dist,dist,dm_events_tab):
	#function to perform Monte Carlo over n_run pseudo-experiments
	#m DM mass, sig cross-section
	# emin and emax limits of ROI zone, nbins number of energy bins
	# neu_arr array containing event dist. for all neutrinos NR and ER
	# err_rej ER rejection efficiency, exp exposure in kg s
	# acc function for acceptance as a function of nuclear recoil energy (keV)
	#dm_events_tab 2D arrat containging n_run DM distributions for various DM densities


	z_tab = np.zeros(n_run) #array to contain z value of discovery on each pseudo-experiment
	k=0
	p=0
	l=0
	y=0
	while k<n_run: #loop over n_run pseudo-experiments
		dm_events = dm_events_tab[k] * sig
		z = one_mc_run(m,sig,emin,emax,nbin,neu_arr,AA,ZZ,er_rej,acc,exp,speed_dist,dist,dm_events) #run the MC
		print('MC run #{}, z={}'.format(k,z))
		if not isNaN(z): #if minimisation converged properly
			z_tab[k] = z
			if z_tab[k]<3:
				l+=1
			k+=1
			p=0
		else: #p index to avoid looping indefinetly if minimization not working
			p+=1
		if p>2:
			#print('Minimization algorithm not properly converging for sigma = {}'.format(sig))
			break
		if l>=n_run/10: #already more than 10% of MC pseudo experiment below 3sigma
			result = False
			break
		if z>100: #Z too high
			#print('z too high ({}), go to next iteration'.format(z))
			result = True
			y=1
			break

	#print('z values for this iteration:{}'.format(z_tab))

	if p<=3 and l<n_run/10 and y==0:
		score = np.sum(z_tab>3)/n_run

		if score>0.9: # if more than 90% of pseudo-experiments allow discovery of DM with a significance of at least 3 sigmas
			result = True
		else:
			result = False
	elif p>3:
		result = False #tends no not converge below the floor

	#histogram of z values for this m/cross-section pair
	"""
	
	plt.hist(z_tab, bins='auto')  # arguments are passed to np.histogram
	plt.title("Distribution of significance of detection over simulated experiments")
	plt.show()
	"""


	z_av = np.mean(z_tab)
	#print( 'average z:{}, above floor? {}'.format(z_av,result))
	return result , z_av


def e_theo(m,AA):
	mn = AA*m_p
	red = m*mn/(m+mn)
	e = 2*red**2*(v_0+v_esc+2*err_v0+2*err_vesc)**2/mn
	return e

def roi_def(m,AA,emin,emax,nbin):
	#function to define energy bins according to detector material, min and max energy and maximal number  of  wanted  binsn (nbin)
	e_max_theo = e_theo(m,AA)

	if e_max_theo < emax:
		#cut ROI if maximum energy possible for dm with this mass is below maximal energy
		emax = e_max_theo

	if (emax-emin)/nbin < 0.1:
		#smaller  number of energy bins if size of bins < 0.1 keV
		nbin = int((emax-emin)//0.1)


def dicho_search_sig(search_min,search_max,niter,n_run,m,emin,emax,nbin,neu_arr,AA,ZZ,er_rej,acc,const,speed_dist,dist):
	# dichotomic search to find minimal cross section where at least 90% of pseudo-experiments allow discovery of DM with significance of 3 sigmas
	# for a fixed exposure const

	e_max_theo = e_theo(m,AA)
	#e_max_theo = emax

	if e_max_theo<emin:
		search_mid = 0 #mass too small to be detected given detector threshold

	else:
		"""
		if e_max_theo < emax:
			emax = e_max_theo

		if (emax-emin)/nbin < 0.1:
			nbin = int((emax-emin)//0.1)
			#print(nbin)
		"""
		dm_events_tab = np.zeros((n_run,nbin))

		if 'pp' in speed_dist:
			beta = np.random.randint(6,size=n_run)
			eta = np.random.randint(8,size=n_run)

			for i in range(n_run):

				dm_params = np.array([v_esc,v_0,beta[i],eta[i]])
				dm_events_tab[i,:] = dm_bins(emin,emax,nbin,m,acc,AA,ZZ,dm_params,speed_dist,dist)
		else:
			dm_params = np.array([v_esc,v_0,0,0])
			dm_tab = dm_bins(emin,emax,nbin,m,acc,AA,ZZ,dm_params,speed_dist,dist)

			for i in range(n_run):
				dm_events_tab[i,:] = dm_tab

		sig_min = 1*10**(search_min)
		sig_max = 1*10**(search_max)

		n = 0

		av_tab = np.zeros(niter)
		sig_tab = np.zeros(niter)

		search_mid = (search_max + search_min)/2.

		true_count = 0


		while n<niter:

			print('iteration in bissection search:{}, cross-section{}'.format(n,search_mid+4))

			sc , av_tab[n] = n_mc_run(n_run,m,1*10**(search_mid),emin,emax,nbin,neu_arr,AA,ZZ,er_rej,acc,const,speed_dist,dist,dm_events_tab)
			sig_tab[n] = search_mid
			#print(sc)

			true_count+=int(sc==True)

			if sc:
				search_max = search_mid
			elif not sc:
				search_min = search_mid

			search_mid = (search_max + search_min)/2.

			n+=1
		if true_count==0: #if algorithm never converges to solution, mass probably not detectable within ROI (under threshold)
			search_mid = 0

		#print(av_tab)

	return search_mid


def dicho_search_exp(exp_min,exp_max,niter,n_run,m,emin,emax,nbin,neu_arr,AA,ZZ,er_rej,acc,const,speed_dist,dist):
	# dichotomic search to find minimal exposure where at least 90% of pseudo-experiments allow discovery of DM with significance of 3 sigmas
	# for a fixed cross-section const

	n = 0


	av_tab = np.zeros(niter)
	exp_tab = np.zeros(niter)

	exp_mid = (exp_max + exp_min)/2.

	under = 0


	while n<niter:

		print('iteration in bissection search:{}, exposure:{}'.format(n,exp_mid))

		sc , av_tab[n] = n_mc_run(n_run,m,const,emin,emax,nbin,neu_arr,AA,ZZ,er_rej,acc,365*24*3600*1000*10**exp_mid,speed_dist,dist)
		exp_tab[n] = 365*24*3600*1000*10**exp_mid

		if sc:
			exp_max = exp_mid
		elif not sc:
			exp_min = exp_mid

		exp_mid = (exp_max + exp_min)/2.

		#if exp_mid > 1e9*365*24*3600*100:
		#	under = 1


		#print(n)

		n+=1

	#print(av_tab)
	"""

	plt.plot(exp_tab,av_tab,'k.')
	plt.xlabel(r'$\sigma$')
	plt.ylabel('Average exposure needed of detection')
	plt.yscale('log')
	plt.xscale('log')
	plt.title('{} GeV: floor at {} cm^2'.format(m/1e6,10000*10**search_mid))
	plt.show()
	"""

	return 365*24*3600*1000*10**exp_mid

"""
exp_tab = np.logspace(-2,7,30)

search_min = -54
search_max = -47
niter = 10
n_run = 30
m = 10e6
emin = 1
emax = 90
nbin = 50
neu_arr = np.load('dist_neutrino_xenon_long.npy')
AA = A_xenon_s
ZZ = Z_xenon
er_rej = 10**9
acc = acc_100 
speed_dist = 'pp'
dist = empty

for exp in exp_tab:
	cs = dicho_search_sig(search_min,search_max,niter,n_run,m,emin,emax,nbin,neu_arr,AA,ZZ,er_rej,acc,exp*24*3600*365*1000,speed_dist,dist)
	cross_section = 10000*10**cs

	print('###m={} GeV, cs={} cm2, exposure={} Ty (Xenon, no ER)'.format(m/1e6,cross_section,exp))


er_rej = 1000

for exp in exp_tab:
	cs = dicho_search_sig(search_min,search_max,niter,n_run,m,emin,emax,nbin,neu_arr,AA,ZZ,er_rej,acc,exp,speed_dist,dist)
	cross_section = 10000*10**cs

	print('###m={} GeV, cs={} cm2, exposure={} Ty (Xenon, ER)'.format(m,cross_section,exp))


emin = 15
emax = 110
nbin = 50
neu_arr = np.load('dist_neutrino_argon_long.npy')
AA = A_argon_s
ZZ = Z_argon
er_rej = 10**8

for exp in exp_tab:
	cs = dicho_search_sig(search_min,search_max,niter,n_run,m,emin,emax,nbin,neu_arr,AA,ZZ,er_rej,acc,exp,speed_dist,dist)
	cross_section = 10000*10**cs

	print('###m={} GeV, cs={} cm2, exposure={} Ty (Argon, ER)'.format(m,cross_section,exp))
"""