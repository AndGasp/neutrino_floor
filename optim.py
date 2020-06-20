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

order_fit = np.array([1e7,1e10,1e4,1e5,1e14,1e13])
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

	e_tab = np.logspace(np.log10(emin),np.log10(emax),nbin+1)


	dm_events = np.zeros(nbin)
	dm_events_new = np.zeros(nbin)

	if 'reg' in speed_dist: #if want to use SHM

		ee = np.linspace(emin,emax,50000)

		RR = int_new(ee,m,AA,ZZ,acc,iso_1,vesc,v0,1,1,'none',0,0)


		for i in range(nbin-1):
			"""

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
		"""
		tabb = np.linspace(0,100,nbin)
		plt.plot(tabb,dm_events)
		plt.yscale('log')
		plt.show()
		"""

	elif 'pp' in speed_dist: #if want to use SHM++

		ee = np.linspace(emin,emax,20000)
		dR = np.zeros(20000)

		#start1 = time.time()
		for k in range(20000):
			dR[k] = int_irreg(ee[k],m,AA,acc,beta_ind,eta_ind)

		for i in range(nbin):
			ind_min = np.argmin(np.abs(ee-e_tab[i]))
			ind_max = np.argmin(np.abs(ee-e_tab[i+1])) 	
			dm_events_new[i] = int_simpson(ee[ind_min:ind_max],dR[ind_min:ind_max])
		#time1 = time.time() - start1
		"""
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


def neutrino_bin(emin,emax,nbin,neu_array):
	#Compute number of neutrino events in each bin per kg per second for given detector properties (NR and ER)
	#emin and emax ROI limits, nbin number of energy bins
	#neu_array 6xn array containing normalized distribution of neutrino events as a function of energy for all solar, atm and DSNB neutrinos, NR and ER

	e_neu = neu_array[0,:]

	new_array = np.zeros((6,nbin))

	#fit distributions for different neutrino fluxes
	fit_hep=scipy.interpolate.interp1d(e_neu,neu_array[1,:],kind='linear',bounds_error = False, fill_value = 0)
	fit_b8=scipy.interpolate.interp1d(e_neu,neu_array[2,:],kind='linear',bounds_error = False, fill_value = 0)
	fit_atm=scipy.interpolate.interp1d(e_neu,neu_array[3,:],kind='linear',bounds_error = False, fill_value = 0)
	fit_dsnb=scipy.interpolate.interp1d(e_neu,neu_array[4,:],kind='linear',bounds_error = False, fill_value = 0)
	fit_pp=scipy.interpolate.interp1d(e_neu,neu_array[5,:],kind='linear',bounds_error = False, fill_value = 0)
	fit_be7=scipy.interpolate.interp1d(e_neu,neu_array[6,:],kind='linear',bounds_error = False, fill_value = 0)

	ee = np.linspace(emin,emax,50000)
	
	events_hep = fit_hep(ee)
	events_b8 = fit_b8(ee)
	events_atm = fit_atm(ee)
	events_dsnb = fit_dsnb(ee)
	events_pp = fit_pp(ee)
	events_be7 = fit_be7(ee)


	e_tab = np.logspace(np.log10(emin),np.log10(emax),nbin+1) 

	for i in range(nbin):
		"""

		start1 = time.time()

		r_neu[i] = scint.quad(fit_neu_nr,e_tab[i],e_tab[i+1])[0] + scint.quad(fit_neu_er,e_tab[i],e_tab[i+1])[0]

		time1 = time.time()- start1
		"""
		#start2 = time.time()

		ind_min = np.argmin(np.abs(ee-e_tab[i]))
		ind_max = np.argmin(np.abs(ee-e_tab[i+1]))
			
		new_array[0,i] = int_simpson(ee[ind_min:ind_max],events_hep[ind_min:ind_max])
		new_array[1,i] = int_simpson(ee[ind_min:ind_max],events_b8[ind_min:ind_max])
		new_array[2,i] = int_simpson(ee[ind_min:ind_max],events_atm[ind_min:ind_max])
		new_array[3,i] = int_simpson(ee[ind_min:ind_max],events_dsnb[ind_min:ind_max])
		new_array[4,i] = int_simpson(ee[ind_min:ind_max],events_pp[ind_min:ind_max])
		new_array[5,i] = int_simpson(ee[ind_min:ind_max],events_be7[ind_min:ind_max])

		#time2 = time.time()-start2
		#print('comparison between integral method neutrino:{} / {} vs {} / {}'.format(r_neu[i],time1,r_neu_new[i],time2))

		#print(i)


	return new_array


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

	gauss_1_log = np.log(np.sqrt(2*np.pi*sig_values_new**2))-((params_0-av_values_new)**2/(2*sig_values_new**2)) #gaussian term to penalize neutrino fluxes with values far away from expected value
	gauss_1_log_tot = np.sum(gauss_1_log)

	gauss_0_log = np.log(np.sqrt(2*np.pi*sig_values_new**2))-((neu_fluxes_mc-av_values_new)**2/(2*sig_values_new**2))
	gauss_0_log_tot = np.sum(gauss_0_log)

	z_sq = -2*p_0_log_tot-2*gauss_0_log_tot+2*p_1_log_tot+2*gauss_1_log_tot

	if z_sq<0 or z_sq!=z_sq: #error in computation, happens when under neutrino floor
		z_sq = 0

	z= np.sqrt(z_sq)

	

	if details:

		bini = np.arange(len(tab_dm))
		
		#plot number of events in each bin for each category
		plt.plot(bini,tab_dm*exp,'r.',ms=10,label='Simulated DM signal')
		plt.plot(bini,tab_neu_0*exp,'b.',ms=12,label='Simulated neutrino signal')
		plt.plot(bini,(tab_neu_0+tab_dm)*exp,'k.',ms=10,label='Simulated total signal')
		#plt.plot(bini,tab_neu_mc*exp,'r.',ms=6,label='Best fit for neutrino only events')
		plt.yscale('log')
		plt.legend(loc='best')
		plt.show()
		

		#plot number of events in each bin for each category
		#plt.plot(bini,tab_dm*exp,'r.',ms=6,label='Simulated DM signal')
		#plt.plot(bini,tab_neu_0*exp,'b.',ms=12,label='Simulated neutrino signal')
		plt.plot(bini,(tab_neu_0+tab_dm)*exp,'k.',ms=12,label='Simulated total signal')
		plt.plot(bini,tab_neu_mc*exp,'r.',ms=10,label='Best fit for neutrino only hypothesis')
		plt.yscale('log')
		plt.legend(loc='best')
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

	gauss_0_log = np.log(np.sqrt(2*np.pi*sig_values_new**2))-((neu_fluxes_mc-av_values_new)**2/(2*sig_values_new**2))

	gauss_0_log_tot = np.sum(gauss_0_log)


	return p_0_log_tot + gauss_0_log_tot


def optimize_this(flux_fit,tab_dm,tab_neu_0,neu_arr,err_rej,exp):
	#fonction to MINIMIZE
	# fluxe neutrino fluxes for background only events that maximisee background-only likelihood
	# tab_dm distribution of DM events
	# tab_neu_0 distribution of simulated neutrino events


	fluxes = flux_fit*order_fit #multilply value of fluxe (between 0 and 1) by order of magnitude (did this for stability of the minimisation algorithm)
	tab_mc = neu_arr[0,:]*fluxes[0] + neu_arr[1,:]*fluxes[1] + neu_arr[2,:]*fluxes[2] + neu_arr[3,:]*fluxes[3] + neu_arr[4,:]*fluxes[4]/er_rej +neu_arr[5,:]*fluxes[5]/er_rej

	like = -like_0(tab_dm,tab_neu_0,tab_mc,exp,fluxes) #compute likelihood

	return like



def one_mc_run(m,sig,emin,emax,nbin,neu_arr,AA,ZZ,er_rej,acc,exp,speed_dist,dist,dm_events,fac_dm=0.001,fac_nr=1,fac_er=1):
	#run this function for every iteration of the Monte-Carlo

	#randomly pick DM and neu parameters from normal dist.
	#fac factor to include or not DM parameters as nuisance parameters in MC
	dens = np.random.normal(loc=dens_dm,scale=err_dens)

	#6 neutrino sources that will be considered (pp and be7 for neutrino-electron recoils only, other for NR only)
	f_pp = np.random.normal(loc=pp_flux,scale=fac_er*0.01*pp_flux)
	f_hep = np.random.normal(loc=hep_flux,scale=0.16*hep_flux)
	f_8b = np.random.normal(loc=b8_flux,scale=0.16*b8_flux)
	f_be7 = np.random.normal(loc=be7_flux,scale=fac_er*0.105*be7_flux)
	f_atm = np.random.normal(loc=atm_flux,scale=0.2*atm_flux)
	f_dsnb = np.random.normal(loc=dsnb_flux,scale=0.5*dsnb_flux)

	fluxes_neu = np.array([f_hep,f_8b,f_atm,f_dsnb,f_pp,f_be7]) #neutrino fluxes
	fluxes_neu_fit = np.array([f_hep,f_8b,f_atm,f_dsnb,f_pp,f_be7]) #fluxes that are fitted (remove CNO and pep because negligeable)
	params_0 = np.array([f_hep,f_8b,f_atm,f_dsnb,f_pp,f_be7]) #initial guess for minimisation algo.

	#compute observed signal distribution from these parameters
	tab_dm = dm_events * dens

	#add all sources of neutrinos to get total events in every bin for actual neutrino fluxes
	tab_neu_tot = neu_arr[0,:]*f_hep + neu_arr[1,:]*f_8b + neu_arr[2,:]*f_atm + neu_arr[3,:]*f_dsnb + neu_arr[4,:]*f_pp/er_rej +neu_arr[5,:]*f_be7/er_rej


	#minimize - likelihood of background hypothesis 
	res = scipy.optimize.minimize(optimize_this,fluxes_neu_fit/order_fit,bounds=bnds,args=(tab_dm,tab_neu_tot,neu_arr,er_rej,exp))
	flux_optim = res.x*order_fit #optimized neutrino fluxes

	#add all sources of neutrinos to get total events in every bin for neutrino fluxes that fit DM signal best
	tab_mc = neu_arr[0,:]*flux_optim[0] + neu_arr[1,:]*flux_optim[1] + neu_arr[2,:]*flux_optim[2] + neu_arr[3,:]*flux_optim[3] + neu_arr[4,:]*flux_optim[4]/er_rej + neu_arr[5,:]*flux_optim[5]/er_rej #compute event in energy bins for optimized neutrino flux
	
	#compute z value for this MC iteration
	z = likelihood(tab_dm,tab_neu_tot,tab_mc,exp,params_0,flux_optim)

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
		dm_events = dm_events_tab * sig
		z = one_mc_run(m,sig,emin,emax,nbin,neu_arr,AA,ZZ,er_rej,acc,exp,speed_dist,dist,dm_events) #run the MC
		#print('MC run #{}, z={}'.format(k,z))
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

		if 'pp' in speed_dist:
			beta = np.random.randint(6,size=n_run)
			eta = np.random.randint(8,size=n_run)

			for i in range(n_run):

				dm_params = np.array([v_esc,v_0,beta[i],eta[i]])
				dm_events_tab[i,:] = dm_bins(emin,emax,nbin,m,acc,AA,ZZ,dm_params,speed_dist,dist)
		else:
			dm_params = np.array([v_esc,v_0,0,0])
			dm_events_tab = dm_bins(emin,emax,nbin,m,acc,AA,ZZ,dm_params,speed_dist,dist) #DM events in every bin / by dm flux

		neu_arr_new = neutrino_bin(emin,emax,nbin,neu_arr) #array containing neutrino events per bin for each neutrino source, 


		sig_min = 1*10**(search_min)
		sig_max = 1*10**(search_max)

		n = 0

		av_tab = np.zeros(niter)
		sig_tab = np.zeros(niter)

		search_mid = (search_max + search_min)/2.

		true_count = 0


		while n<niter:

			#print('iteration in bissection search:{}, cross-section{}'.format(n,search_mid+4))

			sc , av_tab[n] = n_mc_run(n_run,m,1*10**(search_mid),emin,emax,nbin,neu_arr_new,AA,ZZ,er_rej,acc,const,speed_dist,dist,dm_events_tab)
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

		#print('iteration in bissection search:{}, exposure:{}'.format(n,exp_mid))

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


	return 365*24*3600*1000*10**exp_mid

def empty(e):
	return 1


#==========================TEST=================================================
"""
m_tab = np.logspace(0.5,4,40)

search_min = -55
search_max = -47
niter = 12
n_run = 1000
emin = 4.9
emax = 90
nbin = 50
neu_arr = np.load('dist_neutrino_xenon_new.npy')
AA = A_xenon_s
ZZ = Z_xenon
er_rej = 10**2
acc = acc_100 
speed_dist = 'reg'
dist = empty
exp = 1e4

n_try = 1
print('%12 iteration in bi, 1000 mc runs, Xe with 100 rej and 1000TY')

for m in m_tab:
	for n in range(n_try):
		cs = dicho_search_sig(search_min,search_max,niter,n_run,m*1e6,emin,emax,nbin,neu_arr,AA,ZZ,er_rej,acc,exp*24*3600*365*1000,speed_dist,dist)
		cross_section = 10**(cs+4)
		print('###m={} GeV, cs={} cm2, exposure={} Ty'.format(m,cross_section,exp))

"""
