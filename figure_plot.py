# coding=utf-8
# Andrea Gaspert 05/2019
# code that takes in information about neutrino floor vs DM mass in the form of a npy file, plots it and saves the plot.

import numpy as np
import matplotlib.pyplot  as plt
import sys
from scipy import interpolate


#sys.path.insert(0, './data')


#MODES NOT UPDATED YET< TO DO!!!

####### MODES ##################################

#Mode 'reg' is to plot a regular neutrino floor, for a given fixed exposure, with the option of showing experimental data or not
#Mode 'exp' is to plot the neutrino floor vs exposure for a fixed DM mass

#############################################


#neu_floor = np.load('floor_xenon_ER.npy') #xenon floor with ER
def plot_floor(neu_floor, save_name, exp_show=True):
	#neu_floor numpy array containing floor with data-type [('mass', float), ('cs', float), ('exp', float)]
	m = neu_floor['mass']
	cs = neu_floor['cs']
	exp = neu_floor['exp']

	if exp_show == True:
		#EXCLUSION CURVES

		#read files
		read_1=np.genfromtxt('data/xenon1t_2018.txt')
		read_2=np.genfromtxt('data/pandax_2017.txt')
		read_3=np.genfromtxt('data/new_g_2017.txt')
		read_4=np.genfromtxt('data/lux_2017.txt')
		read_5=np.genfromtxt('data/deap_2018.txt')
		read_6=np.genfromtxt('data/dark_side_2018.txt')
		read_7=np.genfromtxt('data/cresst_2015.txt')
		read_8=np.genfromtxt('data/cdmslite2_2015.txt')
		read_9=np.genfromtxt('data/darkside_lowmass2018.txt')

		read_neu=np.genfromtxt('data/neu_floor.txt')


		dama1_data = np.genfromtxt('data/dama1.txt')
		dama2_data = np.genfromtxt('data/dama2.txt')
		dama3_data = np.genfromtxt('data/dama3.txt')
		dama4_data = np.genfromtxt('data/dama4.txt')

		m_1, ex_1 = read_1[:,0],read_1[:,1]*1e-47
		m_2, ex_2 = read_2[:,0],read_2[:,1]*1e-47
		m_3, ex_3 = read_4[:,0],read_4[:,1]*1e-47
		m_4, ex_4 = read_5[:,0],read_5[:,1]*1e-47
		m_5, ex_5 = read_6[:,0],read_6[:,1]*1e-47
		m_6, ex_6 = read_7[:,0],read_7[:,1]*1e-47
		m_7, ex_7 = read_8[:,0],read_8[:,1]*1e-47
		m_8, ex_8 = read_3[:,0],read_3[:,1]*1e-47
		m_9, ex_9 = read_9[:,0],read_9[:,1]*1e-47

		m_neu , ex_neu = read_neu[:,0],read_neu[:,1]*1e-47

		m_dama1 , ex_dama1 = dama1_data[:,0],dama1_data[:,1]*1e-47
		m_dama2 , ex_dama2 = dama2_data[:,0],dama2_data[:,1]*1e-47
		m_dama3 , ex_dama3 = dama3_data[:,0],dama3_data[:,1]*1e-47
		m_dama4 , ex_dama4 = dama4_data[:,0],dama4_data[:,1]*1e-47

		x2_dama = np.linspace(7,16.5,100)
		x1_dama = np.linspace(29,115,100)

		fit1 = interpolate.interp1d(m_dama1,ex_dama1)	
		y1_dama = fit1(x1_dama)
		fit2 = interpolate.interp1d(m_dama2,ex_dama2)
		y2_dama = fit2(x1_dama)
		fit3 = interpolate.interp1d(m_dama3,ex_dama3)
		y3_dama = fit3(x2_dama)
		fit4 = interpolate.interp1d(m_dama4,ex_dama4)
		y4_dama = fit4(x2_dama)

		mass_tab = [m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9]
		ex_tab = [ex_1,ex_2,ex_3,ex_4,ex_5,ex_6,ex_7,ex_8,ex_9]
		color_tab = ['darkorange','brown','darkorange','steelblue','steelblue','olivedrab','olivedrab','olivedrab','peru']
		text_tab = ['XENON1T 2018','PandaX-II 2017','LUX 2017', 'DEAP-3600 2018','Dark-Side 2018','CRESST 2015','CDMSLite2 2015','NEWS-G 2017','Dark-Side Low Mass 2018']
		style_tab = ['-',':','--','-',':','--',':','-','-']
		pos_x = 1000
		ini = 1.5e-37
		pos_y = [ini,ini/4,ini/16,ini/4**3,ini/4**4,ini/4**5,ini/4**6,ini/4**7]

		for i in range(8):
			plt.plot(mass_tab[i],ex_tab[i],linewidth=2,color=color_tab[i],linestyle=style_tab[i],label=text_tab[i])


		plt.fill_between(x1_dama,y1_dama,y2_dama,color='salmon',label='DAMA/LIBRA')
		plt.fill_between(x2_dama,y3_dama,y4_dama,color='salmon')

	plt.plot(m,cs,'r.',ms=12, label = 'Neutrino floor Xenon with ER')
	plt.xlabel(r'WIMP mass $[GeV/c^2]$',fontsize=14)
	plt.ylabel(r'Nucleon-wimp cross-section $[cm^{2}]$', fontsize=14)
	plt.xlim([0.001,1e5])
	plt.ylim([1e-50,1e-40])
	plt.xscale('log')
	plt.yscale('log')
	plt.legend(loc='best')

	plt.savefig(save_name+'.png')