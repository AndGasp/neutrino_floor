# coding=utf-8
# Andrea Gaspert 05/2019
# code that takes in information about neutrino floor vs DM mass in the form of a npy file, plots it and saves the plot.

import numpy as np
import matplotlib.pyplot  as plt
import sys
import scipy.interpolate

sys.path.insert(0, './data')

####### MODES ##################################

#Mode 'reg' is to plot a regular neutrino floor, for a given fixed exposure, with the option of showing experimental data or not
#Mode 'exp' is to plot the neutrino floor vs exposure for a fixed DM mass
#

#############################################
exp_show = True

neu_ar_floor = np.load('floor_final_Ar_ER_1e6_reg.npy') #Argon new floor with ER
neu_xe_floor = np.load('floor_final_Xe_ER_1e4_reg.npy') #Xenon new floor with ER

def ready_to_plot(file_name,m_tab):
	#computes fit of neutrino floor from computation results in file_name for masses in m_tab 
	
	neu_floor = np.load(file_name)
	m = neu_floor['mass']
	cs = neu_floor['cs']
	ind = np.where(cs==10000)
	m = np.delete(m,ind)
	cs = np.delete(cs,ind)
	exp = neu_ar_floor['exp']

	print(m,cs)


	fit_floor = scipy.interpolate.interp1d(m,cs,kind='cubic',bounds_error = False, fill_value = 0)

	cs_new = fit_floor(m_tab)
	ind_bad = np.where(cs_new==0)
	m_new = np.delete(m_tab,ind_bad)
	cs_new = np.delete(cs_new,ind_bad)

	return m_new,cs_new


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

	fit1 = scipy.interpolate.interp1d(m_dama1,ex_dama1)	
	y1_dama = fit1(x1_dama)
	fit2 = scipy.interpolate.interp1d(m_dama2,ex_dama2)
	y2_dama = fit2(x1_dama)
	fit3 = scipy.interpolate.interp1d(m_dama3,ex_dama3)
	y3_dama = fit3(x2_dama)
	fit4 = scipy.interpolate.interp1d(m_dama4,ex_dama4)
	y4_dama = fit4(x2_dama)

	mass_tab = [m_1,m_2,m_3,m_4,m_5,m_6,m_9]
	ex_tab = [ex_1,ex_2,ex_3,ex_4,ex_5,ex_6,ex_9]
	color_tab = ['darkorange','brown','darkorange','steelblue','steelblue','olivedrab','olivedrab','olivedrab','peru']
	text_tab = ['XENON1T 2018','PandaX-II 2017','LUX 2017', 'DEAP-3600 2018','Dark-Side 2018','NEWS-G 2017','Dark-Side Low Mass 2018']
	style_tab = ['-',':','--','-',':','--',':','-','-']
	pos_x = 1000
	ini = 1.5e-37
	pos_y = [ini,ini/4,ini/16,ini/4**3,ini/4**4,ini/4**5,ini/4**6,ini/4**7]
	"""
	plt.figure(figsize=(15,15))
	plt.fill_between(x1_dama,y1_dama,y2_dama,color='peachpuff',label='DAMA/LIBRA')
	plt.fill_between(x2_dama,y3_dama,y4_dama,color='peachpuff')
	plt.plot(x1_dama,y1_dama,linewidth=18,color='peachpuff')
	plt.plot(x1_dama,y2_dama,linewidth=18,color='peachpuff')
	plt.plot(x2_dama,y3_dama,linewidth=18,color='peachpuff')
	plt.plot(x2_dama,y4_dama,linewidth=18,color='peachpuff')
	"""

	for i in range(7):
		plt.plot(mass_tab[i],ex_tab[i],linewidth=4,color=color_tab[i],linestyle=style_tab[i],label=text_tab[i])

m_tab = np.logspace(0,5,100)

m_ar,cs_ar = ready_to_plot('floor_final_Ar_ER_1e6_reg.npy',m_tab)
m_xe,cs_xe = ready_to_plot('floor_final_Xe_ER_1e4_reg.npy',m_tab)
m_xe_noer,cs_xe_noer = ready_to_plot('floor_final_Xe_noER_1e4_reg.npy',m_tab)
m_gaia,cs_gaia = ready_to_plot('floor_final_Xe_ER_1e4_pp_lowmass.npy',m_tab)

plt.plot(m_ar,cs_ar,'r--',linewidth= 6, label = 'Neutrino floor Ar ({:.1E} TY)'.format(1000000))
plt.plot(m_xe,cs_xe,'k--', linewidth=6, label = 'Neutrino floor Xe ({:.1E} TY)'.format(10000))

plt.xlabel(r'WIMP mass $[GeV/c^2]$',fontsize=30)
plt.ylabel(r'Nucleon-wimp cross-section $[cm^{2}]$', fontsize=30)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.xlim([0.5,3e4])
plt.ylim([1e-51,1e-38])
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best',fontsize=18)
plt.legend(bbox_to_anchor=(0.55, 0.6), prop={'size': 24})

fig_name = 'floor_final_money.png'
plt.savefig(fig_name)
plt.show()


"""
mm =np.linspace(1,200,250)
xx = np.ones(250)*1e-45
plt.figure(figsize=(15,9))
plt.plot(m_xe,cs_xe,'k--',linewidth=8,label = 'SHM')
plt.plot(m_tab,cs_fit_g,color='deepskyblue',linestyle='-',linewidth=8,label=r'Profiling over anisotropy')
plt.fill_between(mm,xx,fit_1(mm),color='grey',alpha=0.2)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlim([3,180])
plt.ylim([5e-50,4e-46])
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'WIMP mass $[GeV/c^2]$',fontsize=24)
plt.ylabel(r'Nucleon-wimp cross-section $[cm^{2}]$', fontsize=24)
plt.legend(loc='best',fontsize=18)
plt.text(40, 1e-46, 'Excluded Parameter Space', fontsize=18, rotation=12,color='grey')
plt.savefig('gaia_effect.png')
plt.show()
"""