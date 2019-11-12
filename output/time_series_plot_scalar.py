import numpy as np
import pylab as plt
from numpy import*
import matplotlib .pyplot as plt
import time

plt.ion()
plt.rcParams['xtick.major.size'] = 9

plt.rcParams['xtick.major.width'] = 1

plt.rcParams['xtick.minor.size'] = 5

plt.rcParams['xtick.minor.width'] = 1

plt.rcParams['ytick.major.size'] = 9

plt.rcParams['ytick.major.width'] = 1

plt.rcParams['ytick.minor.size'] = 5

plt.rcParams['ytick.minor.width'] = 1

font = {'family' : 'serif', 'weight' : 'normal', 'size' : 15}
plt.rc('font', **font)

#===============FIRST PLOT++++++++++++++++++++++++++++++++++++++++++++++

#Load Data==========================================================
data = np.loadtxt('TimeSeries.dat',comments='#')
t = data[:,0]
Re = data[:,1]
Nu = data[:,2]
Div = data[:,3]
dt = data[:,4]

#Define Figure and subplots=========================================
fig = plt.figure(1)
f_ax1 = fig.add_subplot(2,2,1)
f_ax2 = fig.add_subplot(2,2,2)
f_ax3 = fig.add_subplot(2,2,3)
f_ax4 = fig.add_subplot(2,2,4)
fig.suptitle(r"$Time\:Plots$" ,fontsize=30)


#Sub-Plot 1=============================================================
f_ax1.set_title(r" $timestep\: vs\: time$" ,fontsize=25)
f_ax1.plot(t, dt, "b", lw=2.0, label=r'$dt\: vs\: time$')
	
f_ax1.set_xlabel(r"$t$", fontsize = 20)
f_ax1.set_ylabel(r"$dt$", fontsize = 20)
f_ax1.legend(loc = 0,fontsize=15)

#Sub-Plot 2=============================================================
f_ax2.set_title(r" $Reynolds\:No\:vs\:time$",fontsize=25)
f_ax2.plot(t, Re, "b", lw=2.0, label=r'$Re\: vs\: time$')


f_ax2.set_xlabel("$t$", fontsize = 20)
f_ax2.set_ylabel("$Re$", fontsize = 20)
f_ax2.legend(loc = 0,fontsize=15)

#Sub-Plot 3=============================================================
f_ax3.set_title(r" $Nusselt\:vs\: time$",fontsize=25)
f_ax3.plot(t, Nu, "b", lw=2.0, label=r'$Nu\:vs\:time$')

f_ax3.set_xlabel("$t$", fontsize = 20)
f_ax3.set_ylabel("$Nu$", fontsize = 20)
f_ax3.legend(loc = 0,fontsize=15)

#Sub-Plot 4=============================================================
f_ax4.set_title(r" $Divergence\:vs\:time$",fontsize=25)
f_ax4.plot(t, Div, "b", lw=2.0, label=r'$Div\:vs\:time$')

f_ax4.set_xlabel("$t$", fontsize = 20)
f_ax4.set_ylabel("$Div$", fontsize = 20)
f_ax4.legend(loc = 0,fontsize=15)

#fig.tight_layout()
figManager = plt.get_current_fig_manager()
figManager.resize(*figManager.window.maxsize())
plt.show()

#LOOP===================================================================
while (1) :
	
	#Load Data==========================================================
	data = np.loadtxt('TimeSeries.dat',comments='#')
	t = data[:,0]
	Re = data[:,1]
	Nu = data[:,2]
	Div = data[:,3]
	dt = data[:,4]

	#Define Figure and subplots=========================================
	fig = plt.figure(1)
	f_ax1 = fig.add_subplot(2,2,1)
	f_ax2 = fig.add_subplot(2,2,2)
	f_ax3 = fig.add_subplot(2,2,3)
	f_ax4 = fig.add_subplot(2,2,4)
	fig.suptitle(r"$Time\:Plots$" ,fontsize=30)


	#Sub-Plot 1=============================================================
	f_ax1.plot(t, dt, "b", lw=2.0)
	
	#Sub-Plot 2=============================================================
	f_ax2.plot(t, Re, "b", lw=2.0)

	#Sub-Plot 3=============================================================
	f_ax3.plot(t, Nu, "b", lw=2.0)

	#Sub-Plot 4=============================================================
	f_ax4.plot(t, Div, "b", lw=2.0)
	
	plt.pause(2)
	plt.show()
	

