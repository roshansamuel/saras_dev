import numpy as np
import h5py
#from pylab import*
#import matplotlib as m
import matplotlib.pyplot as plt
from matplotlib import ticker, colors
from mpl_toolkits.mplot3d import axes3d
from matplotlib import gridspec
from scipy.integrate import simps, trapz
#from mayavi import mlab

def hdf5_reader(filename,dataset):
    file_V1_read = h5py.File(filename)
    dataset_V1_read = file_V1_read["/"+dataset]
    V1=dataset_V1_read[:,:]
    return V1

kappa = 0.0001
d = 0.1

A=20
font = {'family' : 'serif', 'weight' : 'normal', 'size' : A}
plt.rc('font', **font)
B = 19
print("CHECK")

t = 43
ypos = 17
def plot_density():
#  U = hdf5_reader("Soln_{0:009.4f}.h5".format(t), "Vx")
#  V = hdf5_reader("Soln_{0:009.4f}.h5".format(t), "Vy")
#  W = hdf5_reader("Soln_{0:009.4f}.h5".format(t), "Vz")

  
  
  
  
  fName = "SolnIP_0005.0000.h5"
  slc = 64
  g = 16

  U = hdf5_reader("U.V1r.h5", "U.V1r")#[:,slc,:]
  V = hdf5_reader("U.V2r.h5", "U.V2r")#[:,slc,:] 
  W = hdf5_reader("U.V3r.h5", "U.V3r")#[:,slc,:]
  
  T = hdf5_reader("T.Fr.h5", "T.Fr")#[:,slc,:]
  
  [Nx, Nz] = T.shape

  print(T.shape, U.shape, V.shape, W.shape) 
  
  L = 1.0
  x = np.linspace(0,1,Nx)# loadtxt("mesh.d")[:,0]
  z = np.linspace(0,1,Nz) #loadtxt("mesh.d")[:,0]  
  
  
  
  fig, axes = plt.subplots(figsize=(10,10))
  levels = []

  
  density = axes.pcolor(x, z, np.transpose(T), cmap='jet')#, norm = colors.LogNorm())
  Q = axes.quiver(x[::g], z[::g], np.transpose(U[::g,::g]), np.transpose(W[::g,::g]), units='width')
  
  
  axes.set_aspect(1)
  axes.set_xticks([0, 0.5, 1.0])
  axes.set_yticks([0, 0.5, 1.0])
  axes.set_xlabel('x/d')
  axes.set_ylabel('z/d')
  axes.tick_params(axis='x', which='major', pad=10)
  cb1 = fig.colorbar(density, fraction=0.05, ax=axes)#, ticks=[1e-2, 1e0, 1e2]) ###### TICKS FOR THE COLORBARS ARE DEFINED HERE
  cb1.ax.tick_params(labelsize=A)
  fig.tight_layout()
  plt.show()
  
plot_density()
