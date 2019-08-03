#This code interpolates all the fields to staggered gridpoints.


import h5py
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt

def hdf5_reader(filename,dataset):
    file_V1_read = h5py.File(filename)
    dataset_V1_read = file_V1_read["/"+dataset]
    V1=dataset_V1_read[:,:,:]
    return V1



Ra = 1e7
Pr=1.0

kappa = np.sqrt(1.0/(Ra*Pr))
nu = np.sqrt(Pr/Ra)

lx, ly, lz = 1.0, 1.0, 1.0

U = hdf5_reader("Soln_0027.0000.h5", "Vx")
V = hdf5_reader("Soln_0027.0000.h5", "Vy")
W = hdf5_reader("Soln_0027.0000.h5", "Vz")

T = hdf5_reader("Soln_0027.0000.h5", "T")
P = hdf5_reader("Soln_0027.0000.h5", "P")

U_st = np.zeros_like(P)
V_st = np.zeros_like(P)
W_st = np.zeros_like(P)

Nx_s = len(P[:,0,0])
Nx_c = len(U[:,0,0])
Ny_s = len(P[0,:,0])
Ny_c = len(V[0,:,0])
Nz_s = len(P[0,0,:])
Nz_c = len(W[0,0,:])

#print Nx_s, Nx_c

x = np.linspace(0,1,Nx_s)
y = np.linspace(0,1,Ny_s)
z = np.linspace(0,1,Nz_s)

U_st[1:Nx_c,:,:] = 0.5*(U[0:Nx_c-1,:,:]+U[1:Nx_c,:,:])
U_st[0,:,:] = 0.0
U_st[-1,:,:] = 0.0

V_st[:,1:Ny_c,:] = 0.5*(V[:,0:Ny_c-1,:]+V[:,1:Ny_c,:])
V_st[:,0,:] = 0.0
V_st[:,-1,:] = 0.0

W_st[:,:,1:Nz_c] = 0.5*(W[:,:,0:Nz_c-1]+W[:,:,1:Nz_c])
W_st[:,:,0] = 0.0
W_st[:,:,-1] = 0.0


f1 = h5py.File("SolnI_0027.0000.h5", "w")
dset1 = f1.create_dataset("Vx", data = U_st)
dset2 = f1.create_dataset("Vy", data = V_st)
dset3 = f1.create_dataset("Vz", data = W_st)
dset4 = f1.create_dataset("T", data = T)
dset5 = f1.create_dataset("P", data = P)
f1.close()

UzT = W_st*T_st


av_UzT = np.mean(UzT)
Nu = 1.0 + av_UzT/(kappa/lz)
print "Nu", Nu

Uh_rms = U_st**2 + V_st**2
Uh_z_avg = np.zeros([Nz_s])
for i in range(Nz_s):
    Uh_z_avg[i] = np.sqrt(np.mean(Uh_rms[:,:,i]))

plt.plot(Uh_z_avg, z)
plt.show()



