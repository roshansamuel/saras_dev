#This code interpolates all the fields back from staggered gridpoints
import h5py
import numpy as np


def hdf5_reader(filename,dataset):
    file_V1_read = h5py.File(filename)
    dataset_V1_read = file_V1_read["/"+dataset]
    V1=dataset_V1_read[:,:,:]
    return V1

fileName = "SolnIP_0300.0000.h5"

U = hdf5_reader(fileName, "Vx")
V = hdf5_reader(fileName, "Vy")
W = hdf5_reader(fileName, "Vz")

T = hdf5_reader(fileName, "T")
P = hdf5_reader(fileName, "P")


[Nx, Ny, Nz] = T.shape

U_p = (U[0:Nx-1,:,:] + U[1:Nx,:,:])/2.0
V_p = (V[:,0:Ny-1,:] + V[:,1:Ny,:])/2.0
W_p = (W[:,:,0:Nz-1] + W[:,:,1:Nz])/2.0

f1 = h5py.File("SolnP_0300.0000.h5", "w")
dset1 = f1.create_dataset("Vx", data = U_p)
dset2 = f1.create_dataset("Vy", data = V_p)
dset3 = f1.create_dataset("Vz", data = W_p)
dset4 = f1.create_dataset("T", data = T)
dset5 = f1.create_dataset("P", data = P)
f1.close()


