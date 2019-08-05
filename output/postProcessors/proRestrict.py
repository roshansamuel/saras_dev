import numpy as np
import h5py




def hdf5_reader(filename,dataset):
    file_V1_read = h5py.File(filename)
    dataset_V1_read = file_V1_read["/"+dataset]
    V1=dataset_V1_read[:,:,:]
    return V1



#Reduces the size of the array to a lower level, 2^(n-1)+1.
def restrict(function):
    global sInd, sLst

    
    restricted = np.zeros([rx + 1, ry + 1, rz + 1])

    for i in range(2, rx-1):
        for j in range(2, ry-1):
            for k in range(2, rz-1):
                restricted[i, j, k] = function[2*i - 1, 2*j - 1, 2*k - 1]

    return restricted


#Increases the size of the array to a higher level, 2^(n+1)+1.
def prolong(function):
    [lx, ly, lz] = np.shape(function)
    rx, ry, rz = 2*(lx-1), 2*(ly-1), 2*(lz-1)
    prolonged = np.zeros([rx + 1, ry + 1, rz + 1])
    print "Initializing the prolonging operation"
    for i in range(0, rx+1, 2):
        print i, "out of", rx
        for j in range(0, ry+1, 2):
            for k in range(0, rz+1, 2):
                prolonged[i, j, k] = function[i/2, j/2, k/2]
    
    print "Prolonging along x direction"
    for i in range(1, rx, 2):
        print i, "out of", rx-1
        for j in range(0, ry+1, 2):
            for k in range(0, rz+1, 2):
                prolonged[i, j, k] = (prolonged[i-1, j, k] + prolonged[i+1, j, k])/2.0

    print "Prolonging along y direction"
    for i in range(0, rx+1):
        print i, "out of", rx
        for j in range(1, ry, 2):
            for k in range(0, rz+1):
                prolonged[i, j, k] = (prolonged[i, j-1, k] + prolonged[i, j+1, k])/2.0

    print "Prolonging along z direction"
    for i in range(0, rx+1):
        print i, "out of", rx
        for j in range(0, ry+1):
            for k in range(1, rz, 2):
                prolonged[i, j, k] = (prolonged[i, j, k-1] + prolonged[i, j, k+1])/2.0
                
    return prolonged



T = hdf5_reader("Soln_0005.0000.h5", "T")
T_p = prolong(T)
#print T_p.shape
print np.mean(T), np.mean(T_p)
#print T[25, 38, 50]
print T_p[121, 79, 121]
print (T[60, 39, 60] + T[61, 39, 60] + T[60, 39, 61] + T[61, 39, 61] + T[60, 40, 60] + T[61, 40, 60] + T[60, 40, 61] + T[61, 40, 61])/8.0





