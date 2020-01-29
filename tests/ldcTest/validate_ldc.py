#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import h5py as hp
import yaml as yl

# Pyplot-specific directives
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = 'cm'
plt.rcParams["font.weight"] = "medium"

ptFile = False

def init():
    global U, W
    global Nx, Nz
    global figSize
    global xLen, zLen

    Nx = 0
    Nz = 0
    xLen = 0.0
    zLen = 0.0

    U = np.zeros([1, 1, 1])
    W = np.zeros([1, 1, 1])

    figSize = (12, 6)


def parseYAML(paraFile):
    global Nx, Nz
    global xLen, yLen, zLen

    yamlFile = open(paraFile, 'r')
    yamlData = yl.load(yamlFile)

    Nx = 2**yamlData["Mesh"]["X Index"] + 1
    Nz = 2**yamlData["Mesh"]["Z Index"] + 1

    xLen = yamlData["Program"]["X Length"]
    zLen = yamlData["Program"]["Z Length"]


def loadGhia():
    global u_ghia, v_ghia

    u_ghia = np.loadtxt("u_profile_ghia.dat", comments='#')
    v_ghia = np.loadtxt("v_profile_ghia.dat", comments='#')


def loadData(timeVal):
    global Nx, Nz
    global U, W

    fileName = "output/Soln_{0:09.4f}.h5".format(float(timeVal))
    print("Processing file " + fileName + "\n")

    try:
        f = hp.File(fileName, 'r')
    except:
        print("Could not open file " + fileName + "\n")
        exit()

    # Initialize and read staggered grid data
    U = np.zeros([Nx, Nz])
    U = np.array(f['Vx'])

    W = np.zeros([Nx, Nz])
    W = np.array(f['Vz'])


def getVorticity():
    global Nx, Nz

    wy = np.zeros((Nx-2, Nz-2))

    hx = 1.0/(Nx - 1)
    hz = 1.0/(Nz - 1)

    for i in range(1, Nx-1):
        for k in range(1, Nz-1):
            dux_dz = (U[i, k+1] - U[i, k-1])*0.5/hz
            duz_dx = (W[i+1, k] - W[i-1, k])*0.5/hx

            wy[i-1, k-1] = dux_dz - duz_dx

    return wy


def plotProfile():
    global U, W
    global Nx, Nz
    global u_ghia, v_ghia

    wy = getVorticity()

    if ptFile:
        plt.switch_backend('agg')

    # Plot a data frame
    fig, axes = plt.subplots(1, 2, figsize=figSize)

    uProfile = U[int(Nx/2), :]
    profAxis = np.linspace(0.0, zLen, Nz)
    axes[0].plot(u_ghia[:,2], u_ghia[:,1], marker='*', markersize=10, linestyle=' ', label='Ghia et al')
    axes[0].plot(uProfile, profAxis, linewidth=2, label='SARAS')
    axes[0].set_xlim([-0.6, 1.1])
    axes[0].tick_params(labelsize=20)
    axes[0].legend(fontsize=20)
    axes[0].set_title(r"$v_x$ at $z=0.5$", fontsize=22)

    vProfile = W[:, int(Nz/2)]
    profAxis = np.linspace(0.0, xLen, Nx)
    axes[1].plot(v_ghia[:,1], v_ghia[:,2], marker='*', markersize=10, linestyle=' ', label='Ghia et al')
    axes[1].plot(profAxis, vProfile, linewidth=2, label='SARAS')
    axes[1].set_ylim([-0.62, 0.42])
    axes[1].tick_params(labelsize=20)
    axes[1].legend(fontsize=20)
    axes[1].set_title(r"$v_z$ at $x=0.5$", fontsize=22)

    plt.gca().set_aspect('auto')
    plt.tight_layout()

    if ptFile:
        plt.savefig("test.png")
    else:
        plt.show()


if __name__ == "__main__":
    init()

    parseYAML("input/parameters.yaml")

    loadData(30.0)

    loadGhia()

    plotProfile()

