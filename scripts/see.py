#!/usr/bin/env python3
import glob

import h5py
import matplotlib.pyplot as plt

files = glob.glob('hdf5_out_2d_*.h5')


def printname(name):
    print(name)


fig = plt.figure()
for file in files:
    # file = 'hdf5_out_2d_00000.h5'

    f = h5py.File(file, 'r')

    print(f['Bx'])

    f.visit(printname)
    bx = f['Bx']
    ax = plt.subplot(2, 3, 1)
    im = ax.pcolor(f['Velx'])
    ax = plt.subplot(2, 3, 2)
    im = ax.pcolor(f['Vely'])
    ax = plt.subplot(2, 3, 3)
    im = ax.pcolor(f['Bx'])
    ax = plt.subplot(2, 3, 4)
    im = ax.pcolor(f['By'])
    ax = plt.subplot(2, 3, 5)
    im = ax.pcolor(f['Density'])
    ax = plt.subplot(2, 3, 6)
    im = ax.pcolor(f['Energy'])
    fig.canvas.draw()
    plt.savefig(file + '.png', dpi=100)
