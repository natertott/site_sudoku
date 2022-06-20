from scipy import optimize
#from pathLib import Path
#from typing import List
import numpy as np
import os

# imports
import sys
import glob

from objective_functions import total_agreement
from utilities import electron_number

# append the path to your GSASII installation directory here
sys.path.append('/Users/natej/g2full/GSASII/')
import GSASIIscriptable as G2sc

def main():

    # define some directories for keeping things clean
    workPath = os.getcwd()
    gpxPath = workPath + '/workDir/'
    dataPath = workPath + '/dataDir/'
    exportPath = workPath + '/exportDir/'

    # define a path to the cif that you want to use as the original
    cifPath = dataPath + 'I-4m2.cif'

    # get the histogram csv files
    fpaths = glob.glob(dataPath + '*.csv')
    # only want the first one for testing
    file = fpaths[0]

    # define the species in the material and the number of electrons they each have
    species = ['Nb', 'Co', 'Sn']
    z = np.array([41, 27, 50])  # number of electrons

    # define the composition of the material (theoretical or measured by WDS)
    expected_compo = [0.33, 0.33, 0.33]  # in atomic percent

    # calculate the expected electron count based on the species and composition
    expected_e = electron_number(expected_compo, z)

    # define all the site you want to place atoms in the unit cell
    sites = [[0.00000, 0.00000, 0.00000],
             [0.50000, 0.50000, 0.50000],
             [0.00000, 0.00000, 0.50000],
             [.50000, 0.50000, 1.00000],
             [0.00000, 0.50000, 0.25000],
             [0.50000, 1.00000, 0.75000],
             [0.00000, 0.50000, 0.75000],
             [0.50000, 1.00000, 1.25000]]

    n_sites = len(sites)

    # create an initial set of occupancies for each atom at each site
    # for starters place one atom of each type at the first n many sites
    n = len(sites)
    m = len(species)

    init = np.zeros([n, m])
    # initial site doping
    # put 1/n percent of each species on each site
    for ii in range(n):
        for jj in range(m):
            init[ii, jj] = np.round(1 / n, 4)

    init = np.reshape(init, m * n)

    # define some bounds for the solution
    # constrain the occupancies to be between 0 and 2.5
    bnds = m * n * [[0, 2.5]]

    result = optimize.minimize(
        lambda x: total_agreement(x,expected_e,expected_compo,z,species,sites),
        x0 = init,  # Initial guess
        method='Powell',
        bounds = bnds,
        options = {'maxiter':1}
    )

    print(result)

if __name__ == '__main__':
    main()