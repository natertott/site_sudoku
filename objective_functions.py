import numpy as np
import os, sys
from utilities import write_cif

# append the path to your GSASII installation directory here
sys.path.append('/Users/natej/g2full/GSASII/')
import GSASIIscriptable as G2sc

# define a path to the cif that you want to use as the original
workPath = os.getcwd()
dataPath = workPath + '/dataDir/'
cifPath = dataPath + 'I-4m2.cif'

def total_agreement(occupancy_guess,
                    expected_el_count,
                    expected_composition,
                    z,
                    species,
                    sites,
                    electron_weight = 0.5,
                    composition_weight = 1,
                    n = 8,
                    m = 3,
                    reitveld_weight = 5):
    """Compute the overall agreement across many factors

    Smaller weight values mean to consider disagreement for that factor more heavily.

    Args:
        occupancies: Flattened array of the site occupanices of site x elements
        expected_el_count: Expected number of electrons on each site
        expected_composition: Expected overall composition for each element
        actual_pattern: Path to the actual diffraction pattern
        z: Atomic number of each element
        electron_weight: Expected level of disagreement in number of atoms
        composition_weight: Expected level of disagreement in composition fraction
        reitveld_weight:
    Returns:
        Overall agreement across all of the factors, weighted
    """

    # n and m are defined in a separate cell for now
    # for working out the initial solution
    # ideally n and m would be passed along into the function call
    occupancy_reshape = np.reshape(occupancy_guess, [n, m])

    return composition_agreement(occupancy_reshape, expected_composition) / composition_weight + \
           electron_agreement(occupancy_reshape, expected_el_count, z) / electron_weight + \
           rietveld_agreement(occupancy_guess, species, sites)['wR'] / reitveld_weight


# now we try a rietveld refinement

def rietveld_agreement(occupancies, species, sites,
                       blank_cif_path = cifPath,):
    """
    Measures the agreement between an experimental diffraction histogram and a simulated histogram. The simulated histogram
    is simulated based on a provided blank cif file, a list of elements in the unit cell, and the occupancies of each
    element on each site in the unit cell.
    :param occupancies: N by M array of atomic occupancies at each site. Entry [0,0] is the first element on the first
    site, [0,1] is the first element on the second site, [1,0] is the second element on the first site, etc.
    :param actual_pattern_path: Path to a histogram csv file for the experimental histogram
    :param z: List of electron numbers for each element.
    :param blank_cif_path: path to the blank cif file used for simulation
    :return wR: the weighted residual between the experimental and simulated histograms
    :return data: the experimental and simulated histogram arrays
    :return occ: the occupancies of the unit cell
    """

    # write the cif file with the occupancies
    ncif = write_cif(dataPath, blank_cif_path, species, occupancies, sites, title='mix')

    # hard coding a few things for the sake of testing
    projname = 'test_gpx'

    # create a project
    gpx = G2sc.G2Project(newgpx=os.path.join(workPath, projname + '.gpx'))

    # add the histograms
    inst = '/Users/natej/Documents/TECCA/sudoku/basin_hopping/data/BL1-5_Feb2022.instprm'  # path to instrument parameter file
    histPath = '/Users/natej/Documents/TECCA/sudoku/basin_hopping/data/NP018_pad2.csv'
    thist = gpx.add_powder_histogram(histPath, inst)  # hist is the path to the histogram file

    # set the number of refinement cycles
    gpx.data['Controls']['data']['max cyc'] = 1  # not in API

    # now we want to add a phase
    # for now define a single phase path
    p = os.path.join(dataPath, ncif)  # path to the .cif file with the phase
    p0 = gpx.add_phase(p, phasename='I-4m2', histograms=[0])

    # set the microstrain for the phase
    p0.setSampleProfile(0, 'micro', 'isotropic', 1000)

    # the instrument parameters are read in from the instrument parameters file
    # we need to refine the unit cell, background, and scale factor
    refdict = [{'set': {
        'Background': {'no.coeffs': 9, 'refine': True},
        'Sample Parameters': {'Scale': 1.3805996}}
    }]

    gpx.do_refinements(refdict)

    # get the weighted residual
    wR = thist.get_wR()

    return {'wR': wR}


def composition_agreement(occupancies, expected_composition,n=8,m=3):
    """Compute the agreement between predicted and actual composition

    Args:
        occupancies: Flattened array of the site occupanices of site x elements
        expected_composition: Expected overall composition for each element
        z: Atomic number of each element
    Returns:
        Degree of mismatch
    """

    occupancies = np.reshape(occupancies, [n, m])

    # Measure the overall composition
    pred_comp = occupancies.sum(axis=0)

    diff = np.subtract(pred_comp, expected_composition)
    return np.power(diff, 2).mean()


def electron_agreement(occupancies, expected_el_count, z, n=8,m=3):
    """Compute the electron count agreement

    Args:
        occupancies: Flattened array of the site occupanices of site x elements
        expected_el_count: Expected number of electrons on each site
        z: Atomic number of each element
    Returns:
        Degree of mismatch
    """

    occupancies = np.reshape(occupancies, [n, m])

    print(z)
    pred_count = 0.01*np.sum(occupancies, axis=0)*z  # calculate the number of electrons

    # Compute the difference
    diff = np.subtract(expected_el_count, pred_count)
    return np.power(diff, 2).mean()