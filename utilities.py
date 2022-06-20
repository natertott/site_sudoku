import os
import numpy as np

def write_cif(dirt, ciffile, species, occ, sites, title='mix'):
    r = open(os.path.join(dirt, ciffile), 'r')

    # now create a new cif to write to
    newcif = title + '.cif'
    w = open(os.path.join(dirt, newcif), 'w')

    # doing this kind of stupidly, probably
    for ind, row in enumerate(r):
        if '# ATOMIC COORDINATES AND DISPLACEMENT PARAMETERS' in row:
            break
        else:
            w.write(row)

    for ind, row in enumerate(r):
        if 'loop_' in row:
            break
        else:
            continue

    # okay so we wrote all the existing lines
    # now we want to write the unit cell information
    w.write('# ATOMIC COORDINATES AND DISPLACEMENT PARAMETERS\n')
    w.write('loop_\n')
    w.write('   _atom_site_label\n')
    w.write('   _atom_site_type_symbol\n')
    w.write('   _atom_site_fract_x\n')
    w.write('   _atom_site_fract_y\n')
    w.write('   _atom_site_fract_z\n')
    w.write('   _atom_site_occupancy\n')
    w.write('   _atom_site_adp_type\n')
    w.write('   _atom_site_U_iso_or_equiv\n')
    w.write('   _atom_site_site_symmetry_multiplicity\n')

    # number of atoms you have
    n = len(species)
    # number of sites you have
    m = len(sites)

    ocounter = 0
    # iterate by the atom type
    for ii in range(n):
        ncounter = 0

        # iterate by the site
        for jj in range(m):
            s = species[ii]  # the species of the atom
            x = sites[jj]  # the site coordinates
            name = str(s) + str(ncounter)  # the name of the atom
            o = occ[ocounter]  # the occupancy of the atom on the site

            # write the line
            newline = f'{name}\t{s}\t{x[0]}\t{x[1]}\t{x[2]}\t{o}\t' + 'Uiso\t0.0100\t2\n'
            w.write(newline)

            ocounter += 1
            ncounter += 1

    for ind, row in enumerate(r):
        if 'loop_' in row:
            w.write('\n')
            w.write(row)
            break
        else:
            continue

    for s, o in zip(species, occ):
        total = np.sum(o)
        newline = f'    {s}\t{str(total)}\n'
        w.write(newline)

    for ind, row in enumerate(r):
        if '_' in row or '#' in row:
            w.write('\n')
            w.write(row)
            break
        else:
            continue

    for ind, row in enumerate(r):
        w.write(row)

    # finally close the files
    r.close()
    w.close()

    # return the location of the new cif file
    return newcif

# calculate the scattering from each site based on the refinement site speices
# that is, multiply the scattering factor by the atomic number
def site_scattering(occ,Z):
    F = []
    for O,charge in zip(occ,Z):
        F.append(np.sum(O)*charge)
    return F

# calculate the total number of electrons in the material
# based on how many of each species are in the material (composition) and the number of electrons per species
def electron_number(compo,Z):
    e = []
    for C, charge in zip(compo,Z):
        e.append(0.01*C*charge)
    return(e)