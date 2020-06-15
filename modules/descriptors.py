#!/usr/bin/env python3

"""
With this module one can
    - Define descriptor Hyperparameters
    - Create MBTR for all structures
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import ase.data
from dscribe.descriptors import CoulombMatrix,MBTR,SOAP
from itertools import combinations_with_replacement # for SOAP plotting NOTE: I am not sure if this is the correct way of doing it

# MBTR parameters
mk1={
        "geometry": {"function": "atomic_number"},
        "grid": {"min": 0, "max": 18, "n": 100, "sigma": 0.6}, # org "sigma" : 0.1
    }
mk2={
        "geometry": {"function": "inverse_distance"},
        "grid": {"min": 0, "max": 2.5, "n": 100, "sigma": 0.08}, # org "sigma#: 0.1, "scale": 0.5, "cutoff": 1e-3
        "weighting": {"function": "exponential", "scale": 0.6, "cutoff": 1e-3},
    }
mk3={
        "geometry": {"function": "cosine"},
        "grid": {"min": -1, "max": 1, "n": 100, "sigma": 0.03}, # org "sigma#: 0.1, "scale": 0.5, "cutoff": 1e-3
        "weighting": {"function": "exponential", "scale": 0.2, "cutoff": 1e-3},
    }

# k1 not relevant for clustering clusters
# sigma: gaussian broadening
# cutoff: to speed up calcs
# scale: weighs the atoms that are further away. 1 is the equal weight
# n: keep at 100 unless want to reduce n_features

# SOAP parameters
srcut = 6.0  # A cutoff for local region in angstroms. Should be bigger than 1 angstrom.
snmax = 8    # The number of radial basis functions.
slmax = 8    # The maximum degree of spherical harmonics.
ssigma = 0.4 # org 1.0, The standard deviation of the gaussians used to expand the atomic density.
srbf = 'gto'   # The radial basis functions to use.
            # The available options are:
            # ”gto”: Spherical gaussian type orbitals defined as 𝑔𝑛𝑙(𝑟)=∑𝑛max𝑛′=1𝛽𝑛𝑛′𝑙𝑟𝑙𝑒−𝛼𝑛′𝑙𝑟2
            # ”polynomial”: Polynomial basis defined as 𝑔𝑛(𝑟)=∑𝑛max𝑛′=1𝛽𝑛𝑛′(𝑟−𝑟cut)𝑛′+2

# 1. rcut, sets the cutoff for atoms whose gaussian densities will be included in the integral.
# 2. nmax, sets the number of orthogonal radial basis functions to use.
# 3. lmax, sets the number of angular momentum terms, so l = 0, 1, ...., lmax

# Create MBTR or CM
def setupDescs(structs,indexs,level,descname,chemsyms_uniques,n_atoms,steve,v):
    """
    Setup descriptor and run it for ASE structures.
    Return DataFrame with given strictures as descriptors
    """
    # choose the descriptor
    if descname == "CM":
        desc = CoulombMatrix(n_atoms_max=n_atoms, flatten=True)
        # permutation = 'sorted_l2' is default
        n_feat = desc.get_number_of_features()

    if descname == "MBTR":
        desc = MBTR(species=chemsyms_uniques,k1=mk1,k2=mk2,k3=mk3,periodic=False,normalization="l2_each",flatten=True)
        n_feat = desc.get_number_of_features()

    if descname == "SOAP":
        desc = SOAP(species=chemsyms_uniques,periodic=False,rcut=srcut,nmax=snmax,lmax=slmax,average=True) # Averaging for global
        n_feat = desc.get_number_of_features()

    # Create descriptors
    descs = desc.create(structs,n_jobs=steve) # Parallel

    # Create a DF of returned `list` of `arrays` of descs
    descs_df = pd.DataFrame(descs, index=indexs)

    if v: print("""🔘 Created {}-descriptors for all {} {}-structures.
    Number of features in {}: {}""".format(descname,structs.shape[0],level,descname,n_feat))

    return descs_df,n_feat

def plotDescs(structs,indexs,level,descname,chemsyms,n_atoms,steve,v,path_output,save=True):
    """
    Plot descriptors
    """
    # choose the descriptor
    if descname == "CM":
        desc = CoulombMatrix(n_atoms_max=n_atoms, flatten=False, permutation='none') # permutation = 'sorted_l2' is default
        n_feat = desc.get_number_of_features()
        # Create descriptors
        descs = desc.create(structs,n_jobs=steve) # Parallel
        # Plot CM of zero_cluster and save it to outputs-folder
        sns.heatmap(descs,cmap='Spectral',robust=True,xticklabels=chemsyms,yticklabels=chemsyms)
        plt.title("CM of {}".format(indexs))
        if save: plt.savefig("{}/{}_CM.png".format(path_output,indexs[:-4]))

    if descname == "MBTR":
        desc = MBTR(species=list(set(chemsyms)),k1=mk1,k2=mk2,k3=mk3,periodic=False,normalization="l2_each",flatten=False)
        n_feat = desc.get_number_of_features()
        descs = desc.create(structs, n_jobs=steve) # Parallel
        # Create the mapping between an index in the output and the corresponding chemical symbol
        n_elements = len(desc.species)
        # dict({index_of_atom_type:Z_of_atom_type})
        imap = desc.index_to_atomic_number
        # dict({index_of_atom_type:atom_type_symbol})
        smap = {index: ase.data.chemical_symbols[number] for index, number in imap.items()}

        # Plot k=1
        x = np.linspace(0, 1, 100)          # las number defines the resolution of x-axis
        x1 = desc.get_k1_axis()             # from fullmetalfelix/ML-CSC-tutorial
        fig, ax = plt.subplots()
        for i in range(n_elements):
            plt.plot(x1, descs["k1"][i, :], label="{}".format(smap[i]))
        ax.set_xlabel("Charge")
        ax.set_xlabel("Atomic number")#, size=20) # from fullmetalfelix/ML-CSC-tutorial
        ax.set_ylabel("k1 values (arbitrary units)")#, size=20)
        plt.legend()
        plt.title("MBTR k1 of {}".format(indexs))
        if save: plt.savefig("{}/{}_MBTR_k1.png".format(path_output,indexs[:-4]))

        # Plot k=2
        x = np.linspace(0, 0.5, 100)  # Kato mitä tää on docsista
        x2 = desc.get_k2_axis()            # from fullmetalfelix/ML-CSC-tutorial
        fig, ax = plt.subplots()
        for i in range(n_elements):
            for j in range(n_elements):
                if j >= i:
                    plt.plot(x2, descs["k2"][i, j, :], label="{}-{}".format(smap[i], smap[j]))
        ax.set_xlabel("Inverse distance (1/angstrom)")#, size=20) # How to make not inverse?
        ax.set_ylabel("k2 values (arbitrary units)")#, size=20)
        plt.legend()
        plt.title("MBTR k2 of {}".format(indexs))
        if save: plt.savefig("{}/{}_MBTR_k2.png".format(path_output,indexs[:-4]))

        # Plot k=3
        x = np.linspace(0, 0.5, 100)  # Kato mitä tää on docsista
        x3 = desc.get_k3_axis()            # from fullmetalfelix/ML-CSC-tutorial
        fig, ax = plt.subplots()
        for i in range(n_elements):
            for j in range(n_elements):
                if j >= i:
                    for k in range(n_elements):
                        if k >= j and smap[k]=="S":
                            plt.plot(x3, descs["k3"][i, j, k, :], label="{}-{}-{}".format(smap[i], smap[j], smap[k]))
        ax.set_xlabel("cos(angle)")#, size=20)
        ax.set_ylabel("k3 values (arbitrary units)")#, size=20)
        plt.legend()
        plt.title("MBTR k3 of {}".format(indexs))
        if save: plt.savefig("{}/{}_MBTR_k3.png".format(path_output,indexs[:-4]))

    if descname == "SOAP":
        desc = SOAP(species=list(set(chemsyms)),periodic=False,rcut=srcut,nmax=snmax,lmax=slmax,average=False) # Averaging for global
        n_feat = desc.get_number_of_features()
        descs = desc.create(structs, n_jobs=steve)

        # Plot SOAPs for all atom pairs
        chemsyms_combos = list(combinations_with_replacement(desc.species,2))
        for combo in chemsyms_combos:
            # The locations of specific element combinations can be retrieved like this.
            pairloc = desc.get_location(combo)
            # These locations can be directly used to slice the corresponding part from an
            # SOAP output for e.g. plotting.
            plt.plot(descs[0, pairloc], label="{}-{}".format(combo[0],combo[1]))
        plt.legend()
        #plt.xlim(20,40)
        plt.xlabel("N of features for an atom pair")
        plt.ylabel("Output value of SOAPs")
        plt.title("SOAPs of {}".format(indexs))
        if save: plt.savefig("{}/{}_SOAP.png".format(path_output,indexs[:-4]))

    if v: print("🔘 Plotting {} done.".format(descname))
