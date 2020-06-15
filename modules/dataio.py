#!/usr/bin/env python3

"""
With this module one can
  - Make directories if they don't exist with makedir function
  - Create folder for outputs inside the project folder
  - Read in resultsXXX.dat into a DF
  - Read in structures from movieXXX.xyz with ASE
  - Choose a structure as representative
  - Get the number of atoms in the structure
  - Get the list of atom names in the structure
  - Get a unique list of atom names
"""

import os
import pandas as pd
import ase.io

# Get current working directory
wrkdir = os.getcwd()

def makedir(dirname):
    """Function to make directories if they don't exist"""
    try:
        os.mkdir(dirname)
        print('Folder created in {}'.format(dirname))
    except OSError as error:
        print(error)

def init_files(wrkdir):
    path_output = wrkdir + "/plots"
    #print('$PATH_output: {}'.format(path_output))
    makedir(path_output)
    return path_output

def read_data(level,norm):
    """
    Read in XTB of DFT structure info from given filepath.
    Return a dataframe with columns defined in this function.
    Option to normalise and convert to kcal/mol.
    """
    global minEindex

    if level == "XTB":
        # Read in this file
        datafile = wrkdir + "/resultsXTB.dat"
        # Use these columns from the datafiles
        cols = ["path2xyz","GyRadius","XTBEnergy","Dipole"]        # Use these when reading in the XTB data
    elif level == "DFT":
        datafile = wrkdir + "/resultsDFT_freq.dat"
        # Use these columns from the datafiles
        cols = ["path2xyz","GyRadius","DFTEnergy","mrGibbs"]       # Use these when reading in the DFT data

    results_df = pd.read_csv(datafile,sep='\s+',header=None,names=cols)  # Read ResultsXTB into a DF
    results_df["xyzs"] = [x.split('/')[-1] for x in results_df.path2xyz]   # Get only the structure names
    results_df.set_index("xyzs",inplace=True,drop=False)                       # Make the structure names into index-col

    minEindex = results_df.iloc[:,2].idxmin(axis='columns')
    if norm:
        minE = results_df.loc[minEindex][2]
        results_df.iloc[:,2] -= minE
        results_df.iloc[:,2] *= 627.509

    return results_df

def read_xyz(level,indexs):
    """
    Read in structure data from given filepath. Return a dataframe that has all structures as ASE objects
    """
    if level == "XTB":
        # Read in these files
        xyzfile = wrkdir + "/movieXTB.xyz"
    elif level == "DFT":
        xyzfile = wrkdir + "/movieDFT_freq.xyz"
    structures = ase.io.read(xyzfile,index=':')
    structures = pd.DataFrame(structures,index=indexs,columns=['ase'])
    return structures

def get_structure(xyz_df,index_or_name):
    """
    Get the elements of and example structure
    """
    zero_structure = xyz_df.ase[index_or_name]
    if type(index_or_name) == int:
        zero_name = xyz_df.index[index_or_name]

    chemsyms = zero_structure.get_chemical_symbols()
    chemsyms_uniques = list(set(chemsyms))
    n_atoms = len(chemsyms)
    print("🔘 The example structure {} have {} atoms of these elements:".format(zero_name,n_atoms),chemsyms_uniques)
    return zero_structure, zero_name, chemsyms, chemsyms_uniques, n_atoms
