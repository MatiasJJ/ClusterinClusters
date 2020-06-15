#!/usr/bin/env python3
"""
Main program for clustering-based structure selection for JKCS

### This script should be run inside the same folder where JKCS is run ###

Before running make sure that input4Clustering.csv -file is filled as you prefer.
"""

# Imports
import os                   # for finding files
import pandas as pd
import distutils.util       # for parsing string to boolean
from timeit import default_timer as timer # for timing

# Modules
import modules.own_colormap as cmap
import modules.dataio as dataio
import modules.descriptors as dd
import modules.selection as sel
import modules.visualize as viz

# Init timer
start = timer()
print("🥳 Clusterin'Clusters start! 🥳")

# Parse arguments
args = pd.read_csv("input4Clustering.csv",index_col='parameter',sep='|')

# save args into nice variables
n_jobs = int(args.loc['n'][0])
n_clusters_init = int(args.loc['k'][0])
n_clusters_out = int(args.loc['c'][0])
M = int(args.loc['m'][0])
sampl = distutils.util.strtobool(args.loc['r'][0])
n_structures_out = int(args.loc['s'][0])
level = args.loc['l'][0]
normEd = distutils.util.strtobool(args.loc['e'][0])
verbose = distutils.util.strtobool(args.loc['v'][0])
descname = args.loc['d'][0]
plotD = distutils.util.strtobool(args.loc['pd'][0])
plotC = distutils.util.strtobool(args.loc['pc'][0])

# use custom colormap
#print(cmap.own_cmap(n_clusters_init))   # Print the names of the colors
#cmap.visualise_colors(n_clusters_init)  # Plot the colors as barplot

# Get current working directory
wrkdir = os.getcwd()
#print('$WRKDIR: {}'.format(wrkdir))

if verbose: print("🥳 Make folder for plots.")
# Make output-folder for plots
path_output = dataio.init_files(wrkdir)
if verbose: print("⏱  Time elapsed: {} sec".format(timer()-start))

if verbose: print("🥳 Read in {}-data and structures.".format(level))
data_df = dataio.read_data(level,normEd)
structure_index = data_df.index
xyz_df = dataio.read_xyz(level,structure_index)

if verbose: print("🥳 Here's some info for you:")
zero_structure, zero_name, chemsyms, chemsyms_uniques, n_atoms = dataio.get_structure(xyz_df,0)
if verbose: print("⏱  Time elapsed: {} sec".format(timer()-start))

if verbose: print("🥳 Set up descriptors:")
descs_df, n_feat = dd.setupDescs(xyz_df.ase, structure_index, level, descname, chemsyms_uniques, n_atoms, n_jobs, verbose)
if verbose: print("⏱  Time elapsed: {} sec".format(timer()-start))

if plotC:
    if verbose: print("🥳 Make dendrogram")
    viz.makeDend(descs_df, structure_index, level, descname, n_jobs, verbose, path_output)
    if verbose: print("⏱  Time elapsed: {} sec".format(timer()-start))

if verbose: print("🥳 Start clustering:")
kmeans_results = sel.calcKmeans(descs_df, n_clusters_init, n_jobs, verbose)
data_df["{}_{}_klabel".format(level, descname)] = kmeans_results
if verbose: print("⏱  Time elapsed: {} sec".format(timer()-start))

if plotC:
    if verbose: print("🥳 Calculate t-SNE:")
    tsne_df = viz.makeTsne_2D(descs_df, structure_index, level, descname, n_clusters_init, n_jobs, verbose, path_output, kmeans_results, cmap.own_cmap(n_clusters_init))
    data_df["{}_{}_tsne1".format(level, descname)] = tsne_df.iloc[:,0]
    data_df["{}_{}_tsne2".format(level, descname)] = tsne_df.iloc[:,1]
    viz.plotTsneE_3D(data_df, level, descname, verbose, path_output, 5*M, cmap.own_cmap(n_clusters_init))
    if verbose: print("⏱  Time elapsed: {} sec".format(timer()-start))

if verbose: print("🥳 Calc cluster average energies and perform structure selection:")
EAvg_df = sel.calcEAvg(data_df, level, descname, n_clusters_init)
best_df = sel.getBestClusters(data_df, level, descname, n_clusters_init, n_clusters_out, n_structures_out, sampl)
if verbose: print("⏱  Time elapsed: {} sec".format(timer()-start))

if plotD:
    if verbose: print("🥳 Plotting example descriptor into plots-folder")
    best_index = best_df.iloc[:,2].idxmin(axis='columns')
    dd.plotDescs(xyz_df.ase.loc[best_index], best_index, level, descname, chemsyms, n_atoms, n_jobs, verbose, path_output)
    viz.struct2img(xyz_df.ase.loc[best_index], best_index, path_output)
    if verbose: print("⏱  Time elapsed: {} sec".format(timer()-start))

best_df.to_csv("{}/selected{}.csv".format(wrkdir,level),sep='\t',index=False)
if verbose:
    print("🔘 The list of selected structures are found in {}/selected{}_{}.csv".format(wrkdir,level,descname))

print("🥳 Clusterin'Clusters done! 🥳")
if verbose: print("⏱  Total time elapsed: {} sec".format(timer()-start))
