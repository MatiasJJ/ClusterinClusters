#!/usr/bin/env python3

"""
With this module one can
    - Calculate k-means labels to find the similar structures
    - Calculate average energies of structures for all clusters
    - Filter out higher energy structures by choosing only the lowest average energy clusters
"""

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

def calcKmeans(desdf,n_clusters,steve,v):
    """
    Runs k-means for specified descriptor set
    returns cluster assignments for all structures in a same order that structures are in DataFrames
    """
    kmeans = KMeans(n_clusters=n_clusters,n_init=12,n_jobs=steve,algorithm='full')#,random_state=42) # n_init means amount of different random initializations
    # Compute k-means for descs
    kmeans.fit(desdf)
    # Predict the closest cluster each sample in desdf belongs to
    kmeans_results = kmeans.predict(desdf)

    if v: print("""🔘 Clustering done""")

    return kmeans_results

def calcEAvg(datadf,level,descname,n_clusters):
    """
    Calculate average energies for all clusters.
    Return a DataFrame of cluster labels and average energies.
    """
    # Group clusters into group-object
    cluster_groups = datadf.groupby(by="{}_{}_klabel".format(level,descname))
    cluster_labels = list(range(n_clusters))
    EAvg_df = pd.DataFrame(index=cluster_labels,columns=["EAvg"])

    for grp in range(n_clusters):
        cluster = cluster_groups.get_group(grp)
        EAvg_df.loc[grp,:] = cluster["{}Energy".format(level)].mean()
    return EAvg_df.astype(np.float64)

def getBestClusters(datadf,level,descname,n_clusters,nk_out,ns_out,sampl):
    """
    Take nk_out best cluster-labels and return the list of structures in those clusters.
    If sampl=True then return ns_out randomly sampled structures from the nk_out best clusters.
    """

    EAvg_df = calcEAvg(datadf, level, descname, n_clusters)
    best_k = EAvg_df.nsmallest(nk_out,["EAvg"]).index.tolist()
    # Group clusters into group-object
    cluster_groups = datadf.groupby(by="{}_{}_klabel".format(level,descname))

    best_list = []
    for grp in best_k:
        cluster = cluster_groups.get_group(grp)
        indices = cluster.index.tolist()
        best_list.extend(indices)

    best_df = datadf[datadf.index.isin(best_list)]

    if sampl:
        best_df = best_df.sample(ns_out)
        print("🔘 Structure selection done with random sampling.")
    else:
        print("🔘 Structure selection done.")

    return best_df
