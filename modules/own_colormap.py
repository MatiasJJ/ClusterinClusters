#!/usr/bin/env python3

"""
Define a colormap for k-means visualisation.
The colors should be distinct.
"""

import matplotlib.pyplot as plt
import numpy as np

def own_cmap(n_clusters):
    colors = [
        "dimgrey",
        "red",
        "indianred",
        "darkorange",
        "gold",
        "greenyellow",
        "darkgreen",
        "limegreen",
        "cyan",
        "dodgerblue",
        "royalblue",
        "blue",
        "blueviolet",
        "plum",
        "deeppink",
        "magenta",
        "purple",
        "saddlebrown",
        "lightsalmon",
        "k"]
        # "mediumblue"
        # "silver"


    if n_clusters > len(colors):
        n_last = n_clusters - len(colors)
        print("Last {} clusters don't have colors!".format(n_last))

    return colors[:n_clusters]

def visualise_colors(n_clusters):
    """Visually inspect colors"""
    arrrr = np.arange(5,5+n_clusters)
    for i in range(arrrr.shape[0]):
        plt.bar(arrrr[i],arrrr[i],color=own_cmap(n_clusters)[i])
    plt.xticks(ticks=arrrr,labels=list(arrrr-5))
    plt.show()
