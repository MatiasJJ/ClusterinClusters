#!/usr/bin/env python3

"""
With this module one can
    - Calculate dendrograms for the structures
    - Visualize clustering results with t-SNE in 2D
"""

import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import scipy.cluster.hierarchy as sch
from sklearn.manifold import TSNE
from mpl_toolkits.mplot3d import Axes3D
from ase.io import write

def makeDend(desdf,indexs,level,descname,steve,v,path_output,onlyShow=False):
    #plt.ioff() # to avoid showing plots in Jupyter
    plt.figure(figsize=(16,6))
    plt.title("{} {}'s Dendrogram".format(level,descname))
    plt.xlabel('Structures as {}s'.format(descname))
    plt.ylabel('Difference')
    #plt.grid(True)
    dendrogram = sch.dendrogram(sch.linkage(desdf, method = 'ward'))

    if onlyShow:
        plt.show()
    else:
        plt.savefig("{}/{}_{}_dend.png".format(path_output,level,descname),dpi=400)
        plt.close()
        if v: print("🔘 A dendrogram plotted in {}/{}_{}_dend.png".format(path_output,level,descname))

def makeTsne_2D(desdf,indexs,level,descname,n_clusters,steve,v,path_output,clusterlabels,own_cmap,clustermethod="kmeans",pp=30.0,onlyShow=False):
    tsne = TSNE(n_components=2,perplexity=pp).fit_transform(desdf)

    tsne_df = pd.DataFrame(tsne,index=indexs,columns=['tsne1','tsne2'])
    tsne_df["klabel"] = clusterlabels

    # Group t_SNE results into group-object
    tsne_ID_groups = tsne_df.groupby(by='klabel')
    # Plot one group at a time
    for grp in range(n_clusters):
        plot_group = tsne_ID_groups.get_group(grp)
        plot_center = plot_group.mean()
        # Plot t-SNE results for one group
        plt.scatter(plot_group.tsne1,plot_group.tsne2, c=own_cmap[grp], alpha=0.3, s=5)
        # and the centers
        plt.scatter(plot_center.tsne1, plot_center.tsne2, c=own_cmap[grp], alpha=0.9, s=200, marker='*', label=grp)
    plt.legend(loc='upper left')
    plt.title("{} {}'s {} t-SNE".format(level,descname,clustermethod))
    plt.ylabel('t-SNE component 2')
    plt.xlabel('t-SNE component 1')

    if onlyShow:
        plt.show()
    else:
        plt.savefig("{}/{}_{}_p{}_tsne.png".format(path_output,level,descname,int(pp)),dpi=400)
        plt.close()

    if v: print("🔘 Finished t-SNE 2D for visualization for all {} {}'s""".format(tsne.shape[0],descname))

    return tsne_df

def plotTsneE_3D(datadf,level,descname,v,path_output,Efilter,own_cmap,clustermethod="kmeans",onlyShow=False,selected=False):
    """
    Visualise clustering results with t-SNE and Energy as the third axis.
    Saves an interactive visualisation as .html (open in browser)
    """

    if not selected:
        datadf = datadf[datadf["{}Energy".format(level)]<=Efilter].sort_values(by="{}_{}_klabel".format(level,descname))

    klabels = datadf["{}_{}_klabel".format(level,descname)].astype(str)
    fig = px.scatter_3d(data_frame=datadf,
              z="{}Energy".format(level),
              x="{}_{}_tsne1".format(level,descname),
              y="{}_{}_tsne2".format(level,descname),
              color=klabels,
              color_discrete_sequence=own_cmap,
              size="GyRadius",
              opacity=0.9,
              #symbol="symbol",   # Use if needed in Jupyter
              hover_name=datadf.index,
              title="{}'s' t-SNE + {}Energy".format(descname,level),
              #range_z=[-36,-20],
              width= 1200,
              height= 900)


    if onlyShow:
        fig.show()
    elif selected:
        fig.write_html("{}/{}_{}_EtSNE_selected.html".format(path_output,level,descname))
    else:
        fig.write_html("{}/{}_{}_EtSNE.html".format(path_output,level,descname))

def struct2img(atoms,indexs,path_output):
    """
    input: ASE structure
    output: Structure visualization
    """
    write('{}/{}.png'.format(path_output,indexs[:-4]),atoms)
