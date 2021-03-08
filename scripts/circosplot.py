'''
slight edit of original code from ponnhide
https://github.com/ponnhide/pyCircos
'''

# to allow importing to work correctly (in a dirty way)
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import constants
from scripts import calculate_beta_correlation
from scripts.pyCircos import Gcircle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
import re
import pandas as pd


def get_links_and_nodes(corr_df, pivot_corr_df, pattern, sign_threshold):
    # calculate mean corr across methods
    links = pivot_corr_df['corr'].mean(axis=1).reset_index().rename(columns={0:'corr'}).copy() 
    # only get corr if in gwas_group_dict
    links = links[links[['gwasx','gwasy']].apply(lambda x: x.str.contains(pattern)).all(axis=1)] 
    nodes = pd.concat([links['gwasx'], links['gwasy']]).drop_duplicates().to_list()
    sign_index = pivot_corr_df['pval'][
        pivot_corr_df['pval'] < calculate_beta_correlation.get_pthres(corr_df)
    ].dropna(thresh=sign_threshold).index # only keep correlations significant in at least 2 methods
    links.set_index(['gwasx','gwasy'], inplace=True)
    sign_links = links[links.index.isin(sign_index)].reset_index()
    return sign_links, nodes

def get_gwas_links(pivot_corr_df, pattern, gwas_name):
    gwas_links = pivot_corr_df[(pivot_corr_df['gwasx']==gwas_name)
                               |(pivot_corr_df['gwasy']==gwas_name)].copy()
    gwas_links = gwas_links[(gwas_links['gwasx'].str.contains(pattern))|
                            (gwas_links['gwasy'].str.contains(pattern))]
    gwas_links['gwasx'],gwas_links['gwasy'] = np.where(gwas_links['gwasy']==gwas_name,
                                                       (gwas_links['gwasx'],gwas_links['gwasy']),
                                                       (gwas_links['gwasy'],gwas_links['gwasx'])
                                                      )
    gwas_links.sort_values('gwasx', inplace=True)
    return gwas_links

def chord_plot_pandas(gcircle, row, gwas_group_dict, bottom, cmap, norm):

    nodes = []
    for i in ['x','y']:
        node = row[f"gwas{i}"]
        for k,v in gwas_group_dict.items():
            if bool(re.search('|'.join(v),node)):
                node = f"{k}, {node}"
        nodes.append(node)
    gcircle.chord_plot([nodes[0],0,0,bottom],[nodes[1],0,0,bottom],color=cmap(norm(row['corr'])),alpha=.3)

    
def annotation_layer(gcircle, links, nodes, cmap, norm, bottom=1100, axis_color='k'):
    step = 0.001
    theta = np.linspace(0,(2-(2*step))*np.pi, len(nodes))
    scale = 500
    original_bottom = bottom
    width_scale = 5
    for m in constants.METHODS:
        corr_values = links[f'corr_{m}'] * scale #scales beta values so its visible

        gcircle.ax.bar(theta, corr_values, color=cmap(norm(corr_values)),
#                        edgecolor=axis_color,
                       edgecolor='none',
                       alpha=.5,
#                        alpha=1,
#                        width=step*np.pi*width_scale,
                       width = step*15,
                       bottom=bottom+scale, zorder=9,
                      )
        x = np.linspace(0,2*np.pi, len(nodes))
        y = np.array([bottom+scale]*len(x))
        gcircle.ax.plot(x, y, alpha=1, color=axis_color, lw=.5, zorder=8)
        gcircle.ax.plot(x, y-scale, alpha=.75, color=axis_color, lw=.25, zorder=8)
        gcircle.ax.plot(x, y+scale, alpha=.75, color=axis_color, lw=.25, zorder=8)
        gcircle.ax.plot(x, y-(0.5*scale), alpha=.75, color=axis_color, lw=.25, zorder=8)
        gcircle.ax.plot(x, y+(0.5*scale), alpha=.75, color=axis_color, lw=.25, zorder=8)        
        bottom += 2*scale
#         width_scale += 1
    # color the circosplot bars according to the number of significant correlations per method
    color_methods = ['purple','green','red'] # purple significant in 1 method, green in 2 and red in all 3
    links_m = links.copy()
    links_m['color'] = 'white' #'None'
    for i,c in zip(range(1,len(constants.METHODS)+1),color_methods):
        links_m.loc[links_m.loc[:,links_m.columns.str.contains('pval')].sum(1)==i,'color'] = c
    gcircle.ax.bar(theta, y-original_bottom+scale, color=links_m['color'],
                   edgecolor='none',#links_m['color'],
                   alpha=.25, width=step*np.pi*5, bottom=original_bottom, zorder=0)

def preprocessing(corr_df, gwas_group_dict, gwas_name, sign_threshold):
    pattern = [v for values in gwas_group_dict.values() for v in values]
    pattern = "|".join(pattern)
    pivot_corr_df = corr_df.pivot(index=['gwasx','gwasy'],columns='method')
    links, nodes = get_links_and_nodes(corr_df, pivot_corr_df, pattern, sign_threshold)
    
    # check significance of correlations
    pivot_corr_df.reset_index(inplace=True)
    pivot_corr_df['pval'] = pivot_corr_df['pval'] < calculate_beta_correlation.get_pthres(corr_df) # set True if signifiant
    pivot_corr_df.columns = [' '.join(col).strip().replace(' ','_' ) for col in pivot_corr_df.columns.values]
    
    gwas_links = get_gwas_links(pivot_corr_df, pattern, gwas_name)
    return links, nodes, gwas_links
    
def plot(corr_df, gwas_group_dict, gwas_name, sign_threshold=len(constants.METHODS)-2,
         bottom=1200, ylim=5000, color_map='tab10', figsize=(10,10),
         save=False, filename=None):
    '''
    corr_df: correlation dataframe from calculate_beta_correlation.calculate_celltype_corr
    gwas_group_dict: dictionary of GWAS groups to compare against
    gwas_name: name of the to be analyzed GWAS
    sign_threshold: the significance threshold for the inner chord plot
    '''
    gwas_group_dict = {k:v for k,v in gwas_group_dict.items() if gwas_name not in v}

    links, nodes, gwas_links = preprocessing(corr_df, gwas_group_dict, gwas_name, sign_threshold)
    
    plt.style.use('default')
    cmap = cm.get_cmap('seismic')
    norm = plt.Normalize(-1, 1)
    cmap_group = {group_name:cm.get_cmap(color_map)(i) for i,group_name in enumerate(gwas_group_dict.keys())}

    # make nodes
    gcircle = Gcircle()
    for node in nodes:
        for k,v in gwas_group_dict.items():
            if re.search("|".join(v),node):
                gcircle.add_locus(f"{k}, {node}", 2,
                                  bottom=bottom,
                                  linewidth=1, interspace=0, 
                                  facecolor=cmap_group[k],
                                  edgecolor=cmap_group[k])
    gcircle.set_locus(figsize=figsize) #Create figure object
    gcircle.ax.set_ylim(0,ylim)
    # make chords inside the circle
    links.apply(lambda row: chord_plot_pandas(gcircle, row, gwas_group_dict, bottom, cmap, norm), axis=1)
    # make the barplots around the circle
    annotation_layer(gcircle, gwas_links, nodes, cmap, norm, bottom=bottom, axis_color='k')
    # add legend
    patches = [mpatches.Patch(color=v, label=k) for k,v in cmap_group.items()]
    gcircle.ax.legend(handles=patches, bbox_to_anchor=(1.2, 1), loc=1, frameon=False)
    if save:
        plt.savefig(filename, dpi=150, bbox_inches='tight')