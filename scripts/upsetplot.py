# to allow importing to work correctly (in a dirty way)
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import constants
import pandas as pd
import re
import upsetplot
import matplotlib.pyplot as plt


def create_group(gwas, gwas_group_dict):
    reverse_dict = {v:k for k,v_list in gwas_group_dict.items() for v in v_list}
    for k,v in reverse_dict.items():
        if bool(re.search(k,gwas)):
            return v

def add_count_to_group_dict(df, gwas_group_dict):
    unqiue_gwas = df['gwas'].unique().tolist()
    reverse_dict = {v:k for k,v_list in gwas_group_dict.items() for v in v_list}
    gwas_count = {k:0 for k in reverse_dict.keys()}
    for gwas in unqiue_gwas:
        for k,v in reverse_dict.items():
            if bool(re.search(k,gwas)):
                gwas_count[k] +=1
    gwas_group_dict_copy = {}
    for k, v_list in gwas_group_dict.items():
        tot_count = (sum([gwas_count[v] for v in v_list]))
        gwas_group_dict_copy[f'{k} (N={tot_count})'] = v_list
    return gwas_group_dict_copy

def make_df_sliced(df, gwas_group_dict, sign_threshold):
    pattern = '|'.join([v for v_list in gwas_group_dict.values() for v in v_list])
    df_sliced = df[(df['gwas'].str.contains(pattern))].copy()
    df_sliced = df_sliced[df_sliced[f'pvalue_{constants.PVAL_CORRECTION}']<=0.05]
    df_sliced = df_sliced.groupby(['gwas','specificity_id','annotation']).size().reset_index().rename(
        columns={0:'N_sign'})
    df_sliced = df_sliced[df_sliced['N_sign']>=sign_threshold] # significant in >2 methods

    df_sliced['group'] = df_sliced.apply(lambda x : create_group(x['gwas'],gwas_group_dict),
                                        axis=1)
    return df_sliced


def make_df_upset(df, gwas_group_dict, sign_threshold, sort_categories_by):
    df_sliced = make_df_sliced(df, gwas_group_dict, sign_threshold)
    df_sliced_upset = df_sliced.groupby(['specificity_id','annotation'])['group'].agg(list).reset_index()
    df_sliced_upset['group'] = df_sliced_upset['group'].apply(set)
    df_sliced_upset = df_sliced_upset['group'].value_counts().reset_index()
#     groups = list(df_sliced['group'].unique())
    groups = list(gwas_group_dict.keys())[::-1]
    for g in groups:
        df_sliced_upset[g] = df_sliced_upset['index'].apply(lambda x : g in x)
    df_sliced_upset.drop(columns='index', inplace=True)
    if sort_categories_by is None:
        new_column = groups.copy()
        new_column.append('group')
        df_sliced_upset = df_sliced_upset.reindex(columns=new_column)
    df_sliced_upset.set_index(groups, inplace=True)
    return df_sliced_upset 


def get_shared_celltypes(df, gwas_group_dict, sign_threshold,
                         save_to_excel=False, filename=None):
    df_sliced = make_df_sliced(df, gwas_group_dict, sign_threshold)
    df_sliced = df_sliced.groupby(['specificity_id','annotation'])[['group','gwas']].agg(set).reset_index()
    df_sliced['gwas'] = df_sliced['gwas'].apply(lambda x : ', '.join(list(x)))
    df_sliced['N_groups'] = df_sliced['group'].apply(len)
    df_sliced['group'] = df_sliced['group'].apply(lambda x : ', '.join(list(x)))
    df_sliced.sort_values(['N_groups','group'], ascending=[False,True], inplace=True)
    df_sliced.drop(columns='N_groups', inplace=True)
    if save_to_excel:         
        df_sliced.to_excel(filename, index=False)
    return df_sliced

def plot_upset(df, gwas_group_dict, sign_threshold=len(constants.METHODS)-1,
               save=False, filename=None,
               show_percentages=False, sort_by='cardinality',
               with_lines=True,  element_size=46, show_counts='%d', sort_categories_by = 'cardinality'):
    gwas_group_dict_copy = add_count_to_group_dict(df, gwas_group_dict)
    df_sliced_upset = make_df_upset(df, gwas_group_dict_copy, sign_threshold, sort_categories_by)
    plt.style.use('default')
    upsetplot.plot(df_sliced_upset['group'], show_percentages=show_percentages,
                   sort_by=sort_by, with_lines=with_lines, #cardinality, degree
                   element_size=element_size, 
                   show_counts=show_counts,
                   sort_categories_by=sort_categories_by #cardinality, None
                  )
    if save:
        plt.savefig(filename, dpi=200, bbox_inches='tight')
    plt.show()



if __name__ == "__main__":
    df_all = pd.read_hdf('data/CELLECT_output/data.h5', 'df_all')
    plot_upset(df_all, constants.GWAS_GROUP_DICT, save=True, 
               filename='figures_and_tables/upsetplot.png')
    get_shared_celltypes(df_all, constants.GWAS_GROUP_DICT, save_to_excel=True, 
                         filename='figures_and_tables/upsetplot.xslx')