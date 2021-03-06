# to allow importing to work correctly (in a dirty way)
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import constants
import pandas as pd
import itertools
from scipy.stats import kendalltau
from pathlib import Path

def calculate_kendalltau(dataframe):
    corr_list = []
    for x,y in itertools.combinations(dataframe.columns, 2):
        corr_frame = dataframe.loc[:,[x,y]].fillna(0).copy()
        corr_frame = corr_frame[(corr_frame>0).all(1)]
        corr, pval = kendalltau(corr_frame.iloc[:,0].values, corr_frame.iloc[:,1].values)
        corr_list.append([x,y,corr, pval])
    corr_df = pd.DataFrame(corr_list, columns=['celltypex','celltypey','corr','pval'])
    return corr_df

def correct_pval_correlation(corr_df):
    '''Uses Bonferroni correction for the input correlation pvalues.'''
    corrected_df = corr_df.copy()
    n_test = corrected_df.shape[0]
    corrected_df['pval_bonferroni'] = corrected_df['pval'] * n_test
    corrected_df.loc[corrected_df['pval_bonferroni'] > 1, 'pval_bonferroni'] = 1
    return corrected_df

    
def calculate_es_corr(datasets):
    df_list = []
    for dataset in datasets:
        cellex_file = f'esmu/{dataset}.mu.csv' # change esmu to mu if file not found
        if Path(cellex_file).is_file():
            df_esmu = pd.read_csv(cellex_file, index_col=0)
        elif Path(cellex_file.replace('.mu','.esmu')).is_file():
            df_esmu = pd.read_csv(cellex_file.replace('.mu','.esmu'), index_col=0)
        else:
            print('file not found')
        df_esmu.columns = [f'{dataset}, {ct}' for ct in df_esmu.columns]
        df_list.append(df_esmu)
    merged_es_df = pd.concat(df_list, join='outer', axis=1)
    merged_es_df.sort_index(axis=1, inplace=True)
    es_corr_df = calculate_kendalltau(merged_es_df.fillna(0))
    es_corr_df = correct_pval_correlation(es_corr_df)
    return es_corr_df

    
if __name__ == "__main__":
    df_all = pd.read_hdf('data/CELLECT_output/data.h5', 'df_all')
    datasets = df_all['specificity_id'].unique()
    es_corr_df = calculate_es_corr(datasets)
    es_corr_df.to_hdf('data/CELLECT_output/data.h5', key='es_corr_df')
