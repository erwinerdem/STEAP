# to allow importing to work correctly (in a dirty way)
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import constants
import pandas as pd
import itertools
from scipy.stats import pearsonr

def calculate_pearson(dataframe):
    corr_list = []
    for x,y in itertools.combinations(dataframe.columns, 2):
        corr, pval = pearsonr(dataframe.loc[:,x].values, dataframe.loc[:,y].values)
        corr_list.append([x,y,corr,pval])
    corr_df = pd.DataFrame(corr_list, columns=['gwasx','gwasy','corr','pval'])
    return corr_df

def get_pthres(corr_df):
    return corr_df[corr_df['corr']<0]['pval'].min()

def calculate_celltype_corr(dataframe):
    df_list = []
    for m in constants.METHODS:
        df_method = dataframe[(dataframe.method==m)]
        df_method_pivot = df_method.pivot_table(index='annotation', columns='gwas', values='beta')
        df_method_pivot.dropna(axis='columns',inplace=True) #drop gwas if not complete
        corr_df = calculate_pearson(df_method_pivot)
        corr_df['method'] = m
        df_list.append(corr_df)
    corr_df = pd.concat(df_list)
    return corr_df

    
if __name__ == "__main__":
    df_all = pd.read_hdf('data/CELLECT_output/data.h5', 'df_all')
    corr_df = calculate_celltype_corr(df_all)
    corr_df.to_hdf('data/CELLECT_output/data.h5', key='corr_df')
