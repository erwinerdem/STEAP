import constants
import pandas as pd
from pathlib import Path
from statsmodels.stats.multitest import multipletests


def csv_file(directory):
    file_dict = {}
    for path in Path(directory).rglob('prioritization.csv'):
        full_path = str(path)
        (name, method,__, __) = path.parts[-4:]
        name = name[8:] # CELLECT- is 8 long
        method = method[8:]
        if name not in file_dict:
            file_dict[name] = {}
        file_dict[name].update({method:full_path})
    return file_dict


def make_df(directory):
    '''Make a DataFrama for the .csv files.'''
    file_dict = csv_file(directory)
    df_list_1 = []
    for name,d in file_dict.items():
        df_list_2 = []
        for method, file_path in d.items():
            df = pd.read_csv(file_path)
            df['method'] = method
            df.sort_values(by=['gwas','specificity_id','annotation'], inplace=True)
            df_list_2.append(df)
        df_list_1.extend(df_list_2)
    df_all = pd.concat(df_list_1, ignore_index=True)

    # count the number of methods used (MAGMA/H-MAGMA/LDSC)
    df_all = df_all.merge(df_all.groupby(['gwas','specificity_id','annotation']).size()\
                            .to_frame('n_methods'), on=['gwas','specificity_id','annotation'], how='left')
    # count the number of annotations/celltypes
    df_all.sort_values(by=['gwas','method'], inplace=True)
    df_all.reset_index(inplace=True, drop=True)
    return df_all


def pvalue_correction(dataframe, method='bonferroni'):
    '''
    Available methods can be found at
    https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html
    '''
    df_p = dataframe.pivot_table(values='pvalue', index=['method','gwas','specificity_id',],
                                     columns=['annotation'])
    df_p = df_p.apply(lambda row: multipletests(row.dropna(), method=method)[1],
                   axis=1,
                   result_type='reduce'
                  ).apply(pd.Series).stack().reset_index().drop('level_3', axis=1)\
    .rename(columns={0:f"pvalue_{method}"})
    df_p['annotation'] = dataframe['annotation']
    return pd.merge(dataframe, df_p, on=['gwas','specificity_id','annotation','method'])

if __name__ == "__main__":
    df_all = make_df(constants.CELLECT_OUTDIR)
    df_all = pvalue_correction(df_all, method=constants.PVAL_CORRECTION)
    df_all.to_hdf('data/CELLECT_output/data.h5', key='df', mode='w')

