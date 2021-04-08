# to allow importing to work correctly (in a dirty way)
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import constants
import pandas as pd
import mygene
import gseapy
from pathlib import Path

def get_annot_list(df, name, param_list, rank):
    param = '|'.join(param_list)
    df_gsea = pd.concat([df[(df['gwas'].str.contains(param))][['gwas','specificity_id','annotation']],
                       (df[(df['gwas'].str.contains(param))][f'pvalue_{constants.PVAL_CORRECTION}']<=0.05).astype(int)
                       ], axis=1).rename(columns={f'pvalue_{constants.PVAL_CORRECTION}':'count'}).groupby(['gwas','specificity_id','annotation']).sum()
    n_gwas = df_gsea.index.unique(level='gwas').shape[0]
    # only use count==1 if multiple gwases are used in analysis
    if n_gwas == 1:
        min_count = 0 # if only one gwas in gwas_group than include all cell-types
    else:
        min_count = 1 # if multiple gwas, cell type should be included in at least 2
        min_count = 0
    original_shape = df[(df['gwas'].str.contains(param))].shape
    df_gsea = df_gsea[df_gsea['count']>=2] # only get if significant in >2 methods
    df_gsea['count'] = 1 # reset count to 1
    df_gsea = df_gsea.groupby(['specificity_id','annotation']).sum().reset_index()
    N_total = df[(df['gwas'].str.contains(param))]['gwas'].unique().shape[0]
    df_gsea['freq'] = df_gsea['count'] / N_total
    df_gsea.sort_values('freq', ascending=False, inplace=True)
    df_gsea['rank'] = df_gsea['count'].rank(ascending=False,method='dense').astype(int)
    if rank is None:
        df_gsea = df_gsea[(df_gsea['count']>min_count)]
    else:
        df_gsea = df_gsea[(df_gsea['rank']<=rank)&
                          (df_gsea['count']>min_count)
                         ]
    annot_list = df_gsea[['specificity_id','annotation']].values.tolist()
    print(f"Top {rank} ranked cell-types (N={len(annot_list)}) in {name} GWAS:")
    for row in df_gsea.iterrows():
        a = row[1]
#         print(f"{a[4]}. {a[0]}: {a[1]} (N={a[2]}, freq={a[3]:.2})")
        print(f"{a[4]}. {a[0]}: {a[1]} (N={a[2]})")
    print()
#     top_cell_df = df_gsea[['specificity_id','annotation','freq','rank']].copy()
#     top_cell_df['cell-type'] = (df_gsea['specificity_id']+': '+df_gsea['annotation'])
#     top_cell_df.drop(columns=['specificity_id','annotation'], inplace=True)
#     top_cell_df['group'] = name
#     df_list.append(top_cell_df)
    return annot_list


def get_top_genes(annot_list):
    celltype_genes_dict = {}
    esmu_dir = 'esmu'
    for dataset, celltype in annot_list:
        cellex_file = f'{esmu_dir}/{dataset}.mu.csv' # change esmu to mu if file not found
        if Path(cellex_file).is_file():
            df_esmu = pd.read_csv(cellex_file, index_col=0)
        elif Path(cellex_file.replace('.mu','.esmu')).is_file():
            df_esmu = pd.read_csv(cellex_file.replace('.mu','.esmu'), index_col=0)
        else:
            print('file not found')
        df = df_esmu[celltype]
        df = df.sort_values(ascending=False).head(round(constants.TOP_FREQ*len(df)))
        if df.shape[0] > 500 or df.shape[0] < 15: #https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html
            print(f'\033[93m{df.shape[0]} genes in {dataset}, {celltype} which is outside the recommended range (15<G<500)\033[0m')
        celltype_genes_dict[f'{dataset}, {celltype}'] = df.index.to_list()

    print('\nConverting Ensembl ID to Gene Symbol...')
    for i, (celltype, genes) in enumerate(celltype_genes_dict.items(),1):
        gsea_dir = 'gsea'
        out_file = f"{gsea_dir}/gsea_{celltype.replace(', ','-')}.xlsx"
        if Path(out_file).is_file() and not constants.OVERWRITE_GSEA_ANALYSIS:
            print(f'\n({i}/{len(celltype_genes_dict)}) Converting {celltype} ({len(genes)} genes)')
            print(f'Cell-type already analyzed. Skipping conversion...')
        else:
            print(f'\n({i}/{len(celltype_genes_dict)}) Converting {celltype} ({len(genes)} genes)')
            mg = mygene.MyGeneInfo()
            ginfo = mg.querymany(genes, scopes='ensembl.gene')
            gene_list = [g['symbol'] for g in ginfo if 'symbol' in g]
            celltype_genes_dict[celltype] = gene_list
    return celltype_genes_dict

def summarize_gsea(gsea_dict, correct_pval=True, min_count=0, save_to_excel=False, filename='summarize_gsea.xlsx'):
    total = len(gsea_dict)
    df_list = []
    for k,v in gsea_dict.items():
        df = v.copy()
        df['Celltype'] = k
        df_list.append(df)
    gsea_df = pd.concat(df_list)
    if correct_pval:
        gsea_df = gsea_df[(gsea_df['Adjusted P-value']<=0.05/total)] # Bonferroni correction
    gsea_grouped_df = gsea_df.groupby(['Gene_set','Term'])['Celltype'].agg(list).reset_index()
    gsea_grouped_df[f'Celltype_count (total={total})'] = gsea_grouped_df['Celltype'].apply(lambda x: len(x))
    gsea_grouped_df.sort_values(f'Celltype_count (total={total})', ascending=False, inplace=True)
    gsea_grouped_df['Celltype'] = gsea_grouped_df['Celltype'].apply(lambda x : '; '.join(x))
    gsea_grouped_df = gsea_grouped_df[gsea_grouped_df[f'Celltype_count (total={total})']>min_count]
    if save_to_excel:         
        gsea_grouped_df.to_excel(filename, index=False)
    return gsea_grouped_df

def gsea(df, gwas_group_dict, rank=constants.TOP_ANNOT):
    '''
    GSEA on enriched cell types.
    rank: The rank a cell type should have to be analysed. 
    For example rank=5 means that only cell types ranked in the top 5 most occuring enriched cell types will be analysed.
    If rank=None, then all cell types will be analysed.
    '''
    print(f"Performing GSEA...\n")
    for name, param_list in gwas_group_dict.items():
        annot_list = get_annot_list(df, name, param_list, rank)
        celltype_genes_dict = get_top_genes(annot_list)
 
        print('\nRunning Enrichr...')
        gsea_dict = {}
        gsea_dir = 'gsea'
        for i, (celltype, genes) in enumerate(celltype_genes_dict.items(),1):
            out_file = f"{gsea_dir}/gsea_{celltype.replace(', ','-')}.xlsx"
            print(f'\n({i}/{len(celltype_genes_dict)}) Analyzing {celltype} ({len(genes)} genes)...')
            if Path(out_file).is_file() and not constants.OVERWRITE_GSEA_ANALYSIS:
                print('Cell-type already analyzed. Skipping analysis...')
                gsea_dict[celltype] = pd.read_excel(out_file)
            else:
                df_list = []
                for gene_set in constants.GENE_SET_LIST:
#                     print(f'Gene-set: {gene_set}')
                    try:
                        enr = gseapy.enrichr(gene_list=genes,
                                             gene_sets=gene_set,
                                             outdir=None,
                                             organism='Human',
                                            )
                        df = enr.results[enr.results['Adjusted P-value']<=0.05]
                        df_list.append(df)
                    except Exception as e:
                        print(e)
                        continue
                try:
                    gsea_dict[celltype] = pd.concat(df_list)
                    gsea_dict[celltype].to_excel(out_file, index=False)
                    print('Writing to file...')
                except ValueError:
                    continue
    return gsea_dict

if __name__ == "__main__":
    df_all = pd.read_hdf('data/CELLECT_output/data.h5', 'df')
    gsea(df_all, constants.GWAS_GROUP_DICT)
