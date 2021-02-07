from constants import *
import pandas as pd
import mygene
import gseapy
from pathlib import Path

df_all = pd.read_hdf('../data/CELLECT_output/data.h5', 'df_all')

print(f"Performing GSEA...\n")

df_list = []
for name, param_list in GWAS_GROUP_DICT.items():
    param = '|'.join(param_list)

    df_gsea = pd.concat([df_all[(df_all['gwas'].str.contains(param))][['gwas','specificity_id','annotation']],
                       (df_all[(df_all['gwas'].str.contains(param))]['pvalue_bonferroni'] <= 0.05).astype(int)
                       ],
                      axis=1
                     ).rename(columns={'pvalue_bonferroni':'count'}).groupby(['gwas','specificity_id','annotation']).sum()
    df_gsea = df_gsea[df_gsea['count']>=2]
    df_gsea['count'] = 1
    df_gsea = df_gsea.groupby(['specificity_id','annotation']).sum().reset_index()
    N_total = df_all[(df_all['gwas'].str.contains(param))]['gwas'].unique().shape[0] 
    df_gsea['freq'] = df_gsea['count'] / N_total
    df_gsea.sort_values('freq', ascending=False, inplace=True)
    df_gsea['rank'] = df_gsea['count'].rank(ascending=False,method='dense').astype(int)

    df_gsea = df_gsea[(df_gsea['rank']<=TOP_ANNOT)&
                      (df_gsea['count']>1)
                     ]
    annot_list = df_gsea[['specificity_id','annotation']].values.tolist()
    datasets = df_all['specificity_id'].unique()
    print(f"Top {TOP_ANNOT} ranked cell-types (N={len(annot_list)}) in {name} GWAS:")
    for row in df_gsea.iterrows():
        a = row[1]
        print(f"{a[4]}. {a[0]}: {a[1]} (N={a[2]}, freq={a[3]:.2})")
    print()
    top_cell_df = df_gsea[['specificity_id','annotation','freq','rank']].copy()
    top_cell_df['cell-type'] = (df_gsea['specificity_id']+': '+df_gsea['annotation'])
    top_cell_df.drop(columns=['specificity_id','annotation'], inplace=True)
    top_cell_df['group'] = name
    df_list.append(top_cell_df)

    celltype_genes_dict = {}
    for dataset, celltype in annot_list:
        cellex_file = f'../esmu/{dataset}.mu.csv' # change esmu to mu if file not found
        if Path(cellex_file).is_file():
            df_esmu = pd.read_csv(cellex_file, index_col=0)
        elif Path(cellex_file.replace('.mu','.esmu')).is_file():
            df_esmu = pd.read_csv(cellex_file.replace('.mu','.esmu'), index_col=0)
        else:
            print('file not found')
        df = df_esmu[celltype]
        df = df.sort_values(ascending=False).head(round(TOP_FREQ*len(df)))
        if df.shape[0] > 500 or df.shape[0] < 15: #https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html
            print(f'\033[93m{df.shape[0]} genes in {dataset}, {celltype} which is outside the recommended range (15<G<500)\033[0m')
        celltype_genes_dict[f'{dataset}, {celltype}'] = df.index.to_list()

    print('\nConverting Ensembl ID to Gene Symbol...')
    for i, (celltype, genes) in enumerate(celltype_genes_dict.items(),1):
        print(f'\n({i}/{len(celltype_genes_dict)}) Converting {celltype} ({len(genes)} genes)')
        mg = mygene.MyGeneInfo()
        ginfo = mg.querymany(genes, scopes='ensembl.gene')
        gene_list = [g['symbol'] for g in ginfo if 'symbol' in g]
        celltype_genes_dict[celltype] = gene_list

    print('\nRunning Enrichr...')
    gsea_dict = {}
    for i, (celltype, genes) in enumerate(celltype_genes_dict.items(),1):
        out_file = f"{GSEA_OUTDIR}gsea_{celltype.replace(', ','-')}.xlsx"
        print(f'\n({i}/{len(celltype_genes_dict)}) Analyzing {celltype} ({len(genes)} genes)...')
        if Path(out_file).is_file() and not OVERWRITE_GSEA_ANALYSIS:
            print('Cell-type already analyzed. Skipping analysis...')
            gsea_dict[celltype] = pd.read_excel(out_file)
        else:
            df_list = []
            for gene_set in GENE_SET_LIST:
                print(f'Gene-set: {gene_set}')
                try:
                    enr = gseapy.enrichr(gene_list=genes,
                                         gene_sets=gene_set,
                                         outdir=None,
					 organism='Human',
                                        )
                    df = enr.results[enr.results['Adjusted P-value']<0.05]
                    df_list.append(df)
                except:
                    print('No enrichment found...')
            try:
                gsea_dict[celltype] = pd.concat(df_list)
                gsea_dict[celltype].to_excel(out_file, index=False)
            except ValueError:
                continue
