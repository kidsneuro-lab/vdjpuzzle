import pandas as pd
import os
# Have to be in the directory containing the chung results

if not os.path.exists('Unique_results'):
    os.makedirs('Unique_results')

for file in os.listdir():
    if file in ['TRA.tsv', 'TRB.tsv', 'TRD.tsv', 'TRG.tsv', 'IGH.tsv', 'IGK.tsv', 'IGL.tsv']:
        df = pd.read_csv(file, sep='\t')
        df = df.sort_values('Expression', ascending = False)
        if file in ['IGH.tsv', 'IGK.tsv', 'IGL.tsv']:
            df = df.dropna(subset=['Isotype/Constant'])
            df = df.drop_duplicates(subset=['CellID', 'cdr3aa', 'Isotype/Constant'], keep='first')
            df = (df.groupby(['CellID', 'cdr3aa'], as_index=False).apply(lambda x: x.loc[(x['Isotype/Constant'] == 'IGHM') | (x['Isotype/Constant'] == 'IGHD')] if (x.iloc[[0]]['Isotype/Constant'] == 'IGHM').bool() else
                                                             (x.loc[(x['Isotype/Constant'] == 'IGHM') | (x['Isotype/Constant'] == 'IGHD')] if (x.iloc[[0]]['Isotype/Constant'] == 'IGHD').bool() else
                                                              (x.iloc[[0]])
                                                             )
                                                       )
            )
        else:
            df = df.drop_duplicates(subset=['CellID', 'cdr3aa'], keep='first')
        df = df.sort_values(['CellID', 'cdr3aa', 'Expression'] , ascending = [True, True, False])
        df.to_csv("Unique_results/" + file, sep='\t', index=False)