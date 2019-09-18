import pandas as pd
import os
import sys

if not os.path.exists('final_receptor_results'):
	os.makedirs('final_receptor_results')

for file in os.listdir(sys.argv[1]):
	if file in ['TRA.tsv', 'TRB.tsv', 'TRD.tsv', 'TRG.tsv', 'IGH.tsv', 'IGK.tsv', 'IGL.tsv']:
		print(file)
		if os.stat(sys.argv[1] + '/' + file).st_size == 0:
			print("{} is empty, ignoring".format(file))
			continue

		df = pd.read_csv(sys.argv[1] + '/' + file, sep='\t')
		if len(df.index) == 0 or 'Expression' not in df.columns:
			print("{} has no data, ignoring".format(file))
			continue

		df = df.sort_values('Expression', ascending = False)
		if file in ['IGH.tsv', 'IGK.tsv', 'IGL.tsv']:
			# Remove BCR sequences without an isotype/constant region
			df = df.dropna(subset=['Isotype/Constant'])
			# Remove duplicate sequences with identical isotype/constant region (keeping only the one with the highest expression)
			df = df.drop_duplicates(subset=['CellID', 'cdr3_aa', 'Isotype/Constant'], keep='first')
			# If the top expressed sequence for a cell is IGHM, keep the second highest sequence for that cell that is IGHD
			# If the top expressed sequence for a cell is IGHD, keep the second highest sequence for that cell that is IGHM
			# If the top expressed sequences is not IGHM/IGHD, remove all other sequences for that cell
			df = (df.groupby(['CellID', 'cdr3_aa'], as_index=False).apply(lambda x: x.loc[(x['Isotype/Constant'] == 'IGHM') | (x['Isotype/Constant'] == 'IGHD')] if (x.iloc[[0]]['Isotype/Constant'] == 'IGHM').bool() else
															 (x.loc[(x['Isotype/Constant'] == 'IGHM') | (x['Isotype/Constant'] == 'IGHD')] if (x.iloc[[0]]['Isotype/Constant'] == 'IGHD').bool() else
															  (x.iloc[[0]])
															  )
															 )
			)
		else:
			df = df.drop_duplicates(subset=['CellID', 'cdr3_aa'], keep='first')
			df = df.sort_values(['CellID', 'cdr3_aa', 'Expression'] , ascending = [True, True, False])

		df.to_csv("final_receptor_results/" + file, sep='\t', index=False)
		print("{} successfully filtered".format(file))
