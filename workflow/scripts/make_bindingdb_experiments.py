import argparse
from pathlib import Path

import pandas as pd


def normalize_cid(values: pd.Series) -> pd.Series:
	return pd.to_numeric(values, errors='coerce').astype('Int64')


def parse_affinity(values: pd.Series) -> pd.Series:
	text = values.astype('string').str.strip().str.replace(',', '', regex=False)
	text = text.mask(text.str.contains(r'[<>~]', na=False))
	return pd.to_numeric(text, errors='coerce')


def load_hdd_cid_map(coldata_path: str) -> pd.DataFrame:
	coldata = pd.read_csv(
		coldata_path,
		usecols=['HDD.Compound.ID', 'Pubchem.CID'],
		dtype='string',
	)
	coldata['PubChem CID'] = normalize_cid(coldata['Pubchem.CID'])
	coldata = coldata.dropna(subset=['PubChem CID'])
	coldata = coldata[['HDD.Compound.ID', 'PubChem CID']]
	coldata = coldata.drop_duplicates()
	return coldata


def make_bindingdb_matrix(
	coldata_path: str,
	bindingdb_path: str,
	output_path: str,
) -> None:
	hdd_cid_map = load_hdd_cid_map(coldata_path)
	bindingdb = pd.read_csv(bindingdb_path, dtype='string')
	bindingdb['PubChem CID'] = normalize_cid(bindingdb['PubChem CID'])
	bindingdb['Target Name'] = bindingdb['Target Name'].astype('string').str.strip()
	bindingdb['Ki.numeric'] = parse_affinity(bindingdb['Ki (nM)'])
	bindingdb['Kd.numeric'] = parse_affinity(bindingdb['Kd (nM)'])
	bindingdb['BindingDB.Affinity.nM'] = bindingdb[['Ki.numeric', 'Kd.numeric']].min(
		axis=1, skipna=True
	)
	bindingdb = bindingdb.dropna(
		subset=['PubChem CID', 'Target Name', 'BindingDB.Affinity.nM']
	)

	matched = bindingdb.merge(hdd_cid_map, on='PubChem CID', how='inner')
	if matched.empty:
		matrix = pd.DataFrame(columns=['BindingDB.Target.Name'])
	else:
		target_compound = (
			matched.groupby(['Target Name', 'HDD.Compound.ID'], sort=True)[
				'BindingDB.Affinity.nM'
			]
			.min()
			.reset_index()
		)
		matrix = target_compound.pivot_table(
			index='Target Name',
			columns='HDD.Compound.ID',
			values='BindingDB.Affinity.nM',
			aggfunc='min',
		)
		ordered_columns = [
			hdd_id
			for hdd_id in hdd_cid_map['HDD.Compound.ID']
			if hdd_id in matrix.columns
		]
		matrix = matrix.loc[sorted(matrix.index), ordered_columns]
		matrix = matrix.reset_index(names='BindingDB.Target.Name')

	output = Path(output_path)
	output.parent.mkdir(parents=True, exist_ok=True)
	matrix.to_csv(output, index=False)
	print(  # noqa: T201
		'[make_bindingdb_experiments] '
		f'targets={max(len(matrix), 0)} '
		f'compounds={max(matrix.shape[1] - 1, 0)} output={output}',
		flush=True,
	)


def main_from_snakemake() -> None:
	make_bindingdb_matrix(
		coldata_path=str(snakemake.input.colData),
		bindingdb_path=str(snakemake.input.bindingdb),
		output_path=str(snakemake.output.binding_db),
	)


if __name__ == '__main__':
	if 'snakemake' in globals():
		main_from_snakemake()
	else:
		parser = argparse.ArgumentParser(
			prog='make_bindingdb_experiments',
			description='Generate the HDD BindingDB target affinity matrix',
		)
		parser.add_argument('--coldata', required=True, help='HDD colData CSV path')
		parser.add_argument(
			'--bindingdb', required=True, help='Cleaned BindingDB CSV path'
		)
		parser.add_argument('--output', required=True, help='Output BindingDB CSV path')
		args = parser.parse_args()

		make_bindingdb_matrix(
			coldata_path=args.coldata,
			bindingdb_path=args.bindingdb,
			output_path=args.output,
		)
