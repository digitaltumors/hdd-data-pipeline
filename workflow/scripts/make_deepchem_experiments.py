import argparse
from pathlib import Path

import pandas as pd


def load_smiles_to_cid_map(mapping_path: str) -> pd.DataFrame:
	mapping = pd.read_csv(
		mapping_path,
		sep='\t',
		header=None,
		names=['smiles', 'pubchem_cid'],
		compression='gzip',
		dtype={'smiles': 'string', 'pubchem_cid': 'string'},
		keep_default_na=False,
	)
	mapping['smiles'] = mapping['smiles'].astype(str)
	mapping['pubchem_cid'] = mapping['pubchem_cid'].replace('', pd.NA)
	mapping['Pubchem.CID'] = pd.to_numeric(
		mapping['pubchem_cid'],
		errors='coerce',
	).astype('Int64')
	mapping = mapping.dropna(subset=['Pubchem.CID'])
	mapping = mapping.drop_duplicates(subset=['smiles'], keep='first')
	return mapping[['smiles', 'Pubchem.CID']]


def build_assay_matrix(
	data: pd.DataFrame,
	measurement_cols: list[str],
	index_name: str,
) -> pd.DataFrame:
	data = data[['HDD.Compound.ID'] + measurement_cols].transpose()
	data = data.rename(columns=data.iloc[0])
	data = data.iloc[1:,]
	data = data.reset_index(names=index_name, drop=False)
	return data


def collapse_measurements(
	data: pd.DataFrame,
	measurement_cols: list[str],
) -> pd.DataFrame:
	return (
		data.groupby('HDD.Compound.ID', sort=False)[measurement_cols]
		.first()
		.reset_index()
	)


def process_deepchem_data(
	coldata: pd.DataFrame,
	smiles_to_cid: pd.DataFrame,
	data: pd.DataFrame,
	index_name: str,
	convert_to_int: bool = False,
) -> pd.DataFrame:
	measurement_cols = [
		col_name for col_name in data.columns if col_name not in ['smiles', 'mol_id']
	]
	data = data.copy()
	data['smiles'] = data['smiles'].astype(str)

	if convert_to_int:
		data[measurement_cols] = (
			data[measurement_cols]
			.apply(lambda col: pd.to_numeric(col, errors='coerce'))
			.astype('Int64')
		)

	coldata = coldata.copy()
	coldata['SMILES'] = coldata['SMILES'].astype(str)
	coldata['Pubchem.CID'] = pd.to_numeric(
		coldata['Pubchem.CID'],
		errors='coerce',
	).astype('Int64')
	cid_lookup = coldata[['HDD.Compound.ID', 'Pubchem.CID']].dropna(
		subset=['Pubchem.CID']
	)

	exact_matches = data.merge(
		coldata[['SMILES', 'HDD.Compound.ID']],
		left_on='smiles',
		right_on='SMILES',
		how='inner',
	)
	unmatched_data = data[~data['smiles'].isin(set(coldata['SMILES']))].copy()
	cid_matches = unmatched_data.merge(
		smiles_to_cid,
		on='smiles',
		how='inner',
	)
	cid_matches = cid_matches.merge(
		cid_lookup,
		on='Pubchem.CID',
		how='inner',
	)

	matched = pd.concat(
		[
			exact_matches[['HDD.Compound.ID'] + measurement_cols],
			cid_matches[['HDD.Compound.ID'] + measurement_cols],
		],
		ignore_index=True,
	)
	matched = collapse_measurements(matched, measurement_cols)
	return build_assay_matrix(matched, measurement_cols, index_name)


def main(
	coldata_path: str,
	smiles_to_cid_path: str,
	clintox_input: str,
	tox21_input: str,
	toxcast_input: str,
	sider_input: str,
	clintox_output: str,
	tox21_output: str,
	toxcast_output: str,
	sider_output: str,
) -> None:
	coldata = pd.read_csv(
		coldata_path,
		usecols=['HDD.Compound.ID', 'SMILES', 'Pubchem.CID'],
	)
	smiles_to_cid = load_smiles_to_cid_map(smiles_to_cid_path)
	clintox = pd.read_csv(clintox_input)
	tox21 = pd.read_csv(tox21_input)
	toxcast = pd.read_csv(toxcast_input)
	sider = pd.read_csv(sider_input)
	Path(clintox_output).parent.mkdir(parents=True, exist_ok=True)

	clintox = process_deepchem_data(
		coldata,
		smiles_to_cid,
		clintox,
		'Clinical Tox Result',
	)
	sider = process_deepchem_data(
		coldata,
		smiles_to_cid,
		sider,
		'Side Effect',
	)
	tox21 = process_deepchem_data(
		coldata,
		smiles_to_cid,
		tox21,
		'Tox Assay',
		True,
	)
	toxcast = process_deepchem_data(
		coldata,
		smiles_to_cid,
		toxcast,
		'Tox Assay',
	)

	clintox.to_csv(clintox_output, index=False)
	sider.to_csv(sider_output, index=False)
	tox21.to_csv(tox21_output, index=False)
	toxcast.to_csv(toxcast_output, index=False)


def main_from_snakemake() -> None:
	main(
		coldata_path=str(snakemake.input.colData),
		smiles_to_cid_path=str(snakemake.input.smiles_to_cid),
		clintox_input=str(snakemake.input.clintox),
		tox21_input=str(snakemake.input.tox21),
		toxcast_input=str(snakemake.input.toxcast),
		sider_input=str(snakemake.input.sider),
		clintox_output=str(snakemake.output.clintox),
		tox21_output=str(snakemake.output.tox21),
		toxcast_output=str(snakemake.output.toxcast),
		sider_output=str(snakemake.output.sider),
	)


if __name__ == '__main__':
	if 'snakemake' in globals():
		main_from_snakemake()
	else:
		parser = argparse.ArgumentParser(
			prog='make_deepchem_experiments',
			description='Generate DeepChem experiments matrix',
		)
		parser.add_argument('-c', '--coldata', required=True, help='colData CSV path')
		parser.add_argument(
			'--smiles-to-cid',
			required=True,
			help='PubChem SMILES-to-CID mapping file path',
		)
		parser.add_argument('--clintox-input', required=True, help='Input ClinTox CSV')
		parser.add_argument('--tox21-input', required=True, help='Input Tox21 CSV')
		parser.add_argument('--toxcast-input', required=True, help='Input ToxCast CSV')
		parser.add_argument('--sider-input', required=True, help='Input SIDER CSV')
		parser.add_argument(
			'--clintox-output', required=True, help='Output ClinTox CSV'
		)
		parser.add_argument('--tox21-output', required=True, help='Output Tox21 CSV')
		parser.add_argument(
			'--toxcast-output', required=True, help='Output ToxCast CSV'
		)
		parser.add_argument('--sider-output', required=True, help='Output SIDER CSV')
		args = parser.parse_args()

		main(
			coldata_path=args.coldata,
			smiles_to_cid_path=args.smiles_to_cid,
			clintox_input=args.clintox_input,
			tox21_input=args.tox21_input,
			toxcast_input=args.toxcast_input,
			sider_input=args.sider_input,
			clintox_output=args.clintox_output,
			tox21_output=args.tox21_output,
			toxcast_output=args.toxcast_output,
			sider_output=args.sider_output,
		)
