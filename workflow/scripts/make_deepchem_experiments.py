import argparse
from pathlib import Path

import pandas as pd


def process_deepchem_data(
	colData: pd.DataFrame,
	data: pd.DataFrame,
	index_name: str,
	convert_to_int: bool = False,
) -> pd.DataFrame:
	measurement_cols = [
		col_name for col_name in data.columns if col_name not in ['smiles', 'mol_id']
	]
	if convert_to_int:
		data[measurement_cols] = data[measurement_cols].apply(
			lambda col: pd.to_numeric(col, errors="coerce")
		).astype("Int64")

	data = pd.merge(colData, data, left_on='SMILES', right_on='smiles')

	data = data[['Pubchem CID'] + measurement_cols].transpose()
	data = data.rename(columns=data.iloc[0])
	data = data.iloc[1:,]
	data = data.reset_index(names=index_name, drop=False)
	return data


def main(
	coldata_path: str,
	clintox_input: str,
	tox21_input: str,
	toxcast_input: str,
	sider_input: str,
	clintox_output: str,
	tox21_output: str,
	toxcast_output: str,
	sider_output: str,
) -> None:
	colData = pd.read_csv(coldata_path, usecols=['SMILES', 'Pubchem CID'])
	clintox = pd.read_csv(clintox_input)
	tox21 = pd.read_csv(tox21_input)
	toxcast = pd.read_csv(toxcast_input)
	sider = pd.read_csv(sider_input)
	Path(clintox_output).parent.mkdir(parents=True, exist_ok=True)

	clintox = process_deepchem_data(colData, clintox, 'Clinical Tox Result')
	sider = process_deepchem_data(colData, sider, 'Side Effect')
	tox21 = process_deepchem_data(colData, tox21, 'Tox Assay', True)
	toxcast = process_deepchem_data(colData, toxcast, 'Tox Assay')

	clintox.to_csv(clintox_output, index=False)
	sider.to_csv(sider_output, index=False)
	tox21.to_csv(tox21_output, index=False)
	toxcast.to_csv(toxcast_output, index=False)


def main_from_snakemake() -> None:
	main(
		coldata_path=str(snakemake.input.colData),
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
		parser.add_argument('--clintox-input', required=True, help='Input ClinTox CSV')
		parser.add_argument('--tox21-input', required=True, help='Input Tox21 CSV')
		parser.add_argument('--toxcast-input', required=True, help='Input ToxCast CSV')
		parser.add_argument('--sider-input', required=True, help='Input SIDER CSV')
		parser.add_argument('--clintox-output', required=True, help='Output ClinTox CSV')
		parser.add_argument('--tox21-output', required=True, help='Output Tox21 CSV')
		parser.add_argument('--toxcast-output', required=True, help='Output ToxCast CSV')
		parser.add_argument('--sider-output', required=True, help='Output SIDER CSV')
		args = parser.parse_args()

		main(
			coldata_path=args.coldata,
			clintox_input=args.clintox_input,
			tox21_input=args.tox21_input,
			toxcast_input=args.toxcast_input,
			sider_input=args.sider_input,
			clintox_output=args.clintox_output,
			tox21_output=args.tox21_output,
			toxcast_output=args.toxcast_output,
			sider_output=args.sider_output,
		)
