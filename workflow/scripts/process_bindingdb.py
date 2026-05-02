import argparse
from pathlib import Path

import pandas as pd

CHUNKSIZE = 250_000


def normalize_cid(values: pd.Series) -> pd.Series:
	return pd.to_numeric(values, errors='coerce').astype('Int64')


def process_bindingdb(
	raw_zip: str,
	output_path: str,
	processing: dict,
) -> None:
	keep_cols = list(processing['keep_cols'])
	output_cols = list(processing['output_cols'])
	organism_col = processing['organism_col']
	organisms = set(processing['keep_organisms'])
	assay_col = processing['assay_col']
	cid_col = processing['cid_col']

	output = Path(output_path)
	output.parent.mkdir(parents=True, exist_ok=True)

	rows_written = 0
	chunks_read = 0
	header = True
	for chunk in pd.read_csv(
		raw_zip,
		sep='\t',
		compression='zip',
		usecols=keep_cols,
		dtype='string',
		chunksize=CHUNKSIZE,
		low_memory=False,
	):
		chunks_read += 1
		chunk[cid_col] = normalize_cid(chunk[cid_col])
		filtered = chunk[
			chunk[organism_col].isin(organisms)
			& chunk[assay_col].isna()
			& chunk[cid_col].notna()
		].copy()
		if filtered.empty:
			continue

		filtered = filtered[output_cols]
		filtered.to_csv(output, mode='w' if header else 'a', index=False, header=header)
		rows_written += len(filtered)
		header = False

	if rows_written == 0:
		pd.DataFrame(columns=output_cols).to_csv(output, index=False)

	print(  # noqa: T201
		f'[process_bindingdb] chunks={chunks_read} rows={rows_written} output={output}',
		flush=True,
	)


def main_from_snakemake() -> None:
	process_bindingdb(
		raw_zip=str(snakemake.input.raw_zip),
		output_path=str(snakemake.output.cleaned_data),
		processing=dict(snakemake.params.processing),
	)


if __name__ == '__main__':
	if 'snakemake' in globals():
		main_from_snakemake()
	else:
		parser = argparse.ArgumentParser(
			prog='process_bindingdb',
			description='Filter BindingDB records for the HDD BindingDB assay',
		)
		parser.add_argument('--raw-zip', required=True, help='BindingDB TSV zip path')
		parser.add_argument('--output', required=True, help='Cleaned CSV output path')
		args = parser.parse_args()

		default_processing = {
			'keep_organisms': ['Homo sapiens (Human)', 'Homo sapiens'],
			'organism_col': 'Target Source Organism According to Curator or DataSource',
			'assay_col': 'PubChem AID',
			'cid_col': 'PubChem CID',
			'keep_cols': [
				'PubChem CID',
				'PubChem AID',
				'Target Name',
				'Target Source Organism According to Curator or DataSource',
				'Ki (nM)',
				'Kd (nM)',
				'Curation/DataSource',
			],
			'output_cols': ['PubChem CID', 'Target Name', 'Ki (nM)', 'Kd (nM)'],
		}
		process_bindingdb(
			raw_zip=args.raw_zip,
			output_path=args.output,
			processing=default_processing,
		)
