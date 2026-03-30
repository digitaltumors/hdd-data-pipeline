from collections import defaultdict
from itertools import product
from pathlib import Path
from typing import Dict, List

import pandas as pd
from damply import dirs
from rdkit import Chem
from rdkit.Chem import AllChem


def make_fingerprint_generators(
	radius_list: List[int], dimension_list: List[int]
) -> Dict[str, rdkit.Chem.rdFingerprintGenerator]:
	fingerprint_generators = {}
	for radius, dimension in product(radius_list, dimension_list):
		fpgen = AllChem.GetMorganGenerator(radius=radius, fpSize=dimension)
		fingerprint_generators[f'Morgan({radius},{dimension})'] = fpgen

	return fingerprint_generators


def main(
	coldata_path: str,
	output_paths: List[str],
	radius_list=[2, 3],
	dimension_list=[512, 1024, 2048],
):
	fingerprint_generators = make_fingerprint_generators(radius_list, dimension_list)
	fp_data = {k: defaultdict(list) for k in fingerprint_generators}
	if output_paths:
		output_by_name = {Path(path).name: Path(path) for path in output_paths}
		Path(output_paths[0]).parent.mkdir(parents=True, exist_ok=True)
	else:
		outpath = dirs.PROCDATA / 'experiments' / 'fingerprints'
		outpath.mkdir(parents=True, exist_ok=True)
		output_by_name = {}
		for fp_type in fingerprint_generators:
			fp_str = '.'.join(fp_type.replace('(', '.').replace(')', '').split(','))
			output_by_name[f'{fp_str}.csv'] = outpath / f'{fp_str}.csv'
	colData = pd.read_csv(coldata_path, usecols=['Pubchem CID', 'SMILES'])

	for _, row in colData.iterrows():
		cid = row['Pubchem CID']
		smiles_str = row['SMILES']

		rdk_mol = Chem.MolFromSmiles(smiles_str)

		for fp_type in fingerprint_generators:
			fp_data[fp_type]['CID'].append(cid)
			generator = fingerprint_generators[fp_type]
			fp = generator.GetCountFingerprintAsNumPy(rdk_mol)
			for j in range(fp.shape[0]):
				vname = f'V{j + 1}'
				fp_data[fp_type][vname].append(fp[j])

	for fp_type in fp_data:
		fp_str = '.'.join(fp_type.replace('(', '.').replace(')', '').split(','))

		fp_matrix = pd.DataFrame(fp_data[fp_type])
		fp_matrix = fp_matrix.transpose()

		fp_matrix = fp_matrix.rename(columns=fp_matrix.iloc[0])

		fp_matrix = fp_matrix.iloc[1:,]
		fp_matrix = fp_matrix.reset_index(drop=True)
		output_path = output_by_name[f'{fp_str}.csv']
		fp_matrix.to_csv(output_path, index=False)


def main_from_snakemake() -> None:
	main(
		coldata_path=str(snakemake.input[0]),
		output_paths=[str(path) for path in snakemake.output.fingerprints],
		radius_list=list(snakemake.params.radius_list),
		dimension_list=list(snakemake.params.dim_list),
	)


if __name__ == '__main__':
	if 'snakemake' in globals():
		main_from_snakemake()
	else:
		main(
			coldata_path=str(dirs.PROCDATA / 'colData.csv'),
			output_paths=[],
		)
