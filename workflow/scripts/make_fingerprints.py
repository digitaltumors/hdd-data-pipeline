import os
import shutil
import tempfile
from contextlib import ExitStack
from itertools import product
from pathlib import Path
from typing import Dict, List, Optional, Sequence, TextIO

import pandas as pd
from damply import dirs
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog('rdApp.*')


def fingerprint_stem(radius: int, dimension: int) -> str:
	return f'Morgan.{radius}.{dimension}'


def make_fingerprint_generators(
	radius_list: List[int], dimension_list: List[int]
) -> Dict[str, object]:
	fingerprint_generators = {}
	for radius, dimension in product(radius_list, dimension_list):
		fpgen = AllChem.GetMorganGenerator(radius=radius, fpSize=dimension)
		fingerprint_generators[fingerprint_stem(radius, dimension)] = fpgen

	return fingerprint_generators


def build_output_by_name(
	output_paths: List[str],
	radius_list: List[int],
	dimension_list: List[int],
) -> Dict[str, Path]:
	if output_paths:
		output_by_name = {Path(path).stem: Path(path) for path in output_paths}
		Path(output_paths[0]).parent.mkdir(parents=True, exist_ok=True)
		return output_by_name

	outpath = dirs.PROCDATA / 'experiments' / 'fingerprints'
	outpath.mkdir(parents=True, exist_ok=True)
	return {
		fingerprint_stem(radius, dimension): outpath
		/ f'{fingerprint_stem(radius, dimension)}.mtx'
		for radius, dimension in product(radius_list, dimension_list)
	}


def write_matrix_market(
	output_path: Path,
	body_path: Path,
	num_features: int,
	num_compounds: int,
	num_nonzero: int,
) -> None:
	output_path.parent.mkdir(parents=True, exist_ok=True)
	with output_path.open('w', encoding='utf-8') as handle:
		handle.write('%%MatrixMarket matrix coordinate integer general\n')
		handle.write('%\n')
		handle.write(f'{num_features} {num_compounds} {num_nonzero}\n')
		with body_path.open('r', encoding='utf-8') as body_handle:
			shutil.copyfileobj(body_handle, handle)


def open_body_handles(
	output_by_name: Dict[str, Path],
) -> tuple[ExitStack, Dict[str, TextIO], Dict[str, Path]]:
	stack = ExitStack()
	body_paths = {}
	body_handles = {}
	for name, output_path in output_by_name.items():
		fd, temp_name = tempfile.mkstemp(
			dir=output_path.parent,
			prefix=f'{name}.',
			suffix='.body',
		)
		os.close(fd)
		body_path = Path(temp_name)
		temp_handle = stack.enter_context(body_path.open('w', encoding='utf-8'))
		body_paths[name] = body_path
		body_handles[name] = temp_handle
	return stack, body_handles, body_paths


def stream_fingerprint_bodies(
	coldata: pd.DataFrame,
	fingerprint_generators: Dict[str, object],
	body_handles: Dict[str, TextIO],
) -> tuple[Dict[str, int], int, list[dict[str, object]]]:
	nonzero_counts = {name: 0 for name in fingerprint_generators}
	invalid_smiles = 0
	fingerprint_columns = []

	for _, row in coldata.iterrows():
		smiles_str = row['SMILES']
		if pd.isna(smiles_str):
			invalid_smiles += 1
			continue
		molecule = Chem.MolFromSmiles(smiles_str)
		if molecule is None:
			invalid_smiles += 1
			continue

		compound_index = len(fingerprint_columns) + 1
		fingerprint_columns.append(
			{
				'Fingerprint.Column': compound_index,
				'HDD.Compound.ID': row['HDD.Compound.ID'],
				'SMILES': smiles_str,
			}
		)

		for name, generator in fingerprint_generators.items():
			fingerprint = generator.GetCountFingerprint(molecule)
			nonzero_elements = fingerprint.GetNonzeroElements()
			if not nonzero_elements:
				continue

			handle = body_handles[name]
			for feature_index, count in nonzero_elements.items():
				handle.write(f'{feature_index + 1} {compound_index} {int(count)}\n')
			nonzero_counts[name] += len(nonzero_elements)

	return nonzero_counts, invalid_smiles, fingerprint_columns


def finalize_outputs(
	output_by_name: Dict[str, Path],
	body_paths: Dict[str, Path],
	fingerprint_dims: Dict[str, int],
	num_compounds: int,
	nonzero_counts: Dict[str, int],
) -> None:
	for name, output_path in output_by_name.items():
		write_matrix_market(
			output_path=output_path,
			body_path=body_paths[name],
			num_features=fingerprint_dims[name],
			num_compounds=num_compounds,
			num_nonzero=nonzero_counts[name],
		)


def main(
	coldata_path: str,
	output_paths: List[str],
	fingerprint_columns_path: str,
	radius_list: Optional[Sequence[int]] = None,
	dimension_list: Optional[Sequence[int]] = None,
) -> None:
	radius_list = list(radius_list or [2, 3])
	dimension_list = list(dimension_list or [512, 1024, 2048])
	output_by_name = build_output_by_name(output_paths, radius_list, dimension_list)
	fingerprint_generators = make_fingerprint_generators(radius_list, dimension_list)
	fingerprint_dims = {
		fingerprint_stem(radius, dimension): dimension
		for radius, dimension in product(radius_list, dimension_list)
	}
	coldata = pd.read_csv(coldata_path, usecols=['HDD.Compound.ID', 'SMILES'])
	stack, body_handles, body_paths = open_body_handles(output_by_name)
	invalid_smiles = 0
	fingerprint_column_count = 0
	try:
		with stack:
			nonzero_counts, invalid_smiles, fingerprint_columns = (
				stream_fingerprint_bodies(
					coldata,
					fingerprint_generators,
					body_handles,
				)
			)
		fingerprint_column_count = len(fingerprint_columns)
		Path(fingerprint_columns_path).parent.mkdir(parents=True, exist_ok=True)
		pd.DataFrame(
			fingerprint_columns,
			columns=['Fingerprint.Column', 'HDD.Compound.ID', 'SMILES'],
		).to_csv(fingerprint_columns_path, sep='\t', index=False)
		finalize_outputs(
			output_by_name=output_by_name,
			body_paths=body_paths,
			fingerprint_dims=fingerprint_dims,
			num_compounds=fingerprint_column_count,
			nonzero_counts=nonzero_counts,
		)
	finally:
		for body_path in body_paths.values():
			if body_path.exists():
				body_path.unlink()

	print(  # noqa: T201
		'[make_fingerprints] '
		f'compounds={len(coldata)} '
		f'fingerprinted={fingerprint_column_count} '
		f'invalid_smiles={invalid_smiles}',
		flush=True,
	)


def main_from_snakemake() -> None:
	main(
		coldata_path=str(snakemake.input[0]),
		output_paths=[str(path) for path in snakemake.output.fingerprints],
		fingerprint_columns_path=str(snakemake.output.fingerprint_columns),
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
			fingerprint_columns_path=str(
				dirs.PROCDATA
				/ 'experiments'
				/ 'fingerprints'
				/ 'fingerprint_columns.tsv'
			),
		)
