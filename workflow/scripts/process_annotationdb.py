import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Union

import pandas as pd
import tqdm

ACTIVE_OUTCOME_METHOD = 2
LOGICAL_COLUMNS = ['FDA Approved', 'In LINCS', 'In JUMP-CP', 'In OASIS', 'In GEOM']
MISSING_MEMBERSHIP_ID = '-'


def normalize_membership_value(value: object) -> str | None:
	if pd.isna(value):
		return None
	normalized = str(value).strip()
	if not normalized or normalized.upper() == 'NA':
		return None
	return normalized


def format_membership_ids(values: pd.Series) -> str:
	unique_values = sorted(
		{
			normalized
			for value in values
			if (normalized := normalize_membership_value(value)) is not None
		}
	)
	if not unique_values:
		return MISSING_MEMBERSHIP_ID
	return '|'.join(unique_values)


def build_membership_map(
	frame: pd.DataFrame,
	inchikey_col: str,
	id_col: str,
	filter_col: str | None = None,
	filter_value: str | None = None,
) -> dict[str, str]:
	required_columns = [inchikey_col, id_col]
	if filter_col is not None:
		required_columns.append(filter_col)
	missing_columns = [
		column for column in required_columns if column not in frame.columns
	]
	if missing_columns:
		raise ValueError(
			'Membership table is missing required columns: '
			+ ', '.join(missing_columns)
		)

	membership_frame = frame.copy()
	if filter_col is not None and filter_value is not None:
		expected = filter_value.strip().casefold()
		membership_frame = membership_frame[
			membership_frame[filter_col].map(
				lambda value: (
					(normalize_membership_value(value) or '').casefold() == expected
				)
			)
		]

	membership_frame['_normalized_inchikey'] = membership_frame[inchikey_col].map(
		normalize_membership_value
	)
	membership_frame = membership_frame[
		membership_frame['_normalized_inchikey'].notna()
	]
	if membership_frame.empty:
		return {}

	return (
		membership_frame.groupby('_normalized_inchikey', sort=True)[id_col]
		.agg(format_membership_ids)
		.to_dict()
	)


def append_membership_columns(
	col_data: defaultdict(list),
	membership_map: dict[str, str],
	inchikey: object,
	flag_col: str,
	id_col: str,
) -> None:
	normalized_inchikey = normalize_membership_value(inchikey)
	membership_id = (
		membership_map.get(normalized_inchikey)
		if normalized_inchikey is not None
		else None
	)
	if membership_id is None:
		col_data[flag_col].append(False)
		col_data[id_col].append(MISSING_MEMBERSHIP_ID)
		return

	col_data[flag_col].append(True)
	col_data[id_col].append(membership_id)


def process_single_drug(
	drug_info: Dict[str, Union[int, float, str]],
	drug_details: Dict[str, Union[int, float, str]],
	col_data: defaultdict(list),
	all_bioassays: Dict[str, Dict],
	seen_bioassays: List[int],
	lincs_membership: dict[str, str],
	jump_cp_membership: dict[str, str],
	oasis_membership: dict[str, str],
	geom_membership: dict[str, str],
	blood_brain_perm: pd.DataFrame,
	cids: List,
) -> None:
	# Molecule Name
	drug_name = (
		drug_info.get('name')
		or drug_info.get('mapped_name')
		or drug_details.get('title')
	)
	cid = drug_info.get('cid') or drug_details.get('cid')
	inchikey = drug_info.get('inchikey') or drug_details.get('inchikey')
	smiles_str = drug_info.get('smiles') or drug_details.get('smiles')

	col_data['Molecule Name'].append(drug_name)
	col_data['Pubchem CID'].append(cid)
	col_data['InChIKey'].append(inchikey)
	col_data['SMILES'].append(smiles_str)
	cids.append(cid)
	# store for later use

	col_data['Molecular Formula'].append(drug_details['molecular_formula'])
	col_data['IUPAC Name'].append(drug_details['iupac_name'])
	col_data['ChEMBL ID'].append(drug_details['molecule_chembl_id'])
	col_data['PubChem 2D Fingerprint'].append(drug_details['fingerprint_2d'])

	# MOA and Approval
	mechanisms = drug_details['mechanisms']
	if len(mechanisms) == 0:
		col_data['Mechanism of Action'].append('None')
	else:
		col_data['Mechanism of Action'].append(
			drug_details['mechanisms'][0]['mechanism_of_action']
		)

	col_data['FDA Approved'].append(drug_details['fda_approval'])
	## Add in Molecular Information (molecular weight + Lipinski Filters)
	## 	for the curious: https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five
	col_data['Molecular Weight'].append(drug_details['molecular_weight'])
	col_data['XlogP'].append(drug_details['xlogp'])
	col_data['Hydrogen Bond Donors'].append(drug_details['h_bond_donor_count'])
	col_data['Hydrogen Bond Acceptors'].append(drug_details['h_bond_acceptor_count'])
	col_data['Exact Molecular Mass'].append(drug_details['exact_mass'])
	col_data['DILI Severity'].append(drug_details['toxicity']['dili_severity_grade'])
	col_data['DILI Annotation'].append(drug_details['toxicity']['dili_annotation'])
	col_data['Hepatotoxicity Likelihood (Detailed)'].append(
		drug_details['toxicity']['hepatotoxicity_likelihood_score']
	)
	hls = drug_details['toxicity']['hepatotoxicity_likelihood_score']
	score = pd.NA
	if isinstance(hls, str) and hls:
		parts = hls.split(':', 1)
		if len(parts) > 1:
			score = parts[1].lstrip().split()[0] if parts[1].strip() else pd.NA
	col_data['Hepatotoxiciy Likelihood (Score)'].append(score)

	append_membership_columns(
		col_data,
		lincs_membership,
		inchikey,
		flag_col='In LINCS',
		id_col='LINCS ID',
	)
	append_membership_columns(
		col_data,
		jump_cp_membership,
		inchikey,
		flag_col='In JUMP-CP',
		id_col='JUMP-CP ID',
	)
	append_membership_columns(
		col_data,
		oasis_membership,
		inchikey,
		flag_col='In OASIS',
		id_col='OASIS ID',
	)
	append_membership_columns(
		col_data,
		geom_membership,
		inchikey,
		flag_col='In GEOM',
		id_col='GEOM Source SMILES',
	)

	blood_brain = blood_brain_perm[blood_brain_perm['smiles'] == smiles_str]
	if blood_brain.shape[0] == 0:
		col_data['BBB Permeable'].append('Unknown')
	else:
		col_data['BBB Permeable'].append(blood_brain['p_np'].to_numpy()[0])

	# Cache bioassays for post-processing
	all_bioassays[drug_details['cid']] = drug_details['bioassays']
	seen_bioassays.extend([assay['aid'] for assay in drug_details['bioassays']])


def iter_records(path: str) -> Iterator[dict]:
	with Path(path).open('r', encoding='utf-8') as handle:
		for raw_line in handle:
			line = raw_line.strip()
			if not line:
				continue
			yield json.loads(line)


def format_logical_values(col_data: pd.DataFrame) -> pd.DataFrame:
	for column in LOGICAL_COLUMNS:
		if column not in col_data.columns:
			continue
		col_data[column] = col_data[column].replace({True: 'TRUE', False: 'FALSE'})
	return col_data


def load_membership_maps(
	lincs_file: str,
	jump_cp_file: str,
	oasis_file: str,
	geom_file: str,
) -> tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]:
	lincs_compounds = pd.read_csv(lincs_file, sep='\t')
	jump_cp_compounds = pd.read_csv(jump_cp_file)
	oasis_compounds = pd.read_csv(oasis_file, sep='\t')
	geom_compounds = pd.read_csv(geom_file, sep='\t')

	lincs_membership = build_membership_map(
		lincs_compounds,
		inchikey_col='inchi_key',
		id_col='pert_id',
	)
	jump_cp_membership = build_membership_map(
		jump_cp_compounds,
		inchikey_col='Metadata_InChIKey',
		id_col='Metadata_JCP2022',
	)
	oasis_membership = build_membership_map(
		oasis_compounds,
		inchikey_col='InChIKey',
		id_col='OASIS.ID',
		filter_col='Perturbation.Type',
		filter_value='treatment',
	)
	geom_membership = build_membership_map(
		geom_compounds,
		inchikey_col='InChIKey',
		id_col='GEOM.Source.SMILES',
		filter_col='GEOM.Source.Subset',
		filter_value='drugs',
	)
	return lincs_membership, jump_cp_membership, oasis_membership, geom_membership


def write_bioassay_matrix(
	all_bioassays: Dict[str, Dict],
	cids: list,
	seen_bioassays: list,
	bioassays_path: Path,
) -> None:
	seen_bioassays = sorted(list(set(seen_bioassays)))
	aid_to_idx = {seen_bioassays[i]: i for i in range(len(seen_bioassays))}
	num_assays = len(seen_bioassays)
	bioassay_res = defaultdict(list)

	for cpd in cids:
		assay_subset = all_bioassays[cpd]
		cpd_results = num_assays * ['Not Measured']

		for assay in assay_subset:
			assay_id = assay['aid']
			if assay_id not in aid_to_idx:
				continue
			assay_idx = aid_to_idx[assay_id]
			outcome = (
				'Active'
				if assay['activity_outcome_method'] == ACTIVE_OUTCOME_METHOD
				else 'Inactive'
			)
			cpd_results[assay_idx] = outcome

		bioassay_res[cpd] = cpd_results

	bioassay_res = pd.DataFrame(
		bioassay_res, index=[f'AID_{aid}' for aid in seen_bioassays]
	)
	bioassay_res = bioassay_res.reset_index(drop=False, names='Assay')
	bioassay_res.to_csv(bioassays_path, index=False)


def main(
	input_path: str,
	lincs_file: str,
	jump_cp_file: str,
	oasis_file: str,
	geom_file: str,
	bbbp_file: str,
	coldata_output: str,
	bioassays_output: str,
) -> None:
	col_data, all_bioassays = defaultdict(list), defaultdict(list)
	seen_bioassays, cids = [], []
	(
		lincs_membership,
		jump_cp_membership,
		oasis_membership,
		geom_membership,
	) = load_membership_maps(
		lincs_file=lincs_file,
		jump_cp_file=jump_cp_file,
		oasis_file=oasis_file,
		geom_file=geom_file,
	)
	blood_brain_perm = pd.read_csv(bbbp_file)
	coldata_path = Path(coldata_output)
	bioassays_path = Path(bioassays_output)
	coldata_path.parent.mkdir(parents=True, exist_ok=True)
	bioassays_path.parent.mkdir(parents=True, exist_ok=True)

	error_cids = []

	for record in tqdm.tqdm(iter_records(input_path)):
		drug_info = record.get('drug_info')
		drug_details = record.get('drug_details')
		if drug_info is None or drug_details is None:
			continue

		keys_before = set(col_data.keys())
		coldata_lengths = {k: len(v) for k, v in col_data.items()}
		seen_len = len(seen_bioassays)
		cids_len = len(cids)

		try:
			process_single_drug(
				drug_info,
				drug_details=drug_details,
				col_data=col_data,
				all_bioassays=all_bioassays,
				seen_bioassays=seen_bioassays,
				lincs_membership=lincs_membership,
				jump_cp_membership=jump_cp_membership,
				oasis_membership=oasis_membership,
				geom_membership=geom_membership,
				blood_brain_perm=blood_brain_perm,
				cids=cids,
			)
		except Exception:
			cid = drug_info.get('cid') if isinstance(drug_info, dict) else None
			if cid is not None:
				error_cids.append(cid)
			for k in list(col_data.keys()):
				if k not in keys_before:
					del col_data[k]
				else:
					col_data[k] = col_data[k][: coldata_lengths.get(k, 0)]
			seen_bioassays[:] = seen_bioassays[:seen_len]
			cids[:] = cids[:cids_len]
			if cid is not None:
				all_bioassays.pop(cid, None)

	if error_cids:
		tqdm.tqdm.write(
			f'Warning: {len(error_cids)} compounds failed during processing'
		)

	col_data = pd.DataFrame(col_data)
	col_data = format_logical_values(col_data)
	col_data.to_csv(coldata_path, index=False)
	write_bioassay_matrix(all_bioassays, cids, seen_bioassays, bioassays_path)


def main_from_snakemake() -> None:
	main(
		input_path=str(snakemake.input.raw_data),
		lincs_file=str(snakemake.input.lincs_file),
		jump_cp_file=str(snakemake.input.jump_file),
		oasis_file=str(snakemake.input.oasis_file),
		geom_file=str(snakemake.input.geom_file),
		bbbp_file=str(snakemake.input.bbbp_file),
		coldata_output=str(snakemake.output.colData),
		bioassays_output=str(snakemake.output.bioassays),
	)


if __name__ == '__main__':
	if 'snakemake' in globals():
		main_from_snakemake()
	else:
		parser = argparse.ArgumentParser(
			prog='process_annotationdb',
			description='Generate colData and bioassays from AnnotationDB JSONL',
		)
		parser.add_argument(
			'-i', required=True, help='Input JSONL from fetch_annotationdb'
		)
		parser.add_argument('-l', required=True, help='LINCS compounds TSV')
		parser.add_argument('-j', required=True, help='JUMP-CP compounds CSV')
		parser.add_argument('-o', required=True, help='OASIS HDD membership TSV')
		parser.add_argument('-g', required=True, help='GEOM HDD membership TSV')
		parser.add_argument('-b', required=True, help='Blood brain barrier CSV')
		parser.add_argument('-c', required=True, help='Output colData CSV')
		parser.add_argument('-a', required=True, help='Output bioassays CSV')
		args = parser.parse_args()

		main(
			input_path=args.i,
			lincs_file=args.l,
			jump_cp_file=args.j,
			oasis_file=args.o,
			geom_file=args.g,
			bbbp_file=args.b,
			coldata_output=args.c,
			bioassays_output=args.a,
		)
