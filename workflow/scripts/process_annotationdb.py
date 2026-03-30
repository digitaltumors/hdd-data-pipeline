import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Union

import pandas as pd
import tqdm

GOLD_STANDARD_AIDS = [
	485290,
	1508612,
	1645840,
	1645841,
	1645842,
	492947,
	1030,
	743075,
	743080,
	588795,
	2101,
	602179,
	504327,
	995,
	493208,
	1777,
	1631,
	743094,
	651631,
	504847,
	1159551,
	1259242,
	1259241,
	743012,
	1224868,
	1224880,
	1346977,
	743014,
	1224870,
	1224872,
	1224874,
	1224877,
	1224885,
	1224886,
	743015,
	1224873,
	1224887,
	1224889,
	1224867,
	720516,
	651632,
	651634,
]
ACTIVE_OUTCOME_METHOD = 2


def process_single_drug(
	drug_info: Dict[str, Union[int, float, str]],
	drug_details: Dict[str, Union[int, float, str]],
	col_data: defaultdict(list),
	all_bioassays: Dict[str, Dict],
	seen_bioassays: List[int],
	lincs_compounds: pd.DataFrame,
	jump_cp_compounds: pd.DataFrame,
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
	col_data['Hydrogen Bond Acceptors'].append(
		drug_details['h_bond_acceptor_count']
	)
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

	## Check Against The Broad Data
	# print("pingo herebo")
	# print(lincs_compounds)
	# print(lincs_compounds['inchi_key'].value_counts())
	l1k_subset = lincs_compounds[lincs_compounds['inchi_key'] == inchikey]
	jump_subset = jump_cp_compounds[
		jump_cp_compounds['Metadata_InChIKey'] == inchikey
	]

	if l1k_subset.shape[0] == 0:
		# print("thrig  plibbus")
		col_data['In L1000'].append(False)
		col_data['L1000 ID'].append('-')
	else:
		col_data['In L1000'].append(True)
		col_data['L1000 ID'].append(l1k_subset['pert_id'].to_numpy()[0])

	if jump_subset.shape[0] == 0:
		col_data['In JUMP-CP'].append(False)
		col_data['JUMP-CP ID'].append('-')
	else:
		col_data['In JUMP-CP'].append(True)
		col_data['JUMP-CP ID'].append(jump_subset['Metadata_JCP2022'].to_numpy()[0])

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


def main(
	input_path: str,
	lincs_file: str,
	jump_cp_file: str,
	bbbp_file: str,
	coldata_output: str,
	bioassays_output: str,
) -> None:
	colData, all_bioassays = defaultdict(list), defaultdict(list)
	seen_bioassays, cids = [], []
	lincs_compounds = pd.read_csv(lincs_file, sep='\t')
	jump_cp_compounds = pd.read_csv(jump_cp_file)
	blood_brain_perm = pd.read_csv(bbbp_file)
	coldata_path = Path(coldata_output)
	bioassays_path = Path(bioassays_output)
	coldata_path.parent.mkdir(parents=True, exist_ok=True)
	bioassays_path.parent.mkdir(parents=True, exist_ok=True)

	error_cids = []

	for record in tqdm.tqdm(iter_records(input_path)):
		drug_info = record.get("drug_info")
		drug_details = record.get("drug_details")
		if drug_info is None or drug_details is None:
			continue

		keys_before = set(colData.keys())
		coldata_lengths = {k: len(v) for k, v in colData.items()}
		seen_len = len(seen_bioassays)
		cids_len = len(cids)

		try:
			process_single_drug(
				drug_info,
				drug_details=drug_details,
				col_data=colData,
				all_bioassays=all_bioassays,
				seen_bioassays=seen_bioassays,
				lincs_compounds=lincs_compounds,
				jump_cp_compounds=jump_cp_compounds,
				blood_brain_perm=blood_brain_perm,
				cids=cids,
			)
		except Exception:
			cid = drug_info.get("cid") if isinstance(drug_info, dict) else None
			if cid is not None:
				error_cids.append(cid)
			for k in list(colData.keys()):
				if k not in keys_before:
					del colData[k]
				else:
					colData[k] = colData[k][: coldata_lengths.get(k, 0)]
			seen_bioassays[:] = seen_bioassays[:seen_len]
			cids[:] = cids[:cids_len]
			if cid is not None:
				all_bioassays.pop(cid, None)

	if error_cids:
		tqdm.tqdm.write(
			f'Warning: {len(error_cids)} compounds failed during processing'
		)

	colData = pd.DataFrame(colData)
	colData.to_csv(coldata_path, index=False)

	seen_bioassays = sorted(list(set(seen_bioassays)))
	seen_bioassays = [aid for aid in seen_bioassays if int(aid) in GOLD_STANDARD_AIDS]

	aid_to_idx = {seen_bioassays[i]: i for i in range(len(seen_bioassays))}
	num_assays = len(seen_bioassays)
	bioassay_res = defaultdict(list)

	for cpd in cids:
		assay_subset = all_bioassays[cpd]
		cpd_results = num_assays * ["Not Measured"]

		for assay in assay_subset:
			assay_id = assay["aid"]
			if assay_id not in aid_to_idx:
				continue
			assay_idx = aid_to_idx[assay_id]
			outcome = (
				"Active"
				if assay["activity_outcome_method"] == ACTIVE_OUTCOME_METHOD
				else "Inactive"
			)
			cpd_results[assay_idx] = outcome

		bioassay_res[cpd] = cpd_results

	bioassay_res = pd.DataFrame(
		bioassay_res, index=[f"AID_{aid}" for aid in seen_bioassays]
	)
	bioassay_res = bioassay_res.reset_index(drop=False, names="Assay")
	bioassay_res.to_csv(bioassays_path, index=False)


def main_from_snakemake() -> None:
	main(
		input_path=str(snakemake.input.raw_data),
		lincs_file=str(snakemake.input.lincs_file),
		jump_cp_file=str(snakemake.input.jump_file),
		bbbp_file=str(snakemake.input.bbbp_file),
		coldata_output=str(snakemake.output.colData),
		bioassays_output=str(snakemake.output.bioassays),
	)


if __name__ == "__main__":
	if "snakemake" in globals():
		main_from_snakemake()
	else:
		parser = argparse.ArgumentParser(
			prog="process_annotationdb",
			description="Generate colData and bioassays from AnnotationDB JSONL",
		)
		parser.add_argument("-i", required=True, help="Input JSONL from fetch_annotationdb")
		parser.add_argument("-l", required=True, help="LINCS compounds TSV")
		parser.add_argument("-j", required=True, help="JUMP-CP compounds CSV")
		parser.add_argument("-b", required=True, help="Blood brain barrier CSV")
		parser.add_argument("-c", required=True, help="Output colData CSV")
		parser.add_argument("-a", required=True, help="Output bioassays CSV")
		args = parser.parse_args()

		main(
			input_path=args.i,
			lincs_file=args.l,
			jump_cp_file=args.j,
			bbbp_file=args.b,
			coldata_output=args.c,
			bioassays_output=args.a,
		)
