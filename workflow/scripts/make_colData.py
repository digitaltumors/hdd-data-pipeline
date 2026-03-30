import argparse
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests
import tqdm
from damply import dirs
from process_annotationdb import (
	ACTIVE_OUTCOME_METHOD,
	GOLD_STANDARD_AIDS,
	process_single_drug,
)


def fetch_json(url: str, retries: int = 3, timeout: int = 30):
	for attempt in range(retries):
		try:
			resp = requests.get(url, timeout=timeout)
			resp.raise_for_status()
			return resp.json()
		except Exception:
			if attempt == retries - 1:
				raise
			time.sleep(2 * (attempt + 1))


def main(db_url: str, lincs_file: str, jump_cp_file: str, bbbp_file: str) -> None:
	# set up storage for results
	# default dict allows colData to be easily converted into a dataframe

	colData, all_bioassays = defaultdict(list), defaultdict(list)
	seen_bioassays, cids = [], []

	error_cids = []

	# download all the compounds in annotationdb
	response = fetch_json(db_url)

	# read in jumpcp and lincs
	jump_cp_compounds = pd.read_csv(jump_cp_file)
	lincs_compounds = pd.read_csv(lincs_file)
	blood_brain_perm = pd.read_csv(bbbp_file)
	# print(blood_brain_perm)

	for drug_info in tqdm.tqdm(response):
		keys_before = set(colData.keys())
		coldata_lengths = {k: len(v) for k, v in colData.items()}
		seen_len = len(seen_bioassays)
		cids_len = len(cids)

		try:
			drug_query_url = f'https://annotationdb.bhklab.ca/compound/many?compounds={drug_info["cid"]}&format=json&bioassay=true&mechanism=true&toxicity=true'
			drug_details = fetch_json(drug_query_url)[0]

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
			error_cids.append(drug_info['cid'])
			for k in list(colData.keys()):
				if k not in keys_before:
					del colData[k]
				else:
					colData[k] = colData[k][: coldata_lengths.get(k, 0)]
			seen_bioassays[:] = seen_bioassays[:seen_len]
			cids[:] = cids[:cids_len]
			all_bioassays.pop(drug_info['cid'], None)

	# write and store colData
	colData = pd.DataFrame(colData)
	colData.to_csv(dirs.PROCDATA / 'colData.csv', index=False)

	# process the bioassays since we have the data here.
	seen_bioassays = sorted(list(set(seen_bioassays)))
	seen_bioassays = [aid for aid in seen_bioassays if int(aid) in GOLD_STANDARD_AIDS]

	aid_to_idx = {seen_bioassays[i]: i for i in range(len(seen_bioassays))}
	num_assays = len(seen_bioassays)
	num_cpds = len(cids)
	cid_to_idx = {cids[i]: i for i in range(len(cids))}
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

	outpath = Path(dirs.PROCDATA / 'experiments')
	bioassay_res = pd.DataFrame(
		bioassay_res, index=[f'AID_{aid}' for aid in seen_bioassays]
	)
	bioassay_res.reset_index(drop=False, inplace=True, names='Assay')
	# bioassay_res = bioassay_res[blood_brain_perm['Assay'].isin([f"AID_{aid}" for aid in GOLD_STANDARD_AIDS])]
	bioassay_res.to_csv(outpath / 'bioassays.csv', index=False)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		prog='make_ColData ',
		description='Generate the colData and parse other annotationdb data',
	)
	parser.add_argument('-u', help='route to all compounds in database')
	parser.add_argument('-l', help='lincs')
	parser.add_argument('-j', help='jump cp')
	parser.add_argument('-b', help='blood brain barrier')

	args = parser.parse_args()

	main(db_url=args.u, lincs_file=args.l, jump_cp_file=args.j, bbbp_file=args.b)
