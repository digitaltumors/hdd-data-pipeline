import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import tqdm

script_dir = Path(getattr(globals().get('snakemake'), 'scriptdir', Path(__file__).resolve().parent))
if str(script_dir) not in sys.path:
	sys.path.insert(0, str(script_dir))

import utils


def iter_records(path: str):
	with open(path, "r", encoding="utf-8") as handle:
		for line in handle:
			line = line.strip()
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
	lincs_compounds = pd.read_csv(lincs_file)
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
			utils.process_single_drug(
				drug_info,
				drug_details=drug_details,
				colData=colData,
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
		print(f"Warning: {len(error_cids)} compounds failed during processing")

	colData = pd.DataFrame(colData)
	colData.to_csv(coldata_path, index=False)

	seen_bioassays = sorted(list(set(seen_bioassays)))
	seen_bioassays = [
		aid for aid in seen_bioassays if int(aid) in utils.GOLD_STANDARD_AIDS
	]

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
			outcome = "Active" if assay["activity_outcome_method"] == 2 else "Inactive"
			cpd_results[assay_idx] = outcome

		bioassay_res[cpd] = cpd_results

	bioassay_res = pd.DataFrame(
		bioassay_res, index=[f"AID_{aid}" for aid in seen_bioassays]
	)
	bioassay_res.reset_index(drop=False, inplace=True, names="Assay")
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
		parser.add_argument("-l", required=True, help="LINCS compounds CSV")
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
