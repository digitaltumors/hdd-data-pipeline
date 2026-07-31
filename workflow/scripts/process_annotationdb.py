import argparse
import hashlib
import json
import math
import shutil
from pathlib import Path
from typing import Iterator

import pandas as pd
import tqdm

ACTIVE_OUTCOME_METHOD = 2
BASE62_ALPHABET = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
DEFAULT_HASH_LENGTH = 12
HASH_LENGTH_INCREMENT = 2
MAX_HASH_LENGTH = 32
TOXICITY_SOURCE_COLUMNS = {
	'LTKB': {
		'dili_severity_grade': 'LTKB.DILI.Severity.Grade',
		'dili_annotation': 'LTKB.DILI.Annotation',
	},
	'DILIrank': {
		'dili_severity_grade': 'DILIrank.DILI.Severity.Grade',
		'dili_annotation': 'DILIrank.DILI.Annotation',
	},
	'DILIst': {
		'dili_annotation': 'DILIst.DILI.Annotation',
	},
	'Livertox': {
		'hepatotoxicity_likelihood_score': ('Livertox.Hepatotoxicity.Likelihood.Score'),
		'hepatotoxicity_likelihood_score_reasoning': (
			'Livertox.Hepatotoxicity.Likelihood.Reason'
		),
	},
}
TOXICITY_COLUMNS = [
	column
	for source_columns in TOXICITY_SOURCE_COLUMNS.values()
	for column in source_columns.values()
]
SOURCE_COLUMNS = {
	'jump': 'JUMP.CP.ID',
	'oasis': 'OASIS.ID',
	'geom': 'GEOM.Source.SMILES',
	'lincs': 'LINCS.CMap.Name',
	'ctrp': 'CTRP.Master.CPD.ID',
	'nci60': 'NCI60.NSC',
}
MEMBERSHIP_COLUMNS = {
	'jump': 'In.JUMP.CP',
	'oasis': 'In.OASIS',
	'geom': 'In.GEOM',
	'lincs': 'In.LINCS',
	'ctrp': 'In.CTRP',
	'nci60': 'In.NCI60',
}
LOGICAL_COLUMNS = [
	'In.AnnotationDB',
	'FDA.Approved',
	'In.JUMP.CP',
	'In.OASIS',
	'In.GEOM',
	'In.LINCS',
	'In.CTRP',
	'In.NCI60',
]
OUTPUT_COLUMN_ORDER = [
	'HDD.Compound.ID',
	'Pubchem.CID',
	'InChIKey',
	'SMILES',
	'Molecule.Name',
	'In.AnnotationDB',
	'AnnotationDB.Name',
	'AnnotationDB.SMILES',
	'In.JUMP.CP',
	'JUMP.CP.ID',
	'In.OASIS',
	'OASIS.ID',
	'In.GEOM',
	'GEOM.Source.SMILES',
	'In.LINCS',
	'LINCS.CMap.Name',
	'In.CTRP',
	'CTRP.Master.CPD.ID',
	'In.NCI60',
	'NCI60.NSC',
	'Molecular.Formula',
	'IUPAC.Name',
	'ChEMBL.ID',
	'ATC.Code',
	'PubChem.2D.Fingerprint',
	'Mechanism.of.Action',
	'FDA.Approved',
	'Molecular.Weight',
	'XlogP',
	'Hydrogen.Bond.Donors',
	'Hydrogen.Bond.Acceptors',
	'Exact.Molecular.Mass',
	'LTKB.DILI.Severity.Grade',
	'LTKB.DILI.Annotation',
	'DILIrank.DILI.Severity.Grade',
	'DILIrank.DILI.Annotation',
	'DILIst.DILI.Annotation',
	'Livertox.Hepatotoxicity.Likelihood.Score',
	'Livertox.Hepatotoxicity.Likelihood.Reason',
	'DIRIL.Label.Gong',
	'DIRIL.Label.Shi',
	'DIRIL.Toxicity',
	'DICT.Cardiotoxicity',
	'DICT.Label.Section',
	'DICT.Concern',
	'DICT.Keywords',
	'BBB.Permeable',
]
INTERNAL_COLUMNS = {'_Fallback.Identity'}


def normalize_text(value: object) -> str | None:
	if value is None or pd.isna(value):
		return None
	text = str(value).strip()
	if not text or text.upper() in {'NA', 'NAN', 'NONE', 'NULL'} or text == '-':
		return None
	return text


def normalize_cid(value: object) -> str | None:
	text = normalize_text(value)
	if text is None:
		return None

	try:
		numeric = float(text)
	except ValueError:
		return text

	if math.isfinite(numeric) and numeric.is_integer():
		return str(int(numeric))
	return text


def normalize_inchikey(value: object) -> str | None:
	text = normalize_text(value)
	if text is None:
		return None
	return text.upper()


def normalize_bool(value: object) -> object:
	if value is None or pd.isna(value):
		return pd.NA
	if isinstance(value, bool):
		return value
	text = str(value).strip().upper()
	if text in {'TRUE', 'T', '1', 'YES'}:
		return True
	if text in {'FALSE', 'F', '0', 'NO'}:
		return False
	return pd.NA


def normalize_dili_severity_grade(value: object) -> object:
	text = normalize_text(value)
	if text is None:
		return pd.NA
	normalized = text.upper().replace('_', ' ')
	if normalized in {'N/A', 'NOT APPLICABLE'}:
		return pd.NA
	return text


def first_non_missing(values: list[object]) -> object:
	for value in values:
		normalized = normalize_text(value)
		if normalized is not None:
			return normalized
	return pd.NA


def first_from_columns(row: pd.Series, columns: list[str]) -> object:
	return first_non_missing([row.get(column) for column in columns])


def format_joined_values(values: list[object]) -> object:
	normalized_values = sorted(
		{
			normalized
			for value in values
			if (normalized := normalize_text(value)) is not None
		}
	)
	if not normalized_values:
		return pd.NA
	return '|'.join(normalized_values)


def split_joined_values(value: object) -> set[str]:
	normalized = normalize_text(value)
	if normalized is None:
		return set()
	return {
		item.strip()
		for item in normalized.split('|')
		if item.strip() and item.strip().upper() != 'NA'
	}


def append_joined_value(record: dict, column: str, value: object) -> None:
	values = list(split_joined_values(record.get(column)))
	normalized = normalize_text(value)
	if normalized is not None:
		values.append(normalized)
	record[column] = format_joined_values(values)


def iter_records(path: str) -> Iterator[dict]:
	with Path(path).open('r', encoding='utf-8') as handle:
		for raw_line in handle:
			line = raw_line.strip()
			if not line:
				continue
			yield json.loads(line)


def normalize_toxicity_source(value: object) -> str | None:
	text = normalize_text(value)
	if text is None:
		return None
	lower = text.lower()
	if 'dilirank' in lower or 'drug induced liver injury rank' in lower:
		return 'DILIrank'
	if 'dilist' in lower or 'severity and toxicity' in lower:
		return 'DILIst'
	if 'livertox' in lower:
		return 'Livertox'
	if 'ltkb' in lower or 'benchmark dataset' in lower:
		return 'LTKB'
	return None


def iter_toxicity_records(toxicity: object) -> Iterator[dict]:
	if isinstance(toxicity, dict):
		records = [toxicity]
	elif isinstance(toxicity, list):
		records = toxicity
	else:
		records = []

	for record in records:
		if isinstance(record, dict):
			yield record


def flatten_toxicity_records(toxicity: object) -> dict:
	flattened = {column: pd.NA for column in TOXICITY_COLUMNS}
	for toxicity_record in iter_toxicity_records(toxicity):
		source = normalize_toxicity_source(toxicity_record.get('tox_dataset'))
		if source is None:
			continue
		for source_field, column in TOXICITY_SOURCE_COLUMNS[source].items():
			value = toxicity_record.get(source_field)
			if source_field == 'dili_severity_grade':
				value = normalize_dili_severity_grade(value)
			append_joined_value(flattened, column, value)
	return flattened


def build_annotationdb_record(drug_info: dict, drug_details: dict) -> dict:
	mechanisms = drug_details.get('mechanisms') or []
	toxicity = flatten_toxicity_records(drug_details.get('toxicity'))
	diril_toxicity = drug_details.get('diril_toxicity') or {}
	dict_rank_toxicity = drug_details.get('dict_rank_toxicity') or {}
	molecule_name = first_non_missing(
		[
			drug_info.get('name'),
			drug_info.get('mapped_name'),
			drug_details.get('title'),
		]
	)
	cid = normalize_cid(drug_info.get('cid') or drug_details.get('cid'))
	inchikey = normalize_inchikey(
		drug_info.get('inchikey') or drug_details.get('inchikey')
	)
	smiles = first_non_missing([drug_info.get('smiles'), drug_details.get('smiles')])
	mechanism = pd.NA
	if mechanisms:
		mechanism = first_non_missing(
			[
				mechanism_info.get('mechanism_of_action')
				for mechanism_info in mechanisms
				if isinstance(mechanism_info, dict)
			]
		)

	record = {
		'Pubchem.CID': cid,
		'InChIKey': inchikey,
		'SMILES': smiles,
		'Molecule.Name': molecule_name,
		'In.AnnotationDB': True,
		'AnnotationDB.Name': molecule_name,
		'AnnotationDB.SMILES': smiles,
		'Molecular.Formula': drug_details.get('molecular_formula'),
		'IUPAC.Name': drug_details.get('iupac_name'),
		'ChEMBL.ID': drug_details.get('molecule_chembl_id'),
		'ATC.Code': first_non_missing([drug_details.get('atc_code')]),
		'PubChem.2D.Fingerprint': drug_details.get('fingerprint_2d'),
		'Mechanism.of.Action': mechanism,
		'FDA.Approved': normalize_bool(drug_details.get('fda_approval')),
		'Molecular.Weight': drug_details.get('molecular_weight'),
		'XlogP': drug_details.get('xlogp'),
		'Hydrogen.Bond.Donors': drug_details.get('h_bond_donor_count'),
		'Hydrogen.Bond.Acceptors': drug_details.get('h_bond_acceptor_count'),
		'Exact.Molecular.Mass': drug_details.get('exact_mass'),
		**toxicity,
		'DIRIL.Label.Gong': diril_toxicity.get('label_gong'),
		'DIRIL.Label.Shi': diril_toxicity.get('label_shi'),
		'DIRIL.Toxicity': diril_toxicity.get('toxicity'),
		'DICT.Cardiotoxicity': dict_rank_toxicity.get('cardiotoxicity'),
		'DICT.Label.Section': dict_rank_toxicity.get('label_section'),
		'DICT.Concern': dict_rank_toxicity.get('dict_concern'),
		'DICT.Keywords': dict_rank_toxicity.get('keywords')
		or dict_rank_toxicity.get('keyword'),
		'BBB.Permeable': pd.NA,
	}
	for flag in MEMBERSHIP_COLUMNS.values():
		record[flag] = False
	for source_column in SOURCE_COLUMNS.values():
		record[source_column] = pd.NA
	record['_Fallback.Identity'] = pd.NA
	return record


def load_annotationdb_records(input_path: str) -> tuple[list[dict], dict[str, list]]:
	records: list[dict] = []
	bioassays_by_cid: dict[str, list] = {}
	error_cids: list[str] = []

	for raw_record in tqdm.tqdm(iter_records(input_path), desc='AnnotationDB'):
		drug_info = raw_record.get('drug_info')
		drug_details = raw_record.get('drug_details')
		if not isinstance(drug_info, dict) or not isinstance(drug_details, dict):
			continue

		try:
			record = build_annotationdb_record(drug_info, drug_details)
		except Exception:
			cid = drug_info.get('cid')
			if cid is not None:
				error_cids.append(str(cid))
			continue

		cid = normalize_cid(record.get('Pubchem.CID'))
		if cid is None:
			continue
		records.append(record)
		bioassays_by_cid[cid] = drug_details.get('bioassays') or []

	if error_cids:
		print(  # noqa: T201
			'[process_annotationdb] annotationdb_processing_failures='
			f'{len(error_cids)}',
			flush=True,
		)
	return records, bioassays_by_cid


def build_indexes(records: list[dict]) -> tuple[dict[str, int], dict[str, int]]:
	cid_to_index: dict[str, int] = {}
	inchikey_to_index: dict[str, int] = {}
	for index, record in enumerate(records):
		cid = normalize_cid(record.get('Pubchem.CID'))
		inchikey = normalize_inchikey(record.get('InChIKey'))
		if cid is not None and cid not in cid_to_index:
			cid_to_index[cid] = index
		if inchikey is not None and inchikey not in inchikey_to_index:
			inchikey_to_index[inchikey] = index
	return cid_to_index, inchikey_to_index


def update_indexes_for_record(
	index: int,
	record: dict,
	cid_to_index: dict[str, int],
	inchikey_to_index: dict[str, int],
) -> None:
	cid = normalize_cid(record.get('Pubchem.CID'))
	inchikey = normalize_inchikey(record.get('InChIKey'))
	if cid is not None and cid not in cid_to_index:
		cid_to_index[cid] = index
	if inchikey is not None and inchikey not in inchikey_to_index:
		inchikey_to_index[inchikey] = index


def source_smiles(row: pd.Series, spec: dict) -> object:
	return first_from_columns(row, spec.get('source_smiles_columns', []))


def source_name(row: pd.Series, spec: dict) -> object:
	return first_from_columns(row, spec.get('source_name_columns', []))


def source_identity(
	dataset: str,
	source_key: object,
	row: pd.Series,
) -> str:
	identity_value = normalize_text(source_key) or normalize_text(
		row.get('RDS.RowName')
	)
	if identity_value is None:
		identity_value = str(row.name)
	return f'source:{dataset}:{identity_value}'


def create_source_record(dataset: str, row: pd.Series, spec: dict) -> dict:
	source_key_column = spec['source_key_column']
	source_key = row.get(source_key_column)
	smiles = source_smiles(row, spec)
	molecule_name = source_name(row, spec)
	record = {
		'Pubchem.CID': normalize_cid(row.get('Pubchem.CID')),
		'InChIKey': normalize_inchikey(row.get('InChIKey')),
		'SMILES': smiles,
		'Molecule.Name': first_non_missing([molecule_name, source_key]),
		'In.AnnotationDB': False,
		'AnnotationDB.Name': pd.NA,
		'AnnotationDB.SMILES': pd.NA,
		'Molecular.Formula': pd.NA,
		'IUPAC.Name': pd.NA,
		'ChEMBL.ID': pd.NA,
		'ATC.Code': pd.NA,
		'PubChem.2D.Fingerprint': pd.NA,
		'Mechanism.of.Action': pd.NA,
		'FDA.Approved': pd.NA,
		'Molecular.Weight': pd.NA,
		'XlogP': pd.NA,
		'Hydrogen.Bond.Donors': pd.NA,
		'Hydrogen.Bond.Acceptors': pd.NA,
		'Exact.Molecular.Mass': pd.NA,
		**{column: pd.NA for column in TOXICITY_COLUMNS},
		'DIRIL.Label.Gong': pd.NA,
		'DIRIL.Label.Shi': pd.NA,
		'DIRIL.Toxicity': pd.NA,
		'DICT.Cardiotoxicity': pd.NA,
		'DICT.Label.Section': pd.NA,
		'DICT.Concern': pd.NA,
		'DICT.Keywords': pd.NA,
		'BBB.Permeable': pd.NA,
		'_Fallback.Identity': source_identity(dataset, source_key, row),
	}
	for flag in MEMBERSHIP_COLUMNS.values():
		record[flag] = False
	for source_column in SOURCE_COLUMNS.values():
		record[source_column] = pd.NA
	return record


def find_matching_record_index(
	row: pd.Series,
	cid_to_index: dict[str, int],
	inchikey_to_index: dict[str, int],
) -> int | None:
	cid = normalize_cid(row.get('Pubchem.CID'))
	if cid is not None and cid in cid_to_index:
		return cid_to_index[cid]

	inchikey = normalize_inchikey(row.get('InChIKey'))
	if inchikey is not None and inchikey in inchikey_to_index:
		return inchikey_to_index[inchikey]

	return None


def update_record_from_source(
	record: dict,
	dataset: str,
	row: pd.Series,
	spec: dict,
) -> None:
	flag_column = spec['flag_column']
	source_key_column = spec['source_key_column']
	record[flag_column] = True
	append_joined_value(record, source_key_column, row.get(source_key_column))

	if normalize_text(record.get('SMILES')) is None:
		record['SMILES'] = source_smiles(row, spec)
	if normalize_text(record.get('Molecule.Name')) is None:
		record['Molecule.Name'] = first_non_missing(
			[source_name(row, spec), row.get(source_key_column)]
		)
	if normalize_text(record.get('Pubchem.CID')) is None:
		record['Pubchem.CID'] = normalize_cid(row.get('Pubchem.CID'))
	if normalize_text(record.get('InChIKey')) is None:
		record['InChIKey'] = normalize_inchikey(row.get('InChIKey'))

	if not record.get('In.AnnotationDB', False):
		record['AnnotationDB.Name'] = pd.NA
		record['AnnotationDB.SMILES'] = pd.NA
	elif normalize_text(record.get('AnnotationDB.SMILES')) is None:
		record['AnnotationDB.SMILES'] = row.get('AnnotationDB.SMILES')


def load_sub_dataset_metadata(
	metadata_paths: list[str],
	dataset_names: list[str],
) -> dict[str, pd.DataFrame]:
	if len(metadata_paths) != len(dataset_names):
		message = 'sub_dataset metadata paths and names have different lengths'
		raise ValueError(message)

	frames = {}
	for dataset, path in zip(dataset_names, metadata_paths, strict=True):
		frames[dataset] = pd.read_csv(path, sep='\t', dtype='string')
	return frames


def merge_sub_dataset_records(
	records: list[dict],
	sub_dataset_frames: dict[str, pd.DataFrame],
	sub_dataset_specs: dict,
) -> None:
	cid_to_index, inchikey_to_index = build_indexes(records)

	for dataset, frame in sub_dataset_frames.items():
		spec = sub_dataset_specs[dataset]
		source_key_column = spec['source_key_column']
		if source_key_column not in frame.columns:
			message = (
				f'Sub-dataset {dataset} metadata is missing source key column '
				f'{source_key_column}'
			)
			raise ValueError(message)

		for _, row in frame.iterrows():
			index = find_matching_record_index(row, cid_to_index, inchikey_to_index)
			if index is None:
				record = create_source_record(dataset, row, spec)
				records.append(record)
				index = len(records) - 1
				update_indexes_for_record(
					index, record, cid_to_index, inchikey_to_index
				)

			update_record_from_source(records[index], dataset, row, spec)
			update_indexes_for_record(
				index, records[index], cid_to_index, inchikey_to_index
			)


def base62_encode_digest(digest: bytes) -> str:
	number = int.from_bytes(digest, 'big')
	if number == 0:
		return BASE62_ALPHABET[0]
	encoded = []
	while number:
		number, remainder = divmod(number, len(BASE62_ALPHABET))
		encoded.append(BASE62_ALPHABET[remainder])
	return ''.join(reversed(encoded))


def identity_for_record(record: dict) -> str:
	cid = normalize_cid(record.get('Pubchem.CID'))
	if cid is not None:
		return f'pubchem:{cid}'

	inchikey = normalize_inchikey(record.get('InChIKey'))
	if inchikey is not None:
		return f'inchikey:{inchikey}'

	fallback = normalize_text(record.get('_Fallback.Identity'))
	if fallback is not None:
		return fallback

	message = 'Cannot build HDD.Compound.ID for record without any identity'
	raise ValueError(message)


def compact_hash(identity: str, length: int) -> str:
	digest = hashlib.sha256(identity.encode('utf-8')).digest()
	encoded = base62_encode_digest(digest)
	if len(encoded) < length:
		encoded = encoded.rjust(length, BASE62_ALPHABET[0])
	return f'HDD_{encoded[:length]}'


def assign_hdd_ids(records: list[dict]) -> None:
	identities = [identity_for_record(record) for record in records]
	length_by_index = [DEFAULT_HASH_LENGTH for _ in identities]

	while True:
		ids = [
			compact_hash(identity, length)
			for identity, length in zip(identities, length_by_index, strict=True)
		]
		duplicates = {
			hdd_id
			for hdd_id, count in pd.Series(ids).value_counts().items()
			if count > 1
		}
		if not duplicates:
			break

		for index, hdd_id in enumerate(ids):
			if hdd_id in duplicates:
				length_by_index[index] += HASH_LENGTH_INCREMENT
				if length_by_index[index] > MAX_HASH_LENGTH:
					message = 'Unable to resolve HDD.Compound.ID hash collision'
					raise ValueError(message)

	for record, hdd_id in zip(records, ids, strict=True):
		record['HDD.Compound.ID'] = hdd_id


def add_bbbp_annotations(records: list[dict], bbbp_file: str) -> None:
	bbbp_path = Path(bbbp_file)
	if not bbbp_path.exists():
		return

	bbbp = pd.read_csv(bbbp_path)
	if 'smiles' not in bbbp.columns or 'p_np' not in bbbp.columns:
		return

	bbbp_map = (
		bbbp.dropna(subset=['smiles'])
		.drop_duplicates(subset=['smiles'], keep='first')
		.set_index('smiles')['p_np']
		.to_dict()
	)
	for record in records:
		smiles = normalize_text(record.get('SMILES'))
		if smiles is not None and smiles in bbbp_map:
			record['BBB.Permeable'] = bbbp_map[smiles]


def records_to_coldata(records: list[dict]) -> pd.DataFrame:
	assign_hdd_ids(records)
	coldata = pd.DataFrame(records)
	coldata = coldata.drop(columns=[col for col in INTERNAL_COLUMNS if col in coldata])

	for column in LOGICAL_COLUMNS:
		if column in coldata.columns:
			coldata[column] = coldata[column].map(normalize_bool).astype('boolean')

	for column in SOURCE_COLUMNS.values():
		if column in coldata.columns:
			coldata[column] = coldata[column].map(
				lambda value: (
					pd.NA if normalize_text(value) is None else normalize_text(value)
				)
			)

	ordered = [column for column in OUTPUT_COLUMN_ORDER if column in coldata.columns]
	remaining = [column for column in coldata.columns if column not in ordered]
	coldata = coldata[ordered + remaining]
	coldata = coldata.sort_values('HDD.Compound.ID').reset_index(drop=True)

	if coldata['HDD.Compound.ID'].isna().any():
		message = 'HDD.Compound.ID contains missing values'
		raise ValueError(message)
	if coldata['HDD.Compound.ID'].duplicated().any():
		message = 'HDD.Compound.ID contains duplicate values'
		raise ValueError(message)

	return coldata


def write_bioassay_matrix(
	bioassays_by_cid: dict[str, list],
	coldata: pd.DataFrame,
	bioassays_path: Path,
) -> int:
	seen_assays = sorted(
		{
			assay.get('aid')
			for assays in bioassays_by_cid.values()
			for assay in assays
			if isinstance(assay, dict) and assay.get('aid') is not None
		}
	)
	aid_to_idx = {aid: idx for idx, aid in enumerate(seen_assays)}
	matrix = {}

	for _, row in coldata.iterrows():
		cid = normalize_cid(row.get('Pubchem.CID'))
		hdd_id = row['HDD.Compound.ID']
		if cid is None or cid not in bioassays_by_cid:
			continue

		outcomes = ['Not Measured'] * len(seen_assays)
		for assay in bioassays_by_cid[cid]:
			assay_id = assay.get('aid')
			if assay_id not in aid_to_idx:
				continue
			outcome = (
				'Active'
				if assay.get('activity_outcome_method') == ACTIVE_OUTCOME_METHOD
				else 'Inactive'
			)
			outcomes[aid_to_idx[assay_id]] = outcome
		matrix[hdd_id] = outcomes

	bioassay_res = pd.DataFrame(matrix, index=[f'AID_{aid}' for aid in seen_assays])
	bioassay_res = bioassay_res.reset_index(names='Assay')
	bioassay_res.to_csv(bioassays_path, index=False)
	return len(matrix)


def expected_source_keys(frame: pd.DataFrame, source_key_column: str) -> set[str]:
	if source_key_column not in frame.columns:
		return set()
	return {
		normalized
		for value in frame[source_key_column]
		if (normalized := normalize_text(value)) is not None
	}


def actual_source_keys(
	coldata: pd.DataFrame,
	flag_column: str,
	source_key_column: str,
) -> set[str]:
	if flag_column not in coldata.columns or source_key_column not in coldata.columns:
		return set()
	flagged = coldata[coldata[flag_column].fillna(False).astype(bool)]
	values: set[str] = set()
	for value in flagged[source_key_column]:
		values.update(split_joined_values(value))
	return values


def write_parity_report(
	coldata: pd.DataFrame,
	sub_dataset_frames: dict[str, pd.DataFrame],
	sub_dataset_specs: dict,
	parity_dir: Path,
) -> None:
	if parity_dir.exists():
		shutil.rmtree(parity_dir)
	parity_dir.mkdir(parents=True, exist_ok=True)

	summary_rows = []
	for dataset, frame in sub_dataset_frames.items():
		spec = sub_dataset_specs[dataset]
		source_key_column = spec['source_key_column']
		flag_column = spec['flag_column']
		expected = expected_source_keys(frame, source_key_column)
		actual = actual_source_keys(coldata, flag_column, source_key_column)
		missing_in_hdd = sorted(expected - actual)
		extra_in_hdd = sorted(actual - expected)
		missing_source_id = frame[
			frame[source_key_column].map(normalize_text).isna()
		].copy()

		details = pd.DataFrame(
			[
				{'status': 'missing_in_hdd', source_key_column: value}
				for value in missing_in_hdd
			]
			+ [
				{'status': 'extra_in_hdd', source_key_column: value}
				for value in extra_in_hdd
			]
		)
		if details.empty:
			details = pd.DataFrame(columns=['status', source_key_column])
		details.to_csv(
			parity_dir / f'{dataset}_source_key_parity.tsv', sep='\t', index=False
		)

		missing_source_id.to_csv(
			parity_dir / f'{dataset}_missing_source_key_rows.tsv',
			sep='\t',
			index=False,
		)
		summary_rows.append(
			{
				'dataset': dataset,
				'source_key_column': source_key_column,
				'rds_source_keys': len(expected),
				'hdd_source_keys': len(actual),
				'missing_in_hdd': len(missing_in_hdd),
				'extra_in_hdd': len(extra_in_hdd),
				'missing_source_key_rows': len(missing_source_id),
				'passes': len(missing_in_hdd) == 0 and len(extra_in_hdd) == 0,
			}
		)

	summary = pd.DataFrame(summary_rows)
	summary.to_csv(parity_dir / 'summary.tsv', sep='\t', index=False)
	if not summary['passes'].all():
		failing = ', '.join(summary.loc[~summary['passes'], 'dataset'])
		message = f'Sub-dataset source-key parity failed for: {failing}'
		raise ValueError(message)


def main(
	input_path: str,
	sub_dataset_metadata_paths: list[str],
	sub_dataset_names: list[str],
	sub_dataset_specs: dict,
	bbbp_file: str,
	coldata_output: str,
	bioassays_output: str,
	parity_output: str,
) -> None:
	coldata_path = Path(coldata_output)
	bioassays_path = Path(bioassays_output)
	parity_dir = Path(parity_output)
	coldata_path.parent.mkdir(parents=True, exist_ok=True)
	bioassays_path.parent.mkdir(parents=True, exist_ok=True)

	records, bioassays_by_cid = load_annotationdb_records(input_path)
	sub_dataset_frames = load_sub_dataset_metadata(
		sub_dataset_metadata_paths,
		sub_dataset_names,
	)
	merge_sub_dataset_records(records, sub_dataset_frames, sub_dataset_specs)
	add_bbbp_annotations(records, bbbp_file)
	coldata = records_to_coldata(records)

	coldata.to_csv(coldata_path, index=False)
	bioassay_columns = write_bioassay_matrix(
		bioassays_by_cid,
		coldata,
		bioassays_path,
	)
	write_parity_report(coldata, sub_dataset_frames, sub_dataset_specs, parity_dir)

	print(  # noqa: T201
		'[process_annotationdb] '
		f'coldata_rows={len(coldata)} '
		f'annotationdb_rows={int(coldata["In.AnnotationDB"].sum())} '
		f'bioassay_columns={bioassay_columns} '
		f'output={coldata_path}',
		flush=True,
	)


def main_from_snakemake() -> None:
	main(
		input_path=str(snakemake.input.raw_data),
		sub_dataset_metadata_paths=[
			str(path) for path in snakemake.input.sub_dataset_metadata
		],
		sub_dataset_names=list(snakemake.params.sub_dataset_names),
		sub_dataset_specs=json.loads(snakemake.params.sub_dataset_specs),
		bbbp_file=str(snakemake.input.bbbp_file),
		coldata_output=str(snakemake.output.colData),
		bioassays_output=str(snakemake.output.bioassays),
		parity_output=str(snakemake.output.parity),
	)


if __name__ == '__main__':
	if 'snakemake' in globals():
		main_from_snakemake()
	else:
		parser = argparse.ArgumentParser(
			prog='process_annotationdb',
			description='Generate HDD colData and bioassays from AnnotationDB and sub-dataset MAEs',
		)
		parser.add_argument(
			'-i', required=True, help='Input JSONL from fetch_annotationdb'
		)
		parser.add_argument(
			'--sub-dataset-metadata',
			nargs='+',
			required=True,
			help='Extracted sub-dataset Drug.Metadata TSV files',
		)
		parser.add_argument(
			'--sub-dataset-names',
			nargs='+',
			required=True,
			help='Dataset names matching --sub-dataset-metadata order',
		)
		parser.add_argument(
			'--sub-dataset-specs',
			required=True,
			help='JSON encoded sub_dataset config object',
		)
		parser.add_argument('-b', required=True, help='Blood brain barrier CSV')
		parser.add_argument('-c', required=True, help='Output colData CSV')
		parser.add_argument('-a', required=True, help='Output bioassays CSV')
		parser.add_argument('-p', required=True, help='Output parity report directory')
		args = parser.parse_args()

		main(
			input_path=args.i,
			sub_dataset_metadata_paths=args.sub_dataset_metadata,
			sub_dataset_names=args.sub_dataset_names,
			sub_dataset_specs=json.loads(args.sub_dataset_specs),
			bbbp_file=args.b,
			coldata_output=args.c,
			bioassays_output=args.a,
			parity_output=args.p,
		)
