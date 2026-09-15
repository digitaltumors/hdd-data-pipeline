import argparse
import hashlib
import json
import threading
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Sequence, Tuple
from urllib.parse import urlparse

import requests
import tqdm

DEFAULT_BATCH_SIZE = 250
DEFAULT_WORKERS = 6
DEFAULT_QUERY_DELAY_SECONDS = 0.1
DEFAULT_ENV_FILE = '.env'
DEFAULT_API_KEY_ENV_VAR = 'ANNOTATIONDB_API_KEY'
MIN_QUOTED_VALUE_LENGTH = 2
MAX_QUERY_SIZE = 250
RETRIES = 3
TIMEOUT = 60
TOXICITY_FIELDS = [
	'tox_dataset',
	'dili_severity_grade',
	'dili_annotation',
	'hepatotoxicity_likelihood_score',
	'hepatotoxicity_likelihood_score_reasoning',
]
SPLIT_STATUS_CODES = {413, 422, 429, 500, 502, 503, 504}
SESSION_HEADERS = {
	'accept': 'application/json',
	'user-agent': 'hdd-data-pipeline/1.0',
}

thread_local = threading.local()


def read_env_value(env_file: str | Path, variable: str) -> Optional[str]:
	path = Path(env_file)
	if not path.is_file():
		return None

	for raw_line in path.read_text(encoding='utf-8').splitlines():
		line = raw_line.strip()
		if not line or line.startswith('#'):
			continue
		if line.startswith('export '):
			line = line.removeprefix('export ').lstrip()

		key, separator, value = line.partition('=')
		if not separator or key.strip() != variable:
			continue

		value = value.strip()
		if (
			len(value) >= MIN_QUOTED_VALUE_LENGTH
			and value[0] == value[-1]
			and value[0] in {'"', "'"}
		):
			value = value[1:-1]
		return value or None
	return None


def configure_authorization(api_key_file: str | Path, api_key_name: str) -> bool:
	api_key = read_env_value(api_key_file, api_key_name)
	SESSION_HEADERS.pop('X-API-Key', None)
	if api_key:
		SESSION_HEADERS['X-API-Key'] = api_key

	thread_local.session = None
	return bool(api_key)


def fetch_json(
	session: requests.Session,
	url: str,
	params: Optional[Sequence[Tuple[str, str]]] = None,
	retries: int = RETRIES,
	timeout: int = TIMEOUT,
) -> object:
	for attempt in range(retries):
		try:
			resp = session.get(url, params=params, timeout=timeout)
			resp.raise_for_status()
			try:
				return resp.json()
			except ValueError as exc:
				preview = resp.text[:200].replace('\n', ' ')
				message = (
					f'Expected JSON response from {resp.url}, '
					f'got content-type={resp.headers.get("content-type")!r}: {preview}'
				)
				raise ValueError(message) from exc
		except Exception:
			if attempt == retries - 1:
				raise
			time.sleep(2 * (attempt + 1))


def normalize_details_response(data: object) -> List[dict]:
	if isinstance(data, list):
		return data

	if isinstance(data, dict) and 'compounds' in data and 'bioassays' in data:
		compounds = data.get('compounds') or []
		bioassays = data.get('bioassays') or {}
		if not isinstance(compounds, list) or not isinstance(bioassays, dict):
			message = 'Unexpected streamline response shape from AnnotationDB'
			raise TypeError(message)

		normalized = []
		for compound in compounds:
			if not isinstance(compound, dict):
				continue

			compound_bioassays = []
			for aid in compound.get('bioassays') or []:
				assay = bioassays.get(str(aid), bioassays.get(aid))
				if isinstance(assay, dict):
					compound_bioassays.append(assay)

			normalized.append({**compound, 'bioassays': compound_bioassays})

		return normalized

	message = 'Expected list or streamline response from AnnotationDB /compound/many'
	raise TypeError(message)


def build_details_params(
	cids: List[str],
	golden_bioassay: bool,
	indication: bool,
) -> List[Tuple[str, str]]:
	params = [('compound', cid) for cid in cids]
	params.extend(
		[
			('format', 'json'),
			('bioassay', 'true'),
			('mechanism', 'true'),
			('toxicity', 'true'),
		]
	)
	if golden_bioassay:
		params.append(('golden_bioassay', 'true'))
	if indication:
		params.append(('indication', 'true'))
	return params


def fetch_split_batch(  # noqa: PLR0917
	session: requests.Session,
	details_url: str,
	cids: List[str],
	failed_cids: List[str],
	query_delay_seconds: float,
	golden_bioassay: bool,
	indication: bool,
) -> List[dict]:
	mid = len(cids) // 2
	left = fetch_details_for_cids(
		session,
		details_url,
		cids[:mid],
		failed_cids,
		query_delay_seconds,
		golden_bioassay,
		indication,
	)
	right = fetch_details_for_cids(
		session,
		details_url,
		cids[mid:],
		failed_cids,
		query_delay_seconds,
		golden_bioassay,
		indication,
	)
	return left + right


def fetch_chunked_details(  # noqa: PLR0917
	session: requests.Session,
	details_url: str,
	cids: List[str],
	failed_cids: List[str],
	query_delay_seconds: float,
	golden_bioassay: bool,
	indication: bool,
) -> List[dict]:
	results = []
	for start in range(0, len(cids), MAX_QUERY_SIZE):
		results.extend(
			fetch_details_for_cids(
				session,
				details_url,
				cids[start : start + MAX_QUERY_SIZE],
				failed_cids,
				query_delay_seconds,
				golden_bioassay,
				indication,
			)
		)
	return results


def record_failed_cids(failed_cids: List[str], cids: List[str]) -> List[dict]:
	failed_cids.extend(cids)
	return []


def resolve_failed_fetch(  # noqa: PLR0917
	session: requests.Session,
	details_url: str,
	cids: List[str],
	failed_cids: List[str],
	query_delay_seconds: float,
	golden_bioassay: bool,
	attempt: int,
	allow_split: bool,
	indication: bool,
) -> Optional[List[dict]]:
	if allow_split and len(cids) > 1:
		return fetch_split_batch(
			session,
			details_url,
			cids,
			failed_cids,
			query_delay_seconds,
			golden_bioassay,
			indication,
		)
	if attempt == RETRIES - 1:
		return record_failed_cids(failed_cids, cids)
	return None


def raise_for_split_status(resp: requests.Response) -> None:
	if resp.status_code in SPLIT_STATUS_CODES:
		message = f'AnnotationDB responded with retryable status {resp.status_code}'
		raise requests.HTTPError(message, response=resp)


def emit_status(message: str) -> None:
	tqdm.tqdm.write(message)


def fetch_details_for_cids(  # noqa: PLR0917
	session: requests.Session,
	details_url: str,
	cids: List[str],
	failed_cids: List[str],
	query_delay_seconds: float,
	golden_bioassay: bool,
	indication: bool,
) -> List[dict]:
	if not cids:
		return []

	if len(cids) > MAX_QUERY_SIZE:
		return fetch_chunked_details(
			session,
			details_url,
			cids,
			failed_cids,
			query_delay_seconds,
			golden_bioassay,
			indication,
		)

	params = build_details_params(cids, golden_bioassay, indication)

	for attempt in range(RETRIES):
		details: Optional[List[dict]] = None
		try:
			resp = session.get(details_url, params=params, timeout=TIMEOUT)
			raise_for_split_status(resp)
			resp.raise_for_status()
			details = normalize_details_response(resp.json())
		except requests.HTTPError as exc:
			status = exc.response.status_code if exc.response is not None else None
			details = resolve_failed_fetch(
				session,
				details_url,
				cids,
				failed_cids,
				query_delay_seconds,
				golden_bioassay,
				attempt,
				allow_split=status in SPLIT_STATUS_CODES,
				indication=indication,
			)
		except (requests.RequestException, TypeError, ValueError):
			details = resolve_failed_fetch(
				session,
				details_url,
				cids,
				failed_cids,
				query_delay_seconds,
				golden_bioassay,
				attempt,
				allow_split=True,
				indication=indication,
			)
		finally:
			if query_delay_seconds > 0:
				time.sleep(query_delay_seconds)

		if details is not None:
			return details

		time.sleep(2 * (attempt + 1))

	return record_failed_cids(failed_cids, cids)


def derive_details_url(db_url: str) -> str:
	parsed = urlparse(db_url)
	path = parsed.path
	if path.endswith('/compound/all'):
		path = path.replace('/compound/all', '/compound/many')
	else:
		path = path.rstrip('/') + '/many'
	return parsed._replace(path=path, query='').geturl()


def chunks(items: List[dict], size: int) -> Iterator[List[dict]]:
	for idx in range(0, len(items), size):
		yield items[idx : idx + size]


def normalize_cid(value: object) -> Optional[str]:
	if value is None:
		return None
	return str(value)


def compact_toxicity_records(toxicity: object) -> List[dict]:
	if isinstance(toxicity, dict):
		records = [toxicity]
	elif isinstance(toxicity, list):
		records = toxicity
	else:
		records = []

	return [
		{field: record.get(field) for field in TOXICITY_FIELDS}
		for record in records
		if isinstance(record, dict)
	]


def slim_detail_record(detail: dict) -> dict:
	diril_toxicity = detail.get('diril_toxicity') or {}
	dict_rank_toxicity = detail.get('dict_rank_toxicity') or {}
	mechanisms = detail.get('mechanisms') or []
	bioassays = detail.get('bioassays') or []
	drug_indications = detail.get('drug_indications') or []
	return {
		'cid': detail.get('cid'),
		'molecular_formula': detail.get('molecular_formula'),
		'iupac_name': detail.get('iupac_name'),
		'molecule_chembl_id': detail.get('molecule_chembl_id'),
		'atc_code': detail.get('atc_code'),
		'fingerprint_2d': detail.get('fingerprint_2d'),
		'mechanisms': [
			{'mechanism_of_action': mechanism.get('mechanism_of_action')}
			for mechanism in mechanisms
			if isinstance(mechanism, dict)
			and mechanism.get('mechanism_of_action') is not None
		],
		'fda_approval': detail.get('fda_approval'),
		'molecular_weight': detail.get('molecular_weight'),
		'xlogp': detail.get('xlogp'),
		'h_bond_donor_count': detail.get('h_bond_donor_count'),
		'h_bond_acceptor_count': detail.get('h_bond_acceptor_count'),
		'exact_mass': detail.get('exact_mass'),
		'toxicity': compact_toxicity_records(detail.get('toxicity')),
		'diril_toxicity': {
			'label_gong': diril_toxicity.get('label_gong'),
			'label_shi': diril_toxicity.get('label_shi'),
			'toxicity': diril_toxicity.get('toxicity'),
		},
		'dict_rank_toxicity': {
			'cardiotoxicity': dict_rank_toxicity.get('cardiotoxicity'),
			'label_section': dict_rank_toxicity.get('label_section'),
			'dict_concern': dict_rank_toxicity.get('dict_concern'),
			'keywords': dict_rank_toxicity.get('keywords')
			or dict_rank_toxicity.get('keyword'),
		},
		'bioassays': [
			{
				'aid': assay.get('aid'),
				'assay_name': assay.get('assay_name'),
				'source_name': assay.get('source_name'),
				'source_id': assay.get('source_id'),
				'activity_outcome_method': assay.get('activity_outcome_method'),
				'target_name': assay.get('target_name'),
				'target_protein_accession': assay.get('target_protein_accession'),
			}
			for assay in bioassays
			if isinstance(assay, dict) and assay.get('aid') is not None
		],
		'drug_indications': [
			{
				'molecule_chembl_id': indication.get('molecule_chembl_id'),
				'parent_molecule_chembl_id': indication.get(
					'parent_molecule_chembl_id'
				),
				'child_molecule_chembl_id': indication.get('child_molecule_chembl_id'),
				'drugind_id': indication.get('drugind_id'),
				'max_phase_for_ind': indication.get('max_phase_for_ind'),
				'mesh_id': indication.get('mesh_id'),
				'mesh_heading': indication.get('mesh_heading'),
				'efo_id': indication.get('efo_id'),
				'efo_term': indication.get('efo_term'),
				'clinical_trials_ref_ids': indication.get('clinical_trials_ref_ids'),
				'daily_med_ref_ids': indication.get('daily_med_ref_ids'),
				'ema_ref_ids': indication.get('ema_ref_ids'),
				'fda_ref_ids': indication.get('fda_ref_ids'),
				'usan_ref_ids': indication.get('usan_ref_ids'),
				'inn_ref_ids': indication.get('inn_ref_ids'),
				'inferred_from_parent': indication.get('inferred_from_parent'),
				'inferred_from_child': indication.get('inferred_from_child'),
			}
			for indication in drug_indications
			if isinstance(indication, dict)
		],
	}


def get_thread_session() -> requests.Session:
	session = getattr(thread_local, 'session', None)
	if session is None:
		session = requests.Session()
		session.headers.update(SESSION_HEADERS)
		thread_local.session = session
	return session


def fetch_batch_records(
	details_url: str,
	batch: List[dict],
	query_delay_seconds: float,
	golden_bioassay: bool,
	indication: bool,
) -> Tuple[List[dict], List[str], List[str]]:
	session = get_thread_session()
	failed_cids: List[str] = []
	missing_cids: List[str] = []

	cids = [
		normalize_cid(item.get('cid'))
		for item in batch
		if normalize_cid(item.get('cid')) is not None
	]
	details = fetch_details_for_cids(
		session,
		details_url,
		cids,
		failed_cids,
		query_delay_seconds,
		golden_bioassay,
		indication,
	)

	details_by_cid = {
		normalize_cid(item.get('cid')): slim_detail_record(item)
		for item in details
		if normalize_cid(item.get('cid')) is not None
	}

	records: List[dict] = []
	for drug_info in batch:
		cid = normalize_cid(drug_info.get('cid'))
		if cid is None:
			continue
		drug_details = details_by_cid.get(cid)
		if drug_details is None:
			missing_cids.append(cid)
			continue
		records.append({'drug_info': drug_info, 'drug_details': drug_details})

	return records, missing_cids, failed_cids


def iter_parallel_fetches(  # noqa: PLR0917
	details_url: str,
	batches: Iterable[List[dict]],
	workers: int,
	query_delay_seconds: float,
	golden_bioassay: bool,
	indication: bool,
) -> Iterator[Tuple[List[dict], List[str], List[str]]]:
	batch_iter = iter(batches)
	max_pending = max(workers * 2, 1)

	with ThreadPoolExecutor(max_workers=workers) as executor:
		pending = set()
		while True:
			while len(pending) < max_pending:
				try:
					batch = next(batch_iter)
				except StopIteration:
					break
				pending.add(
					executor.submit(
						fetch_batch_records,
						details_url,
						batch,
						query_delay_seconds,
						golden_bioassay,
						indication,
					)
				)

			if not pending:
				break

			done, pending = wait(pending, return_when=FIRST_COMPLETED)
			for future in done:
				yield future.result()


def main(  # noqa: PLR0912, PLR0915, PLR0917
	db_url: str,
	output_path: str,
	manifest_path: str | None = None,
	details_url: Optional[str] = None,
	batch_size: int = DEFAULT_BATCH_SIZE,
	workers: int = DEFAULT_WORKERS,
	query_delay_seconds: float = DEFAULT_QUERY_DELAY_SECONDS,
	golden_bioassay: bool = True,
	indication: bool = True,
	limit: Optional[int] = None,
	api_key_file: str | Path = DEFAULT_ENV_FILE,
	api_key_name: str = DEFAULT_API_KEY_ENV_VAR,
) -> None:
	if batch_size < 1 or batch_size > MAX_QUERY_SIZE:
		message = f'batch_size must be between 1 and {MAX_QUERY_SIZE}'
		raise ValueError(message)
	if workers < 1:
		message = 'workers must be >= 1'
		raise ValueError(message)
	if query_delay_seconds < 0:
		message = 'query_delay_seconds must be >= 0'
		raise ValueError(message)
	if limit is not None and limit < 1:
		message = 'limit must be >= 1'
		raise ValueError(message)

	has_authorization = configure_authorization(api_key_file, api_key_name)
	if has_authorization:
		emit_status('AnnotationDB authorization key detected; requesting ATC data')
	else:
		emit_status(
			f'AnnotationDB authorization key {api_key_name} is not set; '
			'ATC data will be excluded'
		)

	outpath = Path(output_path)
	outpath.parent.mkdir(parents=True, exist_ok=True)

	session = requests.Session()
	session.headers.update(SESSION_HEADERS)
	compound_list = fetch_json(session, db_url)
	if not isinstance(compound_list, list):
		message = 'Expected list response from /compound/all'
		raise TypeError(message)
	if limit is not None:
		compound_list = compound_list[:limit]

	details_url = details_url or derive_details_url(db_url)
	missing_cids: List[str] = []
	failed_cids: List[str] = []
	written = 0

	total_batches = (len(compound_list) + batch_size - 1) // batch_size
	emit_status(
		'Fetching AnnotationDB details: '
		f'{len(compound_list)} compounds in {total_batches} batches '
		f'(batch_size={batch_size}, workers={workers}, '
		f'golden_bioassay={golden_bioassay}, details_url={details_url})'
	)
	temporary_path = outpath.with_suffix(outpath.suffix + '.partial')
	with temporary_path.open('w', encoding='utf-8') as handle:
		progress = tqdm.tqdm(
			iter_parallel_fetches(
				details_url,
				chunks(compound_list, batch_size),
				workers,
				query_delay_seconds,
				golden_bioassay,
				indication,
			),
			total=total_batches,
			desc='AnnotationDB',
			unit='batch',
			dynamic_ncols=True,
		)
		for records, batch_missing, batch_failed in progress:
			missing_cids.extend(batch_missing)
			failed_cids.extend(batch_failed)
			for record in records:
				handle.write(json.dumps(record, ensure_ascii=True) + '\n')
				written += 1
			progress.set_postfix(
				written=written,
				failed=len(set(failed_cids)),
				missing=len(set(missing_cids)),
			)

	if missing_cids:
		emit_status(
			'Warning: '
			f'{len(set(missing_cids))} CIDs missing from AnnotationDB /compound/many response'
		)
	if failed_cids:
		emit_status(f'Warning: {len(set(failed_cids))} CIDs failed after retries')
	if missing_cids or failed_cids or written != len(compound_list):
		temporary_path.unlink(missing_ok=True)
		message = (
			'AnnotationDB fetch was incomplete: '
			f'index={len(compound_list)}, written={written}, '
			f'missing={len(set(missing_cids))}, failed={len(set(failed_cids))}'
		)
		raise RuntimeError(message)
	temporary_path.replace(outpath)
	manifest_path = manifest_path or f'{output_path}.manifest.tsv'
	manifest = Path(manifest_path)
	manifest.parent.mkdir(parents=True, exist_ok=True)
	manifest_temporary = manifest.with_suffix(manifest.suffix + '.partial')
	hash_value = hashlib.sha256()
	with outpath.open('rb') as handle:
		for block in iter(lambda: handle.read(1024 * 1024), b''):
			hash_value.update(block)
	manifest_temporary.write_text(
		'Source\tIndex.Endpoint\tDetails.Endpoint\tRecord.Count\t'
		'Indication.Requested\tFetched.At.UTC\tSHA256\n'
		f'AnnotationDB\t{db_url}\t{details_url}\t{written}\t'
		f'{str(indication).lower()}\t{datetime.now(timezone.utc).isoformat()}\t'
		f'{hash_value.hexdigest()}\n',
		encoding='utf-8',
	)
	manifest_temporary.replace(manifest)
	emit_status(
		f'Wrote {written} compound records to {outpath} '
		f'(batch_size={batch_size}, workers={workers}, '
		f'golden_bioassay={golden_bioassay})'
	)


def main_from_snakemake() -> None:
	main(
		db_url=snakemake.params.db_url,
		output_path=str(snakemake.output.raw),
		manifest_path=str(snakemake.output.manifest),
		details_url=snakemake.params.details_url,
		batch_size=int(snakemake.params.batch_size),
		workers=int(snakemake.threads),
		query_delay_seconds=float(snakemake.params.query_delay_seconds),
		golden_bioassay=bool(snakemake.params.golden_bioassay),
		indication=bool(snakemake.params.indication),
		api_key_file=snakemake.params.api_key_file,
		api_key_name=snakemake.params.api_key_name,
	)


if __name__ == '__main__':
	if 'snakemake' in globals():
		main_from_snakemake()
	else:
		parser = argparse.ArgumentParser(
			prog='fetch_annotationdb',
			description='Batch fetch AnnotationDB compound details into compact JSONL',
		)
		parser.add_argument('-u', required=True, help='/compound/all endpoint')
		parser.add_argument('-o', required=True, help='Output JSONL path')
		parser.add_argument('--manifest-path', default=None)
		parser.add_argument(
			'--details-url',
			default=None,
			help='Optional /compound/many endpoint override',
		)
		parser.add_argument(
			'--batch-size',
			type=int,
			default=DEFAULT_BATCH_SIZE,
			help=(
				f'Initial request batch size, must be <= {MAX_QUERY_SIZE} '
				f'(default: {DEFAULT_BATCH_SIZE})'
			),
		)
		parser.add_argument(
			'--workers',
			type=int,
			default=DEFAULT_WORKERS,
			help=f'Concurrent request workers (default: {DEFAULT_WORKERS})',
		)
		parser.add_argument(
			'--query-delay-seconds',
			type=float,
			default=DEFAULT_QUERY_DELAY_SECONDS,
			help='Delay between AnnotationDB queries per worker (default: 0.1)',
		)
		parser.add_argument(
			'--golden-bioassay',
			action=argparse.BooleanOptionalAction,
			default=True,
			help='Request only golden bioassays (default: true)',
		)
		parser.add_argument(
			'--limit',
			type=int,
			default=None,
			help='Only fetch the first N compounds, for sanity checks or benchmarking',
		)
		parser.add_argument(
			'--api-key-file',
			default=DEFAULT_ENV_FILE,
			help=f'Optional ignored API key file (default: {DEFAULT_ENV_FILE})',
		)
		parser.add_argument(
			'--api-key-name',
			default=DEFAULT_API_KEY_ENV_VAR,
			help=f'Key name in the API key file (default: {DEFAULT_API_KEY_ENV_VAR})',
		)
		parser.add_argument(
			'--indication',
			action=argparse.BooleanOptionalAction,
			default=True,
			help='Request drug indications (default: true)',
		)
		args = parser.parse_args()

		main(
			db_url=args.u,
			output_path=args.o,
			manifest_path=args.manifest_path,
			details_url=args.details_url,
			batch_size=args.batch_size,
			workers=args.workers,
			query_delay_seconds=args.query_delay_seconds,
			golden_bioassay=args.golden_bioassay,
			indication=args.indication,
			limit=args.limit,
			api_key_file=args.api_key_file,
			api_key_name=args.api_key_name,
		)
