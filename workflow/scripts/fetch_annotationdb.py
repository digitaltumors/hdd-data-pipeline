import argparse
import json
import threading
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Sequence, Tuple
from urllib.parse import urlparse

import requests
import tqdm

DEFAULT_BATCH_SIZE = 250
DEFAULT_WORKERS = 6
DEFAULT_QUERY_DELAY_SECONDS = 0.1
MAX_QUERY_SIZE = 250
RETRIES = 3
TIMEOUT = 60
SPLIT_STATUS_CODES = {413, 422, 429, 500, 502, 503, 504}
SESSION_HEADERS = {
	'accept': 'application/json',
	'user-agent': 'hdd-data-pipeline/1.0',
}

thread_local = threading.local()


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
	return params


def fetch_split_batch(
	session: requests.Session,
	details_url: str,
	cids: List[str],
	failed_cids: List[str],
	query_delay_seconds: float,
	golden_bioassay: bool,
) -> List[dict]:
	mid = len(cids) // 2
	left = fetch_details_for_cids(
		session,
		details_url,
		cids[:mid],
		failed_cids,
		query_delay_seconds,
		golden_bioassay,
	)
	right = fetch_details_for_cids(
		session,
		details_url,
		cids[mid:],
		failed_cids,
		query_delay_seconds,
		golden_bioassay,
	)
	return left + right


def fetch_chunked_details(
	session: requests.Session,
	details_url: str,
	cids: List[str],
	failed_cids: List[str],
	query_delay_seconds: float,
	golden_bioassay: bool,
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
			)
		)
	return results


def record_failed_cids(failed_cids: List[str], cids: List[str]) -> List[dict]:
	failed_cids.extend(cids)
	return []


def resolve_failed_fetch(
	session: requests.Session,
	details_url: str,
	cids: List[str],
	failed_cids: List[str],
	query_delay_seconds: float,
	golden_bioassay: bool,
	attempt: int,
	allow_split: bool,
) -> Optional[List[dict]]:
	if allow_split and len(cids) > 1:
		return fetch_split_batch(
			session,
			details_url,
			cids,
			failed_cids,
			query_delay_seconds,
			golden_bioassay,
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


def fetch_details_for_cids(
	session: requests.Session,
	details_url: str,
	cids: List[str],
	failed_cids: List[str],
	query_delay_seconds: float,
	golden_bioassay: bool,
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
		)

	params = build_details_params(cids, golden_bioassay)

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


def slim_detail_record(detail: dict) -> dict:
	toxicity = detail.get('toxicity') or {}
	mechanisms = detail.get('mechanisms') or []
	bioassays = detail.get('bioassays') or []
	return {
		'cid': detail.get('cid'),
		'molecular_formula': detail.get('molecular_formula'),
		'iupac_name': detail.get('iupac_name'),
		'molecule_chembl_id': detail.get('molecule_chembl_id'),
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
		'toxicity': {
			'dili_severity_grade': toxicity.get('dili_severity_grade'),
			'dili_annotation': toxicity.get('dili_annotation'),
			'hepatotoxicity_likelihood_score': toxicity.get(
				'hepatotoxicity_likelihood_score'
			),
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


def iter_parallel_fetches(
	details_url: str,
	batches: Iterable[List[dict]],
	workers: int,
	query_delay_seconds: float,
	golden_bioassay: bool,
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
					)
				)

			if not pending:
				break

			done, pending = wait(pending, return_when=FIRST_COMPLETED)
			for future in done:
				yield future.result()


def main(
	db_url: str,
	output_path: str,
	details_url: Optional[str] = None,
	batch_size: int = DEFAULT_BATCH_SIZE,
	workers: int = DEFAULT_WORKERS,
	query_delay_seconds: float = DEFAULT_QUERY_DELAY_SECONDS,
	golden_bioassay: bool = True,
	limit: Optional[int] = None,
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
	with outpath.open('w', encoding='utf-8') as handle:
		progress = tqdm.tqdm(
			iter_parallel_fetches(
				details_url,
				chunks(compound_list, batch_size),
				workers,
				query_delay_seconds,
				golden_bioassay,
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
	emit_status(
		f'Wrote {written} compound records to {outpath} '
		f'(batch_size={batch_size}, workers={workers}, '
		f'golden_bioassay={golden_bioassay})'
	)


def main_from_snakemake() -> None:
	main(
		db_url=snakemake.params.db_url,
		output_path=str(snakemake.output.raw),
		details_url=snakemake.params.details_url,
		batch_size=int(snakemake.params.batch_size),
		workers=int(snakemake.threads),
		query_delay_seconds=float(snakemake.params.query_delay_seconds),
		golden_bioassay=bool(snakemake.params.golden_bioassay),
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
			help='Only fetch the first N compounds, for smoke tests or benchmarking',
		)
		args = parser.parse_args()

		main(
			db_url=args.u,
			output_path=args.o,
			details_url=args.details_url,
			batch_size=args.batch_size,
			workers=args.workers,
			query_delay_seconds=args.query_delay_seconds,
			golden_bioassay=args.golden_bioassay,
			limit=args.limit,
		)
