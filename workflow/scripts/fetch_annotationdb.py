# ruff: noqa: ANN201, EM101, PLR0911, T201, TRY003, TRY004, TRY300, TRY301

import argparse
import json
import threading
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple
from urllib.parse import urlparse

import requests
import tqdm

DEFAULT_BATCH_SIZE = 5
DEFAULT_WORKERS = 8
RETRIES = 3
TIMEOUT = 60
SPLIT_STATUS_CODES = {413, 429, 500, 502, 503, 504}

thread_local = threading.local()


def fetch_json(
	session: requests.Session,
	url: str,
	params: Optional[Dict[str, str]] = None,
	retries: int = RETRIES,
	timeout: int = TIMEOUT,
):
	for attempt in range(retries):
		try:
			resp = session.get(url, params=params, timeout=timeout)
			resp.raise_for_status()
			return resp.json()
		except Exception:
			if attempt == retries - 1:
				raise
			time.sleep(2 * (attempt + 1))


def fetch_details_for_cids(
	session: requests.Session,
	details_url: str,
	cids: List[str],
	failed_cids: List[str],
) -> List[Dict[str, object]]:
	if not cids:
		return []

	params = {
		'compounds': ','.join(cids),
		'format': 'json',
		'bioassay': 'true',
		'mechanism': 'true',
		'toxicity': 'true',
	}

	for attempt in range(RETRIES):
		try:
			resp = session.get(details_url, params=params, timeout=TIMEOUT)
			if resp.status_code in SPLIT_STATUS_CODES:
				raise requests.HTTPError(response=resp)
			resp.raise_for_status()
			data = resp.json()
			if not isinstance(data, list):
				raise ValueError('Expected list response from /compound/many')
			return data
		except requests.HTTPError as exc:
			status = exc.response.status_code if exc.response is not None else None
			if status in SPLIT_STATUS_CODES and len(cids) > 1:
				mid = len(cids) // 2
				left = fetch_details_for_cids(
					session, details_url, cids[:mid], failed_cids
				)
				right = fetch_details_for_cids(
					session, details_url, cids[mid:], failed_cids
				)
				return left + right
			if attempt == RETRIES - 1:
				failed_cids.extend(cids)
				return []
			time.sleep(2 * (attempt + 1))
		except (requests.RequestException, ValueError):
			if len(cids) > 1:
				mid = len(cids) // 2
				left = fetch_details_for_cids(
					session, details_url, cids[:mid], failed_cids
				)
				right = fetch_details_for_cids(
					session, details_url, cids[mid:], failed_cids
				)
				return left + right
			if attempt == RETRIES - 1:
				failed_cids.extend(cids)
				return []
			time.sleep(2 * (attempt + 1))

	failed_cids.extend(cids)
	return []


def derive_details_url(db_url: str) -> str:
	parsed = urlparse(db_url)
	path = parsed.path
	if path.endswith('/compound/all'):
		path = path.replace('/compound/all', '/compound/many')
	else:
		path = path.rstrip('/') + '/many'
	return parsed._replace(path=path, query='').geturl()


def chunks(
	items: List[Dict[str, object]], size: int
) -> Iterator[List[Dict[str, object]]]:
	for idx in range(0, len(items), size):
		yield items[idx : idx + size]


def normalize_cid(value: object) -> Optional[str]:
	if value is None:
		return None
	return str(value)


def slim_detail_record(detail: Dict[str, object]) -> Dict[str, object]:
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
				'activity_outcome_method': assay.get('activity_outcome_method'),
			}
			for assay in bioassays
			if isinstance(assay, dict) and assay.get('aid') is not None
		],
	}


def get_thread_session() -> requests.Session:
	session = getattr(thread_local, 'session', None)
	if session is None:
		session = requests.Session()
		thread_local.session = session
	return session


def fetch_batch_records(
	details_url: str, batch: List[Dict[str, object]]
) -> Tuple[List[Dict[str, object]], List[str], List[str]]:
	session = get_thread_session()
	failed_cids: List[str] = []
	missing_cids: List[str] = []

	cids = [
		normalize_cid(item.get('cid'))
		for item in batch
		if normalize_cid(item.get('cid')) is not None
	]
	details = fetch_details_for_cids(session, details_url, cids, failed_cids)

	details_by_cid = {
		normalize_cid(item.get('cid')): slim_detail_record(item)
		for item in details
		if normalize_cid(item.get('cid')) is not None
	}

	records: List[Dict[str, object]] = []
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
	batches: Iterable[List[Dict[str, object]]],
	workers: int,
) -> Iterator[Tuple[List[Dict[str, object]], List[str], List[str]]]:
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
				pending.add(executor.submit(fetch_batch_records, details_url, batch))

			if not pending:
				break

			done, pending = wait(pending, return_when=FIRST_COMPLETED)
			for future in done:
				yield future.result()


def main(
	db_url: str,
	output_path: str,
	batch_size: int = DEFAULT_BATCH_SIZE,
	workers: int = DEFAULT_WORKERS,
	limit: Optional[int] = None,
) -> None:
	if batch_size < 1:
		raise ValueError('batch_size must be >= 1')
	if workers < 1:
		raise ValueError('workers must be >= 1')
	if limit is not None and limit < 1:
		raise ValueError('limit must be >= 1')

	outpath = Path(output_path)
	outpath.parent.mkdir(parents=True, exist_ok=True)

	session = requests.Session()
	compound_list = fetch_json(session, db_url)
	if not isinstance(compound_list, list):
		raise ValueError('Expected list response from /compound/all')
	if limit is not None:
		compound_list = compound_list[:limit]

	details_url = derive_details_url(db_url)
	missing_cids: List[str] = []
	failed_cids: List[str] = []
	written = 0

	total_batches = (len(compound_list) + batch_size - 1) // batch_size
	with outpath.open('w', encoding='utf-8') as handle:
		progress = tqdm.tqdm(
			iter_parallel_fetches(
				details_url, chunks(compound_list, batch_size), workers
			),
			total=total_batches,
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
		print(
			'Warning: '
			f'{len(set(missing_cids))} CIDs missing from /compound/many response',
			flush=True,
		)
	if failed_cids:
		print(
			f'Warning: {len(set(failed_cids))} CIDs failed after retries',
			flush=True,
		)
	print(
		f'Wrote {written} compound records to {outpath} '
		f'(batch_size={batch_size}, workers={workers})',
		flush=True,
	)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		prog='fetch_annotationdb',
		description='Batch fetch AnnotationDB compound details into compact JSONL',
	)
	parser.add_argument('-u', required=True, help='/compound/all endpoint')
	parser.add_argument('-o', required=True, help='Output JSONL path')
	parser.add_argument(
		'--batch-size',
		type=int,
		default=DEFAULT_BATCH_SIZE,
		help=f'Initial request batch size (default: {DEFAULT_BATCH_SIZE})',
	)
	parser.add_argument(
		'--workers',
		type=int,
		default=DEFAULT_WORKERS,
		help=f'Concurrent request workers (default: {DEFAULT_WORKERS})',
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
		batch_size=args.batch_size,
		workers=args.workers,
		limit=args.limit,
	)
