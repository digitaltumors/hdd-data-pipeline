"""Build the HDD Tox21 assay matrix from the NCATS Challenge training archive."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import logging
import re
from collections import Counter, defaultdict
from pathlib import Path
from zipfile import ZipFile

from rdkit import Chem

MEMBER = 'tox21_10k_data_all.sdf'
TRAINING_SHA256 = '024a3ae2690bcd4a593e6e0b10b455470b9bcb1d8f299dd36f220a250181517b'
ASSAYS = (
    'NR-AR', 'NR-AR-LBD', 'NR-AhR', 'NR-Aromatase', 'NR-ER', 'NR-ER-LBD',
    'NR-PPAR-gamma', 'SR-ARE', 'SR-ATAD5', 'SR-HSE', 'SR-MMP', 'SR-p53',
)
TAG = re.compile(r'^>\s+<([^>]+)>')
INCHIKEY = re.compile(r'[A-Z]{14}-[A-Z]{10}-[A-Z]')
LOG = logging.getLogger(__name__)


def sha256_file(path: Path) -> str:
    with path.open('rb') as raw:
        return hashlib.file_digest(raw, 'sha256').hexdigest()


def read_training(path: Path) -> list[dict[str, str]]:
    """Keep sample IDs and source label calls; do not infer a label from absence."""
    sha256 = sha256_file(path)
    if sha256 != TRAINING_SHA256:
        msg = f'Training archive SHA-256 mismatch: {sha256}'
        raise ValueError(msg)
    with ZipFile(path) as archive, archive.open(MEMBER) as member:
        with io.TextIOWrapper(member, encoding='utf-8') as source:
            return parse_sdf(source)


def parse_record(lines: list[str]) -> dict[str, str]:
    """Read the source sample title and explicit challenge assay calls."""
    sample_id = lines[0].strip() if lines else ''
    if not sample_id or 'M  END' not in lines:
        msg = 'SDF record is missing its title or molblock terminator'
        raise ValueError(msg)
    props = {'source_sample_id': sample_id}
    for i in range(lines.index('M  END') + 1, len(lines)):
        tag = TAG.match(lines[i])
        if not tag:
            continue
        name = tag[1]
        if name not in ASSAYS:
            continue
        value = lines[i + 1].strip() if i + 1 < len(lines) else ''
        if name in props or value not in {'0', '1'}:
            msg = f'Unexpected {name} property for {sample_id}'
            raise ValueError(msg)
        props[name] = value
    return props


def parse_sdf(lines: io.TextIOBase) -> list[dict[str, str]]:
    """Parse records without dependencies or normalizing legacy identifiers."""
    records: list[dict[str, str]] = []
    record: list[str] = []
    for line in lines:
        line = line.rstrip('\r\n')
        if line.strip() == '$$$$':
            records.append(parse_record(record))
            record = []
        else:
            record.append(line)
    if any(line.strip() for line in record) or not records:
        msg = 'No complete SDF or unterminated final record'
        raise ValueError(msg)
    return records


def hdd_key_index(path: Path) -> dict[str, set[str]]:
    """Index full InChIKeys, keeping all candidate HDD IDs per key."""
    by_key: dict[str, set[str]] = defaultdict(set)
    with path.open(encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not {'HDD.Compound.ID', 'InChIKey'} <= set(reader.fieldnames or []):
            msg = 'HDD metadata must contain HDD.Compound.ID and InChIKey'
            raise ValueError(msg)
        for row in reader:
            compound_id = row['HDD.Compound.ID']
            key = row['InChIKey']
            if compound_id and key and INCHIKEY.fullmatch(key):
                by_key[key].add(compound_id)
    return by_key


def source_keys(path: Path, records: list[dict[str, str]]) -> list[str]:
    """Calculate keys without changing the NCATS source structures."""
    result = []
    with ZipFile(path) as archive, archive.open(MEMBER) as member:
        for index, mol in enumerate(Chem.ForwardSDMolSupplier(member)):
            if mol is not None and (
                index >= len(records)
                or mol.GetProp('_Name').strip() != records[index]['source_sample_id']
            ):
                msg = 'RDKit structure and source sample IDs do not agree'
                raise ValueError(msg)
            result.append((Chem.MolToInchiKey(mol) or '') if mol is not None else '')
    if len(result) != len(records):
        msg = 'SDF structure and source-record counts do not agree'
        raise ValueError(msg)
    return result


def collapse_calls(
    records: list[dict[str, str]], keys: list[str], by_key: dict[str, set[str]],
) -> tuple[list[str], list[list[str]], Counter[str], int]:
    """Keep unanimous calls; encode discordant or absent calls as NA."""
    grouped: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    statuses: Counter[str] = Counter()
    for record, key in zip(records, keys, strict=True):
        ids = by_key.get(key, set())
        if not key:
            statuses['no_structure_key'] += 1
        elif not ids:
            statuses['not_in_hdd'] += 1
        elif len(ids) != 1:
            statuses['ambiguous_hdd_key'] += 1
        else:
            statuses['unique_hdd_key'] += 1
            calls = grouped[next(iter(ids))]
            for assay in ASSAYS:
                if assay in record:
                    calls[assay].add(record[assay])

    compounds = sorted(grouped)
    matrix = []
    discordant = 0
    for assay in ASSAYS:
        row = [assay]
        for compound in compounds:
            calls = grouped[compound].get(assay, set())
            row.append(next(iter(calls)) if len(calls) == 1 else 'NA')
            discordant += len(calls) > 1
        matrix.append(row)
    return compounds, matrix, statuses, discordant


def build(archive: Path, coldata: Path, assay: Path) -> None:
    records = read_training(archive)
    if len({r['source_sample_id'] for r in records}) != len(records):
        msg = 'NCATS training sample IDs are not unique'
        raise ValueError(msg)
    compounds, matrix, statuses, discordant = collapse_calls(
        records, source_keys(archive, records), hdd_key_index(coldata),
    )
    if not compounds or not all(any(call != 'NA' for call in row[1:]) for row in matrix):
        msg = 'No HDD matches or one or more assays have no retained calls'
        raise ValueError(msg)

    LOG.info('Mapped %d HDD compounds from %d NCATS samples.', len(compounds), len(records))
    LOG.info('NCATS sample mapping status counts: %s', dict(sorted(statuses.items())))
    LOG.info('Masked %d discordant compound–assay pairs.', discordant)
    assay.parent.mkdir(parents=True, exist_ok=True)
    with assay.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(['Tox Assay', *compounds])
        writer.writerows(matrix)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('archive', 'coldata', 'assay'):
        parser.add_argument(f'--{name}', type=Path, required=True)
    args = parser.parse_args()
    build(args.archive, args.coldata, args.assay)


if __name__ == '__main__':
    main()
