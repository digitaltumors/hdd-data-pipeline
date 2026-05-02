from __future__ import annotations

import sys
import tarfile
from pathlib import Path

EXPECTED_ARG_COUNT = 4
MAX_GZIP_COMPRESSION = 9


def should_skip_export_file(path: Path) -> bool:
	if path.name.startswith('.'):
		return True
	return (
		path.parent.name == 'assays'
		and path.suffix == '.csv'
		and path.name.startswith('fingerprint.')
	)


def add_export_tree(archive: tarfile.TarFile, export_dir: Path) -> None:
	root_name = export_dir.name
	for path in sorted(export_dir.rglob('*')):
		if path.is_dir() or should_skip_export_file(path):
			continue
		archive.add(path, arcname=Path(root_name) / path.relative_to(export_dir))


def add_fingerprint_files(
	archive: tarfile.TarFile,
	export_dir: Path,
	fingerprint_dir: Path,
) -> None:
	archive_root = Path(export_dir.name) / 'assays'
	for path in sorted(
		[*fingerprint_dir.glob('*.mtx'), *fingerprint_dir.glob('*.tsv')]
	):
		archive.add(path, arcname=archive_root / path.name)


def main() -> int:
	if len(sys.argv) != EXPECTED_ARG_COUNT:
		sys.stderr.write(
			'Usage: archive_mae_csvs.py <export_dir> <fingerprint_dir> <archive_path>\n'
		)
		return 2

	export_dir = Path(sys.argv[1])
	fingerprint_dir = Path(sys.argv[2])
	archive_path = Path(sys.argv[3])

	if not export_dir.is_dir():
		message = f'Export directory does not exist: {export_dir}'
		raise FileNotFoundError(message)
	if not fingerprint_dir.is_dir():
		message = f'Fingerprint directory does not exist: {fingerprint_dir}'
		raise FileNotFoundError(message)

	archive_path.parent.mkdir(parents=True, exist_ok=True)
	with tarfile.open(
		archive_path,
		mode='w:gz',
		compresslevel=MAX_GZIP_COMPRESSION,
	) as archive:
		add_export_tree(archive, export_dir)
		add_fingerprint_files(archive, export_dir, fingerprint_dir)

	return 0


if __name__ == '__main__':
	raise SystemExit(main())
