from __future__ import annotations

import sys
import tarfile
from pathlib import Path

MIN_ARG_COUNT = 3
MAX_ARG_COUNT = 4
MAX_GZIP_COMPRESSION = 9


def add_export_tree(archive: tarfile.TarFile, export_dir: Path) -> None:
	root_name = export_dir.name
	for path in sorted(export_dir.rglob('*')):
		if path.is_dir() or path.name.startswith('.'):
			continue
		archive.add(path, arcname=Path(root_name) / path.relative_to(export_dir))


def main() -> int:
	if len(sys.argv) not in {MIN_ARG_COUNT, MAX_ARG_COUNT}:
		sys.stderr.write(
			'Usage: archive_mae_csvs.py <export_dir> [legacy_fingerprint_dir] '
			'<archive_path>\n'
		)
		return 2

	export_dir = Path(sys.argv[1])
	archive_path = Path(sys.argv[-1])

	if not export_dir.is_dir():
		message = f'Export directory does not exist: {export_dir}'
		raise FileNotFoundError(message)

	archive_path.parent.mkdir(parents=True, exist_ok=True)
	with tarfile.open(
		archive_path,
		mode='w:gz',
		compresslevel=MAX_GZIP_COMPRESSION,
	) as archive:
		add_export_tree(archive, export_dir)

	return 0


if __name__ == '__main__':
	raise SystemExit(main())
