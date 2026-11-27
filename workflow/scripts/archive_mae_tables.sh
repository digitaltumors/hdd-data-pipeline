#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 || ! ${3:-6} =~ ^[1-9]$ ]]; then
    printf 'Usage: %s TABLE_DIRECTORY ARCHIVE_PATH [GZIP_LEVEL]\n' "$0" >&2
    exit 2
fi

tables=$(cd "$1" && pwd -P)
archive=$2
level=${3:-6}
archive_dir=$(dirname "$archive")
mkdir -p "$archive_dir"
archive_dir=$(cd "$archive_dir" && pwd -P)
archive="$archive_dir/$(basename "$archive")"
if [[ "$archive_dir/" == "$tables/"* || -d "$archive" ]]; then
    printf 'Archive must be a file outside the table directory: %s\n' "$archive" >&2
    exit 2
fi

tmp=$(mktemp "$archive.tmp.XXXXXX")
trap 'rm -f -- "$tmp"' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

# Disable macOS copyfile metadata and ignore sidecars already on disk.
COPYFILE_DISABLE=1 tar --exclude='._*' --exclude='.snakemake_timestamp' \
    -C "$(dirname "$tables")" -cf - -- "$(basename "$tables")" | gzip -n "-$level" >"$tmp"
chmod 0644 "$tmp"
mv -f -- "$tmp" "$archive"
