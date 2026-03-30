#!/usr/bin/env bash

set -euo pipefail

project_root="$(pwd)"

target_dirs=(
  "${project_root}/data/rawdata"
  "${project_root}/data/procdata"
  "${project_root}/data/results"
)

printf '%s\n' "This will delete generated pipeline artifacts under:"
printf '  %s\n' "${target_dirs[@]}"
printf '%s\n' "README.md files and the parent directories will be preserved."
read -r -p "Are you sure? [y/N] " response

case "${response}" in
  [yY]|[yY][eE][sS])
    ;;
  *)
    printf '%s\n' "Cancelled."
    exit 0
    ;;
esac

for dir in "${target_dirs[@]}"; do
  if [[ ! -d "${dir}" ]]; then
    mkdir -p "${dir}"
    continue
  fi

  find "${dir}" -mindepth 1 -type f ! -name 'README.md' -delete
  find "${dir}" -mindepth 1 -type l -delete
  find "${dir}" -depth -mindepth 1 -type d -empty -delete
done

printf '%s\n' "Pipeline artifacts removed."
