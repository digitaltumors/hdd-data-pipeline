# Changelog

## 3.2.0

- Deprecate the DeepChem SIDER and Tox21 sources; remove the SIDER experiment.
- Re-curate Tox21 from the original, checksum-pinned NCATS 2014 Challenge training archive.
- Reconcile calls by HDD compound and assay; mask active/inactive conflicts as missing.
- Keep the existing Tox21 assay rows and export format, and retain DeepChem BBBP as compound metadata.

## 3.1.0

- Add AnnotationDB DailyMed boxed warnings, warnings/precautions and adverse reactions for 1,750 compounds.
- Preserve label identifiers, selection provenance, source dates, null/empty text and retrieval hashes.
- Preserve all 648,122 compounds, 47 existing metadata columns, seven assays and 55,506 indications from the checksum-pinned v3 release.
- Fetch labels alongside existing AnnotationDB annotations and retain a single QC report updated in place.

## 3.0.0

- Add 55,506 keyed ChEMBL compound indication records from the verified
  AnnotationDB release snapshot.
- Remove deprecated ClinTox and ToxCast experiments.
- Replace the non-regenerative PubChem DeepChem mapping dependency with local
  exact-SMILES/full-InChIKey matching, explicit AnnotationDB preference,
  ambiguity drops, collision handling, and match audits.
- Regenerate HDD from configured sources and export TSV/Matrix Market tables.
- Add release validation, source hashing, and v3 QC.
