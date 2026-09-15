# Changelog

## 3.0.0

- Add 55,506 keyed ChEMBL compound indication records from the verified
  AnnotationDB release snapshot.
- Remove deprecated ClinTox and ToxCast experiments.
- Replace the non-regenerative PubChem DeepChem mapping dependency with local
  exact-SMILES/full-InChIKey matching, explicit AnnotationDB preference,
  ambiguity drops, collision handling, and match audits.
- Regenerate HDD from configured sources and export TSV/Matrix Market tables.
- Add release validation, source hashing, and v3 QC.
