import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi

COLData_COLUMNS = [
	'HDD.Compound.ID',
	'Pubchem.CID',
	'InChIKey',
	'SMILES',
	'In.AnnotationDB',
	'AnnotationDB.Name',
]
UNMATCHED_DECISIONS = {
	'dropped_unmatched',
	'dropped_parse_failure',
	'dropped_missing_smiles',
}

RDLogger.DisableLog('rdApp.*')


def normalize_boolean(value: object) -> bool:
	if pd.isna(value):
		return False
	if isinstance(value, bool):
		return value
	return str(value).strip().lower() in {'true', 't', '1', 'yes'}


def structure_identifiers(smiles: object) -> tuple[object, object]:
	if pd.isna(smiles) or not str(smiles).strip():
		return pd.NA, pd.NA
	molecule = Chem.MolFromSmiles(str(smiles))
	if molecule is None:
		return pd.NA, pd.NA
	canonical_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
	return canonical_smiles, inchi.MolToInchiKey(molecule)


def candidate_indexes(
	coldata: pd.DataFrame,
) -> tuple[dict[str, list[int]], dict[str, list[int]]]:
	smiles_index: dict[str, list[int]] = defaultdict(list)
	inchikey_index: dict[str, list[int]] = defaultdict(list)
	for index, row in coldata.iterrows():
		if pd.notna(row['SMILES']):
			smiles_index[str(row['SMILES'])].append(index)
		if pd.notna(row['InChIKey']):
			inchikey_index[str(row['InChIKey'])].append(index)
	return dict(smiles_index), dict(inchikey_index)


def measurement_conflicts(
	source: pd.DataFrame,
	source_indexes: list[int],
	measurement_columns: list[str],
) -> list[str]:
	conflicts = []
	values = source.loc[source_indexes, measurement_columns]
	for column in measurement_columns:
		numeric = pd.to_numeric(values[column], errors='coerce').dropna()
		if not numeric.empty and not numeric.eq(numeric.iloc[0]).all():
			conflicts.append(column)
	return conflicts


def build_assay_matrix(
	source: pd.DataFrame,
	source_audit: pd.DataFrame,
	measurement_columns: list[str],
	index_name: str,
) -> pd.DataFrame:
	retained = source_audit.dropna(subset=['Final.HDD.Compound.ID']).copy()
	retained = retained.sort_values('Source.Row.Index')
	source_indexes = retained['Source.Row.Index'].astype(int).tolist()
	assay_source = source.loc[source_indexes, measurement_columns].copy()
	assay_source.insert(
		0,
		'HDD.Compound.ID',
		retained['Final.HDD.Compound.ID'].astype(str).tolist(),
	)
	if assay_source['HDD.Compound.ID'].duplicated().any():
		message = 'DeepChem assay contains duplicate HDD.Compound.ID values'
		raise ValueError(message)
	assay = assay_source.transpose()
	assay = assay.rename(columns=assay.iloc[0])
	assay = assay.iloc[1:,]
	return assay.reset_index(names=index_name, drop=False)


def choose_candidates(  # noqa: PLR0913, PLR0917
	dataset: str,
	source_index: int,
	source_id: object,
	source_smiles: object,
	coldata: pd.DataFrame,
	smiles_index: dict[str, list[int]],
	inchikey_index: dict[str, list[int]],
) -> tuple[dict[str, object], list[dict[str, object]]]:
	canonical_smiles, source_inchikey = structure_identifiers(source_smiles)
	route = 'none'
	candidate_indices: list[int] = []

	if pd.isna(source_smiles) or not str(source_smiles).strip():
		match_case = 'missing_smiles'
		final_decision = 'dropped_missing_smiles'
		reason = 'Source SMILES is missing'
	else:
		candidate_indices = smiles_index.get(str(source_smiles), [])
		if candidate_indices:
			route = 'exact_smiles'
		elif pd.isna(source_inchikey):
			match_case = 'rdkit_parse_failure'
			final_decision = 'dropped_parse_failure'
			reason = 'RDKit could not parse source SMILES'
		else:
			route = 'full_inchikey'
			candidate_indices = inchikey_index.get(str(source_inchikey), [])

	candidates = coldata.loc[candidate_indices].copy()
	candidate_count = len(candidates)
	annotationdb_mask = candidates['In.AnnotationDB'].map(normalize_boolean)
	annotationdb_count = int(annotationdb_mask.sum())
	provisional_id: object = pd.NA

	if route != 'none':
		route_label = 'exact' if route == 'exact_smiles' else 'inchikey'
		if candidate_count == 0:
			match_case = 'inchikey_unmatched'
			final_decision = 'dropped_unmatched'
			reason = 'Full InChIKey has no HDD candidate'
		elif candidate_count == 1:
			match_case = f'{route_label}_unique'
			final_decision = 'accepted'
			reason = 'Exactly one HDD candidate'
			provisional_id = candidates.iloc[0]['HDD.Compound.ID']
		elif annotationdb_count == 1:
			match_case = f'{route_label}_annotationdb_preference'
			final_decision = 'accepted'
			reason = 'Exactly one candidate is AnnotationDB-backed'
			provisional_id = candidates.loc[annotationdb_mask, 'HDD.Compound.ID'].iloc[
				0
			]
		else:
			match_case = f'{route_label}_ambiguous'
			final_decision = 'dropped_ambiguous'
			reason = 'Multiple candidates remain after AnnotationDB preference'

	source_row_number = source_index + 1
	source_audit = {
		'Dataset': dataset,
		'Source.Row.Index': source_index,
		'Source.Row.Number': source_row_number,
		'Source.ID': source_id,
		'Source.SMILES': source_smiles,
		'RDKit.Canonical.SMILES': canonical_smiles,
		'RDKit.InChIKey': source_inchikey,
		'Match.Route': route,
		'Match.Case': match_case,
		'Candidate.Count': candidate_count,
		'AnnotationDB.Candidate.Count': annotationdb_count,
		'Provisional.HDD.Compound.ID': provisional_id,
		'Final.HDD.Compound.ID': provisional_id,
		'Final.Decision': final_decision,
		'Decision.Reason': reason,
	}
	candidate_audit = []
	for _, candidate in candidates.iterrows():
		candidate_id = candidate['HDD.Compound.ID']
		candidate_audit.append(
			{
				'Dataset': dataset,
				'Source.Row.Number': source_row_number,
				'Source.ID': source_id,
				'Source.SMILES': source_smiles,
				'Match.Route': route,
				'Match.Case': match_case,
				'Candidate.HDD.Compound.ID': candidate_id,
				'Candidate.Pubchem.CID': candidate['Pubchem.CID'],
				'Candidate.InChIKey': candidate['InChIKey'],
				'Candidate.SMILES': candidate['SMILES'],
				'Candidate.In.AnnotationDB': normalize_boolean(
					candidate['In.AnnotationDB']
				),
				'Candidate.AnnotationDB.Name': candidate['AnnotationDB.Name'],
				'Selected.By.Match.Policy': bool(
					pd.notna(provisional_id) and candidate_id == provisional_id
				),
			}
		)
	return source_audit, candidate_audit


def resolve_many_to_one(  # noqa: PLR0913
	dataset: str,
	source: pd.DataFrame,
	source_audit: pd.DataFrame,
	candidate_audit: pd.DataFrame,
	measurement_columns: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
	accepted = source_audit.dropna(subset=['Provisional.HDD.Compound.ID'])
	counts = accepted.groupby('Provisional.HDD.Compound.ID')[
		'Source.Row.Index'
	].nunique()
	collision_ids = counts[counts > 1].index
	collision_rows = []

	for hdd_id in collision_ids:
		group = accepted.loc[
			accepted['Provisional.HDD.Compound.ID'] == hdd_id
		].sort_values('Source.Row.Index')
		source_indexes = group['Source.Row.Index'].astype(int).tolist()
		kept = group.iloc[0]
		discarded = group.iloc[1:]
		conflicts = measurement_conflicts(source, source_indexes, measurement_columns)

		for _, discarded_row in discarded.iterrows():
			audit_index = source_audit.index[
				source_audit['Source.Row.Index'] == discarded_row['Source.Row.Index']
			][0]
			source_audit.loc[audit_index, 'Final.HDD.Compound.ID'] = pd.NA
			source_audit.loc[audit_index, 'Final.Decision'] = (
				'dropped_many_to_one_keep_first'
			)
			source_audit.loc[audit_index, 'Decision.Reason'] = (
				'Another source row mapped to the same HDD ID and appeared earlier'
			)

		collision_rows.append(
			{
				'Dataset': dataset,
				'HDD.Compound.ID': hdd_id,
				'Source.Row.Count': len(group),
				'Kept.Source.Row.Number': int(kept['Source.Row.Number']),
				'Kept.Source.ID': kept['Source.ID'],
				'Discarded.Source.Row.Numbers': ','.join(
					str(int(value)) for value in discarded['Source.Row.Number']
				),
				'Discarded.Source.IDs': ','.join(
					str(value) for value in discarded['Source.ID']
				),
				'Conflicting.Measurement.Count': len(conflicts),
				'Conflicting.Measurements': ','.join(conflicts),
			}
		)

	retained_pairs = set(
		source_audit.dropna(subset=['Final.HDD.Compound.ID'])[
			['Source.Row.Number', 'Final.HDD.Compound.ID']
		].itertuples(index=False, name=None)
	)
	candidate_audit['Retained.After.Collision.Policy'] = candidate_audit.apply(
		lambda row: (
			(row['Source.Row.Number'], row['Candidate.HDD.Compound.ID'])
			in retained_pairs
		),
		axis=1,
	)
	return source_audit, candidate_audit, pd.DataFrame(collision_rows)


def make_summary(
	dataset: str,
	source_audit: pd.DataFrame,
	collision_audit: pd.DataFrame,
) -> pd.DataFrame:
	match_counts = source_audit['Match.Case'].value_counts()
	decision_counts = source_audit['Final.Decision'].value_counts()
	rows = [
		{'Dataset': dataset, 'Metric': 'source_rows', 'Value': len(source_audit)},
		{
			'Dataset': dataset,
			'Metric': 'final_compounds',
			'Value': source_audit['Final.HDD.Compound.ID'].nunique(),
		},
		{
			'Dataset': dataset,
			'Metric': 'collision_groups',
			'Value': len(collision_audit),
		},
	]
	rows.extend(
		{
			'Dataset': dataset,
			'Metric': f'match_case.{match_case}',
			'Value': count,
		}
		for match_case, count in match_counts.items()
	)
	rows.extend(
		{
			'Dataset': dataset,
			'Metric': f'final_decision.{decision}',
			'Value': count,
		}
		for decision, count in decision_counts.items()
	)
	return pd.DataFrame(rows)


def validate_expected_counts(  # noqa: PLR0917
	dataset_key: str,
	source: pd.DataFrame,
	assay: pd.DataFrame,
	source_audit: pd.DataFrame,
	collision_audit: pd.DataFrame,
	expected: dict[str, int] | None,
) -> None:
	if not expected:
		return
	observed = {
		'source_compounds': len(source),
		'features': len(assay),
		'compounds': len(assay.columns) - 1,
		'dropped_ambiguous': int(
			(source_audit['Final.Decision'] == 'dropped_ambiguous').sum()
		),
		'dropped_unmatched': int(
			source_audit['Final.Decision'].isin(UNMATCHED_DECISIONS).sum()
		),
		'dropped_many_to_one': int(
			(source_audit['Final.Decision'] == 'dropped_many_to_one_keep_first').sum()
		),
		'collision_groups': len(collision_audit),
	}
	mismatches = {
		key: (int(expected[key]), observed[key])
		for key in expected
		if key in observed and int(expected[key]) != observed[key]
	}
	if mismatches:
		details = ', '.join(
			f'{key}: expected {values[0]}, observed {values[1]}'
			for key, values in mismatches.items()
		)
		message = f'{dataset_key} DeepChem count mismatch: {details}'
		raise ValueError(message)


def process_dataset(  # noqa: PLR0913, PLR0917
	dataset: str,
	dataset_key: str,
	source: pd.DataFrame,
	coldata: pd.DataFrame,
	index_name: str,
	convert_to_int: bool,
	expected: dict[str, int] | None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
	measurement_columns = [
		column for column in source.columns if column not in {'smiles', 'mol_id'}
	]
	source = source.copy()
	source['smiles'] = source['smiles'].astype('string')
	if convert_to_int:
		source[measurement_columns] = (
			source[measurement_columns]
			.apply(lambda column: pd.to_numeric(column, errors='coerce'))
			.astype('Int64')
		)

	smiles_index, inchikey_index = candidate_indexes(coldata)
	source_rows = []
	candidate_rows = []
	for source_index, source_row in source.iterrows():
		source_id = source_row['mol_id'] if 'mol_id' in source.columns else pd.NA
		source_audit, candidate_audit = choose_candidates(
			dataset,
			source_index,
			source_id,
			source_row['smiles'],
			coldata,
			smiles_index,
			inchikey_index,
		)
		source_rows.append(source_audit)
		candidate_rows.extend(candidate_audit)

	source_audit = pd.DataFrame(source_rows)
	candidate_audit = pd.DataFrame(candidate_rows)
	source_audit, candidate_audit, collision_audit = resolve_many_to_one(
		dataset,
		source,
		source_audit,
		candidate_audit,
		measurement_columns,
	)
	assay = build_assay_matrix(source, source_audit, measurement_columns, index_name)
	summary = make_summary(dataset, source_audit, collision_audit)
	validate_expected_counts(
		dataset_key,
		source,
		assay,
		source_audit,
		collision_audit,
		expected,
	)
	return assay, source_audit, candidate_audit, collision_audit, summary


def main(  # noqa: PLR0913, PLR0917
	coldata_path: str,
	tox21_input: str,
	sider_input: str,
	tox21_output: str,
	sider_output: str,
	match_audit_output: str,
	candidate_audit_output: str,
	collision_audit_output: str,
	summary_output: str,
	expected_assays: dict[str, dict[str, int]] | None = None,
) -> None:
	coldata = pd.read_csv(
		coldata_path,
		sep='\t',
		usecols=COLData_COLUMNS,
		low_memory=False,
	)
	tox21_source = pd.read_csv(tox21_input)
	sider_source = pd.read_csv(sider_input)
	expected_assays = expected_assays or {}

	tox21 = process_dataset(
		'Tox21',
		'tox21',
		tox21_source,
		coldata,
		'Tox Assay',
		True,
		expected_assays.get('tox21'),
	)
	sider = process_dataset(
		'SIDER',
		'sider',
		sider_source,
		coldata,
		'Side Effect',
		False,
		expected_assays.get('sider'),
	)

	tox21_assay, *tox21_audits = tox21
	sider_assay, *sider_audits = sider
	match_audit = pd.concat([tox21_audits[0], sider_audits[0]], ignore_index=True)
	candidate_audit = pd.concat([tox21_audits[1], sider_audits[1]], ignore_index=True)
	collision_audit = pd.concat([tox21_audits[2], sider_audits[2]], ignore_index=True)
	summary = pd.concat([tox21_audits[3], sider_audits[3]], ignore_index=True)

	for output in [
		tox21_output,
		sider_output,
		match_audit_output,
		candidate_audit_output,
		collision_audit_output,
		summary_output,
	]:
		Path(output).parent.mkdir(parents=True, exist_ok=True)
	tox21_assay.to_csv(tox21_output, sep='\t', index=False)
	sider_assay.to_csv(sider_output, sep='\t', index=False)
	match_audit.to_csv(match_audit_output, sep='\t', index=False)
	candidate_audit.to_csv(candidate_audit_output, sep='\t', index=False)
	collision_audit.to_csv(collision_audit_output, sep='\t', index=False)
	summary.to_csv(summary_output, sep='\t', index=False)


def main_from_snakemake() -> None:
	main(
		coldata_path=str(snakemake.input.colData),
		tox21_input=str(snakemake.input.tox21),
		sider_input=str(snakemake.input.sider),
		tox21_output=str(snakemake.output.tox21),
		sider_output=str(snakemake.output.sider),
		match_audit_output=str(snakemake.output.match_audit),
		candidate_audit_output=str(snakemake.output.candidate_audit),
		collision_audit_output=str(snakemake.output.collision_audit),
		summary_output=str(snakemake.output.match_summary),
		expected_assays=dict(snakemake.params.expected_assays),
	)


if __name__ == '__main__':
	if 'snakemake' in globals():
		main_from_snakemake()
	else:
		parser = argparse.ArgumentParser(
			prog='make_deepchem_experiments',
			description='Generate retained DeepChem experiments and match audits',
		)
		parser.add_argument('-c', '--coldata', required=True)
		parser.add_argument('--tox21-input', required=True)
		parser.add_argument('--sider-input', required=True)
		parser.add_argument('--tox21-output', required=True)
		parser.add_argument('--sider-output', required=True)
		parser.add_argument('--match-audit-output', required=True)
		parser.add_argument('--candidate-audit-output', required=True)
		parser.add_argument('--collision-audit-output', required=True)
		parser.add_argument('--summary-output', required=True)
		args = parser.parse_args()
		main(
			coldata_path=args.coldata,
			tox21_input=args.tox21_input,
			sider_input=args.sider_input,
			tox21_output=args.tox21_output,
			sider_output=args.sider_output,
			match_audit_output=args.match_audit_output,
			candidate_audit_output=args.candidate_audit_output,
			collision_audit_output=args.collision_audit_output,
			summary_output=args.summary_output,
		)
