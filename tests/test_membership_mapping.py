import sys
import unittest
from collections import defaultdict
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parents[1] / 'workflow' / 'scripts'
sys.path.insert(0, str(SCRIPT_DIR))

from process_annotationdb import (  # noqa: E402
	MISSING_MEMBERSHIP_ID,
	append_membership_columns,
	build_membership_map,
)


class MembershipMappingTests(unittest.TestCase):
	def test_membership_map_filters_deduplicates_and_sorts_ids(self) -> None:
		frame = pd.DataFrame(
			{
				'InChIKey': ['AAA', 'AAA', 'AAA', 'BBB', None, 'CCC'],
				'OASIS.ID': ['OASIS2', 'OASIS1', 'OASIS1', 'OASIS3', 'OASIS4', None],
				'Perturbation.Type': [
					'treatment',
					'treatment',
					'control',
					'treatment',
					'treatment',
					'treatment',
				],
			}
		)

		membership = build_membership_map(
			frame,
			inchikey_col='InChIKey',
			id_col='OASIS.ID',
			filter_col='Perturbation.Type',
			filter_value='treatment',
		)

		self.assertEqual(membership['AAA'], 'OASIS1|OASIS2')
		self.assertEqual(membership['BBB'], 'OASIS3')
		self.assertEqual(membership['CCC'], MISSING_MEMBERSHIP_ID)
		self.assertNotIn(None, membership)

	def test_missing_required_membership_column_errors(self) -> None:
		frame = pd.DataFrame({'InChIKey': ['AAA']})

		with self.assertRaisesRegex(ValueError, 'OASIS.ID'):
			build_membership_map(
				frame,
				inchikey_col='InChIKey',
				id_col='OASIS.ID',
			)

	def test_append_membership_columns_handles_match_no_match_and_missing_key(
		self,
	) -> None:
		col_data = defaultdict(list)
		membership = {'AAA': 'SRC1'}

		append_membership_columns(
			col_data, membership, 'AAA', 'In GEOM', 'GEOM Source SMILES'
		)
		append_membership_columns(
			col_data, membership, 'BBB', 'In GEOM', 'GEOM Source SMILES'
		)
		append_membership_columns(
			col_data, membership, None, 'In GEOM', 'GEOM Source SMILES'
		)

		self.assertEqual(col_data['In GEOM'], [True, False, False])
		self.assertEqual(
			col_data['GEOM Source SMILES'],
			['SRC1', MISSING_MEMBERSHIP_ID, MISSING_MEMBERSHIP_ID],
		)


if __name__ == '__main__':
	unittest.main()
