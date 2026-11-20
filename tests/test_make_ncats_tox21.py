"""Exercise the production NCATS parser and Tox21 assay aggregation."""

from __future__ import annotations

import sys
import unittest
from io import StringIO
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow' / 'scripts'))

from make_ncats_tox21 import ASSAYS, collapse_calls, parse_sdf  # noqa: E402


class Tox21Test(unittest.TestCase):
    def test_consensus_missing_and_ambiguous_mapping(self):
        records = [
            {ASSAYS[0]: '0', ASSAYS[1]: '1'},
            {ASSAYS[0]: '1', ASSAYS[2]: '0'},
            {ASSAYS[0]: '1', ASSAYS[2]: '0'},
            {ASSAYS[0]: '0'},
            {ASSAYS[0]: '1'},
            {ASSAYS[0]: '1'},
            {ASSAYS[0]: '1'},
        ]
        compounds, matrix, statuses, discordant = collapse_calls(
            records,
            ['first', 'first', 'first', 'second', 'unmatched', 'ambiguous', ''],
            {'first': {'Z'}, 'second': {'A'}, 'ambiguous': {'A', 'Z'}},
        )
        self.assertEqual(compounds, ['A', 'Z'])
        self.assertEqual(matrix[0], [ASSAYS[0], '0', 'NA'])
        self.assertEqual(matrix[1], [ASSAYS[1], 'NA', '1'])
        self.assertEqual(matrix[2], [ASSAYS[2], 'NA', '0'])
        self.assertEqual(matrix[3], [ASSAYS[3], 'NA', 'NA'])
        self.assertEqual(discordant, 1)
        self.assertEqual(statuses, {
            'unique_hdd_key': 4, 'not_in_hdd': 1, 'ambiguous_hdd_key': 1,
            'no_structure_key': 1,
        })

    def test_mismatched_record_count(self):
        with self.assertRaises(ValueError):
            collapse_calls([{}], [], {})

    def test_parser_does_not_turn_an_absent_call_into_zero(self):
        records = parse_sdf(StringIO('s1\n\nM  END\n> <NR-AR>\n0\n\n$$$$\ns2\nM  END\n$$$$\n'))
        self.assertEqual(records[0]['NR-AR'], '0')
        self.assertNotIn('NR-AR', records[1])

    def test_reject_unexpected_labels(self):
        with self.assertRaisesRegex(ValueError, 'Unexpected NR-AR'):
            parse_sdf(StringIO('s1\n\nM  END\n> <NR-AR>\n2\n\n$$$$\n'))

    def test_reject_truncated_input(self):
        with self.assertRaisesRegex(ValueError, 'unterminated'):
            parse_sdf(StringIO('s1\n\nM  END\n'))


if __name__ == '__main__':
    unittest.main()
