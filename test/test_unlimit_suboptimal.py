#!/usr/bin/env python3
"""
Regression test for the suboptimal-results cap under -u/--unlimit.

create_opt_subopt() always kept only the three suboptimal sites with the
fewest off-targets, even when the run was launched with -u/--unlimit. A run
with unlimited results and no optimal site therefore still produced a
suboptimal file of three rows, when the point of -u is to let the user pick
guides from the full set.

Run from the repository root with:  python3 test/test_unlimit_suboptimal.py
"""

import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.utils import create_opt_subopt  # noqa: E402

CONSTRUCT = "amiRNA"


def make_subopt(n: int) -> list:
    """n suboptimal sites, deliberately unsorted by off-target count."""
    return [
        {
            "off_targets": (i * 7) % n + 1,
            "site": {
                "guide": f"GUIDE{i:02d}",
                "star": f"STAR{i:02d}",
                "oligo1": "", "oligo2": "",
                "tf": f'{{"name": "{CONSTRUCT}{i}"}}',
            },
        }
        for i in range(n)
    ]


class UnlimitSuboptimalTests(unittest.TestCase):

    def test_default_keeps_three_suboptimal_results(self):
        _, subopt_count, _, subopt_results = create_opt_subopt(make_subopt(10), [], CONSTRUCT)
        self.assertEqual(subopt_count, 3)
        self.assertEqual(len(subopt_results), 3)

    def test_unlimit_keeps_every_suboptimal_result(self):
        _, subopt_count, _, subopt_results = create_opt_subopt(
            make_subopt(10), [], CONSTRUCT, unlimit=True)
        self.assertEqual(subopt_count, 10)
        self.assertEqual(len(subopt_results), 10)
        self.assertEqual(sorted(subopt_results), list(range(1, 11)))

    def test_unlimit_results_stay_sorted_by_off_targets(self):
        subopt = make_subopt(10)
        create_opt_subopt(subopt, [], CONSTRUCT, unlimit=True)
        counts = [s["off_targets"] for s in subopt]
        self.assertEqual(counts, sorted(counts))


if __name__ == "__main__":
    unittest.main(verbosity=2)
