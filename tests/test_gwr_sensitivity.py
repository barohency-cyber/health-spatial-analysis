from __future__ import annotations

import shutil
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import gwr_sensitivity as gs  # noqa: E402
from spatial_autocorr import Observation  # noqa: E402


def data(size: int = 7, seed: int = 31):
    rng = np.random.default_rng(seed)
    observations = []
    for x in range(size):
        for y in range(size):
            observations.append(
                Observation(
                    community_id=f"P{x}{y}",
                    name=f"P{x}{y}",
                    lon=120.0 + x * 0.01,
                    lat=30.0 + y * 0.01,
                    population=100,
                    access_per_1000=float(x + y),
                )
            )
    n = len(observations)
    access = np.array([item.access_per_1000 for item in observations])
    x = np.column_stack([
        np.ones(n),
        access,
        rng.normal(0.2, 0.03, n),
        rng.normal(1.0, 0.1, n),
    ])
    y = 0.2 + 0.012 * access + rng.normal(0, 0.03, n)
    return observations, y, x


class GwrSensitivityTest(unittest.TestCase):
    def test_bandwidth_table_and_stability(self) -> None:
        observations, y, x = data()
        table, summary = gs.analyze(observations, y, x, min_k=10, max_k=25, step=5)
        self.assertEqual(len(table), 4)
        self.assertIn(summary["selected_k"], [10, 15, 20, 25])
        self.assertIn("min_pearson_with_selected", summary["stability"])
        self.assertGreaterEqual(summary["stability"]["min_same_sign_share"], 0.0)

    def test_rank_correlation(self) -> None:
        left = np.array([1, 2, 3, 4])
        right = np.array([10, 20, 30, 40])
        self.assertAlmostEqual(gs.rank_correlation(left, right), 1.0)
        self.assertAlmostEqual(gs.rank_correlation(left, left[::-1]), -1.0)

    def test_write_outputs(self) -> None:
        observations, y, x = data(size=5)
        table, summary = gs.analyze(observations, y, x, min_k=8, max_k=18, step=5)
        out_dir = ROOT / "tests" / "_tmp_gwr_sensitivity"
        if out_dir.exists():
            shutil.rmtree(out_dir, ignore_errors=True)
        out_dir.mkdir()
        try:
            paths = gs.write_outputs(table, summary, out_dir, "sens", no_plot=True)
            for path in paths.values():
                self.assertTrue(path.exists(), path)
            self.assertIn("aicc", paths["csv"].read_text(encoding="utf-8-sig"))
        finally:
            shutil.rmtree(out_dir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
