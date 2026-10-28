from __future__ import annotations

import shutil
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import gwr  # noqa: E402
from spatial_autocorr import Observation  # noqa: E402


def grid_data(size: int = 7, seed: int = 23):
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
    elderly = rng.normal(0.2, 0.03, n)
    income = rng.normal(1.0, 0.1, n)
    x = np.column_stack([np.ones(n), access, elderly, income])
    # 让可达性效应在东西方向变化，便于检验局部系数
    local_effect = 0.005 + 0.001 * np.array(
        [obs.lon for obs in observations]
    ) - 0.001 * 120.0
    y = 0.2 + local_effect * access + 0.25 * elderly + 0.05 * income + rng.normal(0, 0.02, n)
    return observations, y, x, ["const", "access_per_1000", "elderly_share", "income_index"]


class GwrTest(unittest.TestCase):
    def test_adaptive_weights_select_neighbors(self) -> None:
        distances = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
        weights = gwr.adaptive_weights(distances, k=3)
        self.assertEqual(weights[0], 1.0)
        self.assertGreater(weights[1], 0)
        self.assertEqual(weights[2], 0.0)
        self.assertEqual(weights[3], 0.0)
        self.assertEqual(weights[4], 0.0)

    def test_gwr_bandwidth_and_coefficients(self) -> None:
        observations, y, x, terms = grid_data()
        rows, summary = gwr.fit_gwr(
            observations, y, x, terms, min_k=10, max_k=35, step=5
        )
        self.assertEqual(len(rows), len(observations))
        self.assertGreaterEqual(summary["bandwidth_k"], 10)
        self.assertLessEqual(summary["bandwidth_k"], 35)
        self.assertEqual(len(summary["candidate_bandwidths"]), 6)
        self.assertGreater(summary["access_coef_summary"]["max"], summary["access_coef_summary"]["min"])

    def test_write_outputs(self) -> None:
        observations, y, x, terms = grid_data(size=5)
        rows, summary = gwr.fit_gwr(observations, y, x, terms, min_k=8, max_k=20, step=4)
        out_dir = ROOT / "tests" / "_tmp_gwr"
        if out_dir.exists():
            shutil.rmtree(out_dir, ignore_errors=True)
        out_dir.mkdir()
        try:
            paths = gwr.write_outputs(
                rows, summary, observations, out_dir, "gwr",
                facilities_path=None, no_map=True,
            )
            for path in paths.values():
                self.assertTrue(path.exists(), path)
            self.assertIn("access_coef", paths["coefficients"].read_text(encoding="utf-8-sig"))
        finally:
            shutil.rmtree(out_dir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
