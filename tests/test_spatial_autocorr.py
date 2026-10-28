from __future__ import annotations

import json
import shutil
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import spatial_autocorr as sa  # noqa: E402


def grid_observations(size: int = 5) -> list[sa.Observation]:
    observations = []
    for x in range(size):
        for y in range(size):
            observations.append(
                sa.Observation(
                    community_id=f"P{x}{y}",
                    name=f"P{x}{y}",
                    lon=120.0 + x * 0.01,
                    lat=30.0 + y * 0.01,
                    population=100,
                    access_per_1000=float(x),
                )
            )
    return observations


class SpatialAutocorrTest(unittest.TestCase):
    def test_knn_weights_are_row_standardized(self) -> None:
        observations = grid_observations()
        weights, neighbors = sa.build_knn_weights(observations, k=4)
        for row, nearest in zip(weights, neighbors):
            self.assertEqual(len(nearest), 4)
            self.assertAlmostEqual(sum(row), 1.0)
            self.assertEqual(sum(1 for value in row if value > 0), 4)

    def test_global_moran_detects_gradient(self) -> None:
        observations = grid_observations()
        values = [item.access_per_1000 for item in observations]
        weights, _ = sa.build_knn_weights(observations, k=4)
        result = sa.global_moran(values, weights, permutations=199, seed=7)
        self.assertGreater(result["global_i"], 0.4)
        self.assertLess(result["p_permutation"], 0.05)
        self.assertGreater(result["z_score"], 2.0)

    def test_cluster_classification(self) -> None:
        self.assertEqual(sa.classify_cluster(1.0, 1.0, True), "High-High")
        self.assertEqual(sa.classify_cluster(-1.0, -1.0, True), "Low-Low")
        self.assertEqual(sa.classify_cluster(1.0, -1.0, True), "High-Low")
        self.assertEqual(sa.classify_cluster(-1.0, 1.0, True), "Low-High")
        self.assertEqual(sa.classify_cluster(1.0, 1.0, False), "Not Significant")

    def test_bh_fdr_is_monotone_and_bounded(self) -> None:
        q_values = sa.bh_fdr([0.01, 0.04, 0.03, 0.2])
        self.assertTrue(all(0 <= value <= 1 for value in q_values))
        self.assertLessEqual(q_values[0], q_values[1])
        self.assertLessEqual(q_values[2], q_values[1])

    def test_analyze_and_write_outputs(self) -> None:
        observations = grid_observations()
        rows, summary = sa.analyze(observations, k=4, permutations=99, seed=3, alpha=0.05)
        self.assertEqual(len(rows), 25)
        self.assertEqual(summary["n"], 25)
        self.assertEqual(sum(summary["cluster_counts"].values()), 25)

        out_dir = ROOT / "tests" / "_tmp_lisa"
        if out_dir.exists():
            shutil.rmtree(out_dir, ignore_errors=True)
        out_dir.mkdir()
        try:
            paths = sa.write_outputs(
                observations, rows, summary, out_dir, "lisa",
                facilities_path=None, no_map=True,
            )
            for path in paths.values():
                self.assertTrue(path.exists(), path)
            payload = json.loads(paths["summary"].read_text(encoding="utf-8"))
            self.assertEqual(payload["n"], 25)
        finally:
            shutil.rmtree(out_dir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
