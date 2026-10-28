from __future__ import annotations

import shutil
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import spatial_regression as sr  # noqa: E402
from spatial_autocorr import Observation, build_knn_weights  # noqa: E402


def simulated_data(size: int = 7, seed: int = 11):
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
    weights, _ = build_knn_weights(observations, k=4)
    w = np.asarray(weights, dtype=float)
    n = len(observations)
    access = np.array([item.access_per_1000 for item in observations])
    elderly = rng.normal(0.2, 0.03, n)
    income = rng.normal(1.0, 0.1, n)
    x = np.column_stack([np.ones(n), access, elderly, income])
    # 空间误差过程：u = 0.45 W u + e
    error = np.linalg.solve(np.eye(n) - 0.45 * w, rng.normal(0, 0.06, n))
    y = 0.2 + 0.01 * access + 0.3 * elderly + 0.05 * income + error
    return observations, y, x, w


class SpatialRegressionTest(unittest.TestCase):
    def test_ols_recovers_simple_signal(self) -> None:
        observations, y, x, w = simulated_data()
        result = sr.fit_ols(y, x, w, ["const", "access_per_1000", "elderly_share", "income_index"])
        by_term = {row["term"]: row for row in result["coefficients"]}
        self.assertAlmostEqual(by_term["access_per_1000"]["estimate"], 0.01, delta=0.005)
        self.assertGreater(result["fits"]["pseudo_r2"], 0.2)

    def test_sem_recovers_positive_spatial_error(self) -> None:
        observations, y, x, w = simulated_data()
        result = sr.fit_sem(y, x, w, ["const", "access_per_1000", "elderly_share", "income_index"])
        self.assertGreater(result["spatial_parameter"], 0.2)
        self.assertLess(result["spatial_parameter"], 0.8)
        self.assertGreater(result["log_likelihood"], -1000)

    def test_write_outputs(self) -> None:
        observations, y, x, w = simulated_data(size=5)
        terms = ["const", "access_per_1000", "elderly_share", "income_index"]
        results = [
            sr.fit_ols(y, x, w, terms),
            sr.fit_slm(y, x, w, terms),
            sr.fit_sem(y, x, w, terms),
        ]
        summary = {
            "n": len(y),
            "fits": [result["fits"] for result in results],
            "lr_tests": {},
        }
        out_dir = ROOT / "tests" / "_tmp_regression"
        if out_dir.exists():
            shutil.rmtree(out_dir, ignore_errors=True)
        out_dir.mkdir()
        try:
            paths = sr.write_outputs(observations, y, results, terms, summary, out_dir, "reg")
            for path in paths.values():
                self.assertTrue(path.exists(), path)
            self.assertIn("model", paths["fits"].read_text(encoding="utf-8-sig"))
        finally:
            shutil.rmtree(out_dir, ignore_errors=True)

    def test_chi2_p_value(self) -> None:
        self.assertAlmostEqual(sr._chi2_1_p(3.841), 0.05, places=2)
        self.assertGreater(sr._chi2_1_p(0), 0.99)


if __name__ == "__main__":
    unittest.main()
