from __future__ import annotations

import math
import shutil
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import two_step_fca as fca  # noqa: E402


class TwoStepFcaTest(unittest.TestCase):
    def setUp(self) -> None:
        # 赤道附近，方便用已知距离构造小例子
        self.d1 = fca.DemandPoint("D1", "A", 0.0, 0.0, 1000)
        self.d2 = fca.DemandPoint("D2", "B", 0.02, 0.0, 500)
        self.f1 = fca.Facility("F1", "one", 0.01, 0.0, 100)
        self.f2 = fca.Facility("F2", "two", 0.05, 0.0, 200)

    def test_hand_calculated_binary(self) -> None:
        result = fca.compute_2sfca(
            [self.d1, self.d2], [self.f1, self.f2], threshold_km=4.0, decay="binary"
        )
        # F2 距 D1 约 5.56 km，超出 4 km；F2 只服务 D2。
        # R1 = 100 / (1000 + 500)，R2 = 200 / 500。
        expected_a = 100 / 1500 * 1000
        expected_b = (100 / 1500 + 200 / 500) * 1000
        self.assertAlmostEqual(result.records[0].access_per_1000, expected_a, places=6)
        self.assertAlmostEqual(result.records[1].access_per_1000, expected_b, places=6)
        self.assertEqual(result.records[0].facilities_in_catchment, 1)
        self.assertEqual(result.records[1].facilities_in_catchment, 2)
        self.assertAlmostEqual(
            result.records[0].nearest_facility_km,
            fca.haversine_km(0.0, 0.0, 0.01, 0.0),
            places=9,
        )

    def test_facility_outside_all_catchments_gets_zero_ratio(self) -> None:
        result = fca.compute_2sfca(
            [self.d1, self.d2], [self.f1, self.f2], threshold_km=2.0, decay="binary"
        )
        self.assertEqual(result.facilities_without_demand, 1)
        self.assertEqual(result.facility_ratios[1], 0.0)

    def test_gaussian_weight_is_lower_at_threshold(self) -> None:
        near = fca.distance_weight(1.0, 5.0, "gaussian")
        far = fca.distance_weight(4.5, 5.0, "gaussian")
        self.assertGreater(near, far)
        self.assertAlmostEqual(fca.distance_weight(5.0, 5.0, "gaussian"), math.exp(-2), places=9)
        self.assertEqual(fca.distance_weight(5.01, 5.0, "gaussian"), 0.0)

    def test_summary_and_outputs(self) -> None:
        result = fca.compute_2sfca(
            [self.d1, self.d2], [self.f1, self.f2], threshold_km=4.0, decay="binary"
        )
        summary = fca.summarize(result, threshold_km=4.0, decay="binary")
        self.assertAlmostEqual(summary["total_population"], 1500)
        self.assertAlmostEqual(summary["coverage_rate"], 1.0)
        tmp = ROOT / "tests" / "_tmp_output"
        if tmp.exists():
            shutil.rmtree(tmp, ignore_errors=True)
        tmp.mkdir()
        try:
            paths = fca.write_outputs(result, summary, tmp, prefix="test")
            for path in paths.values():
                self.assertTrue(path.exists(), path)
            self.assertIn("community_id", paths["csv"].read_text(encoding="utf-8-sig"))
            self.assertIn("FeatureCollection", paths["geojson"].read_text(encoding="utf-8"))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_invalid_threshold(self) -> None:
        with self.assertRaises(ValueError):
            fca.compute_2sfca([self.d1], [self.f1], threshold_km=0)


if __name__ == "__main__":
    unittest.main()
