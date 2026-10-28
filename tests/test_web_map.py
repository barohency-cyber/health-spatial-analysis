from __future__ import annotations

import json
import shutil
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import make_web_map as web  # noqa: E402


class WebMapTest(unittest.TestCase):
    def test_generate_html_and_nojekyll(self) -> None:
        access_rows = web.read_rows(ROOT / "output" / "accessibility.csv")
        facility_rows = web.read_rows(ROOT / "data" / "facilities.csv")
        summary = json.loads(
            (ROOT / "output" / "accessibility_summary.json").read_text(encoding="utf-8")
        )
        out_dir = ROOT / "tests" / "_tmp_web"
        if out_dir.exists():
            shutil.rmtree(out_dir, ignore_errors=True)
        out_dir.mkdir()
        try:
            output = out_dir / "index.html"
            web.build_map(access_rows, facility_rows, summary, output, 5.0)
            page = output.read_text(encoding="utf-8")
            self.assertIn("L.map", page)
            self.assertIn("医疗资源空间可达性", page)
            self.assertIn('"records"', page)
            self.assertTrue((out_dir / ".nojekyll").exists())
        finally:
            shutil.rmtree(out_dir, ignore_errors=True)

    def test_class_edges_starts_at_zero(self) -> None:
        edges = web.class_edges([0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0])
        self.assertEqual(edges[0], 0.0)
        self.assertGreaterEqual(len(edges), 3)


if __name__ == "__main__":
    unittest.main()
