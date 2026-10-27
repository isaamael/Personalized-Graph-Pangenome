import importlib.util
import sys
import unittest
from pathlib import Path

root = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
scripts = root / ("scripts" if (root / "scripts").exists() else "submission")
spec = importlib.util.spec_from_file_location("prepare", scripts / "01_prepare_genotypes.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class PreparationTests(unittest.TestCase):
    def setUp(self):
        self.markers = [dict(marker_id=f"m{i}", chr="chr01", pos=p)
                        for i, p in enumerate([20, 100, 200, 300, 301, 400, 600])]
        self.calls = [dict(sample_id="sample", **{f"m{i}": "NA" for i in range(7)})]
        self.segments = [dict(sample="sample", chr="chr1", first_marker=100,
                              last_marker=200, state="mat"),
                         dict(sample="sample", chr="chr1", first_marker=400,
                              last_marker=500, state="pat")]

    def test_midpoint_and_terminal_assignment(self):
        self.calls[0]["m1"] = "-1"
        header, rows = module.recover(self.markers, self.calls, self.segments)
        self.assertEqual(header, ["sample_id"] + [f"m{i}" for i in range(7)])
        self.assertEqual(rows, [["sample", -1, -1, -1, -1, 1, 1, 1]])

    def test_observed_call_disagreement_stops(self):
        self.calls[0]["m1"] = "1"
        with self.assertRaisesRegex(ValueError, "disagrees"):
            module.recover(self.markers, self.calls, self.segments)

    def test_missing_chromosome_stops(self):
        self.markers[0]["chr"] = "chr2"
        with self.assertRaisesRegex(ValueError, "Missing segments"):
            module.recover(self.markers, self.calls, self.segments)


if __name__ == "__main__":
    unittest.main()
