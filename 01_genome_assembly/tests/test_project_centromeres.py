"""Check reference-to-assembly projection through both PAF strands."""
import csv
import runpy
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = "scripts" if (ROOT / "scripts").is_dir() else "submission"
SCRIPTS = [ROOT / name / "python/04_project_centromeres.py"
           for name in (SCRIPT_DIR, SCRIPT_DIR + "_zh")]


def paf_row(qname, qs, qe, strand, tname, ts, te):
    matches = min(qe - qs, te - ts)
    aligned = max(qe - qs, te - ts)
    return f"{qname}\t2000\t{qs}\t{qe}\t{strand}\t{tname}\t2000\t{ts}\t{te}\t{matches}\t{aligned}\t60\n"


class ProjectionTests(unittest.TestCase):
    def test_target_to_query(self):
        cases = [
            ("plus_offset", "chr1", "+", "chr1", 1300, 500, 600, (300, 400, 100)),
            ("minus_offset", "chr1", "-", "chr1", 1300, 500, 600, (800, 900, 100)),
            ("different_names", "assembly1", "+", "reference1", 1300, 500, 600, (300, 400, 100)),
            ("plus_unequal", "chr1", "+", "chr1", 1500, 600, 900, (350, 600, 300)),
            ("minus_unequal", "chr1", "-", "chr1", 1500, 600, 900, (600, 850, 300)),
            ("clipped", "chr1", "+", "chr1", 1300, 200, 450, (100, 250, 150)),
            ("target_only_overlap", "chr1", "+", "chr1", 1300, 1150, 1200, (950, 1000, 50)),
            ("outside_target", "chr1", "+", "chr1", 1300, 100, 200, None),
        ]
        with tempfile.TemporaryDirectory() as directory:
            paf = Path(directory) / "alignment.paf"
            for script in SCRIPTS:
                module = runpy.run_path(str(script))
                for label, qname, strand, tname, te, start, end, expected in cases:
                    with self.subTest(script=script.parent.parent.name, case=label):
                        paf.write_text(paf_row(qname, 100, 1100, strand, tname, 300, te))
                        blocks = module["load_paf"](paf)
                        hit = module["lift_interval"](blocks, tname, start, end)
                        if expected is None:
                            self.assertIsNone(hit)
                        else:
                            self.assertIsNotNone(hit)
                            self.assertEqual((hit["qry_start"], hit["qry_end"], hit["overlap_bp"]), expected)
                            self.assertEqual(hit["qry_chr"], qname)
                            self.assertEqual((hit["ref_chr"], hit["ref_start"], hit["ref_end"]), (tname, start, end))
                            self.assertEqual(hit["strand"], strand)
                paf.write_text(paf_row("decoy", 400, 800, "+", "reference1", 1400, 1800)
                               + paf_row("assembly1", 100, 1100, "+", "reference1", 300, 1300))
                hit = module["lift_interval"](module["load_paf"](paf), "reference1", 500, 600)
                self.assertIsNotNone(hit)
                self.assertEqual((hit["qry_chr"], hit["qry_start"], hit["qry_end"]), ("assembly1", 300, 400))

    def test_cli_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paf = root / "alignment.paf"
            paf.write_text(paf_row("assembly1", 100, 1100, "+", "reference1", 300, 1300))
            positions = root / "centromeres.tsv"
            positions.write_text("clade\tSLL\treference\treference1\t500\t600\t100\n")
            for index, script in enumerate(SCRIPTS):
                out = root / str(index)
                subprocess.run([sys.executable, str(script), "--cen-pos", str(positions),
                                "--species", "SLL", "--accession", "reference", "--paf", str(paf),
                                "--sample", "sample", "--outdir", str(out)], check=True, capture_output=True, text=True)
                self.assertEqual((out / "sample_centromeres_final.bed").read_text(),
                                 "assembly1\t300\t400\tsample_assembly1_cen\n")
                with (out / "sample.centromere_projection.tsv").open() as handle:
                    row = next(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual((row["ref_chr"], row["qry_chr"], row["qry_start"], row["qry_end"], row["status"]),
                                 ("reference1", "assembly1", "300", "400", "OK"))


if __name__ == "__main__":
    unittest.main()
