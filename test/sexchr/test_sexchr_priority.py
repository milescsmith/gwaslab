"""
Sex-chromosome coverage: PAR filter, to_format xymt, filter_bed vs numeric CHR, MT check_ref,
liftover mock for chr 23.
"""
from __future__ import annotations

import glob
import importlib.util
import os
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab import Sumstats
from gwaslab.info.g_Log import Log
from gwaslab.io.io_to_formats import _to_format
from gwaslab.util.util_in_filter_value import _filter_bed, _filter_region

_HERE = os.path.dirname(os.path.abspath(__file__))
_mod_name = "gwaslab_sexchr_nc_constants"
_spec = importlib.util.spec_from_file_location(
    _mod_name,
    os.path.join(_HERE, "nc_constants.py"),
)
_nc = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
sys.modules[_mod_name] = _nc
_spec.loader.exec_module(_nc)
NC_HUMAN_38_X = _nc.NC_HUMAN_38_X

try:
    import sumstats_liftover  # noqa: F401

    _HAVE_SUMSTATS_LIFTOVER = True
except ImportError:
    _HAVE_SUMSTATS_LIFTOVER = False


def _ref_align_digit(status) -> int:
    try:
        v = int(status)
    except (TypeError, ValueError):
        return -1
    return (v // 10) % 10


class TestSexchrParFilter(unittest.TestCase):
    """filter_region(..., 'PAR') uses packaged PAR BED (CHR 23 in BED; numeric X in sumstats)."""

    def test_par_hg19_drops_inside_par1_keeps_outside(self):
        # hg19 PAR1: 60001–2699520 (1-based inclusive for POS in sumstats vs BED conversion in _filter_bed)
        df = pd.DataFrame({"CHR": [23, 23], "POS": [100_000, 50_000_000]})
        out = _filter_region(df, region="PAR", build="19", verbose=False)
        self.assertEqual(len(out), 1)
        self.assertEqual(int(out["POS"].iloc[0]), 50_000_000)

    def test_par_hg38_drops_inside_par1_keeps_outside(self):
        # hg38 PAR1 in packaged BED: chrom 23, 10001–2781479 (0-based start → 1-based POS after filter)
        df = pd.DataFrame({"CHR": [23, 23], "POS": [500_000, 160_000_000]})
        out = _filter_region(df, region="PAR", build="38", verbose=False)
        self.assertEqual(len(out), 1)
        self.assertEqual(int(out["POS"].iloc[0]), 160_000_000)


class TestSexchrToFormatXymt(unittest.TestCase):
    """tofmt / _to_format: xymt_number toggles 23/24/25 vs X/Y/MT for integer CHR."""

    def setUp(self):
        self._tmpdir = tempfile.mkdtemp(prefix="gwaslab_sexchr_fmt_")
        self.meta = {"gwaslab": {"species": "homo sapiens", "study_name": "t"}}

    def tearDown(self):
        shutil.rmtree(self._tmpdir, ignore_errors=True)

    def test_xymt_number_false_writes_x_y_mt(self):
        df = pd.DataFrame(
            {
                "SNPID": ["a", "b", "c"],
                "CHR": pd.Series([23, 24, 25], dtype="Int64"),
                "POS": [1, 1, 1],
                "EA": ["A", "A", "A"],
                "NEA": ["G", "G", "G"],
                "EAF": [0.1, 0.2, 0.3],
                "BETA": [0.0, 0.0, 0.0],
                "SE": [1.0, 1.0, 1.0],
                "P": [0.5, 0.5, 0.5],
                "STATUS": [3800000, 3800000, 3800000],
            }
        )
        base = os.path.join(self._tmpdir, "n")
        _to_format(
            df.copy(),
            path=base,
            fmt="gwaslab",
            tab_fmt="tsv",
            gzip=False,
            meta=self.meta,
            verbose=False,
            output_log=False,
            xymt_number=False,
        )
        paths = glob.glob(os.path.join(self._tmpdir, "*.tsv"))
        self.assertTrue(paths, "expected a .tsv output file")
        out = pd.read_csv(paths[0], sep="\t")
        chrs = {str(x) for x in out["CHR"].tolist()}
        self.assertEqual(chrs, {"X", "Y", "MT"})

    def test_xymt_number_true_keeps_numeric(self):
        df = pd.DataFrame(
            {
                "SNPID": ["a", "b", "c"],
                "CHR": pd.Series([23, 24, 25], dtype="Int64"),
                "POS": [1, 1, 1],
                "EA": ["A", "A", "A"],
                "NEA": ["G", "G", "G"],
                "EAF": [0.1, 0.2, 0.3],
                "BETA": [0.0, 0.0, 0.0],
                "SE": [1.0, 1.0, 1.0],
                "P": [0.5, 0.5, 0.5],
                "STATUS": [3800000, 3800000, 3800000],
            }
        )
        base = os.path.join(self._tmpdir, "y")
        _to_format(
            df.copy(),
            path=base,
            fmt="gwaslab",
            tab_fmt="tsv",
            gzip=False,
            meta=self.meta,
            verbose=False,
            output_log=False,
            xymt_number=True,
        )
        paths = glob.glob(os.path.join(self._tmpdir, "*.tsv"))
        self.assertTrue(paths)
        out = pd.read_csv(paths[0], sep="\t")
        vals = sorted(int(x) for x in out["CHR"].tolist())
        self.assertEqual(vals, [23, 24, 25])


class TestSexchrFilterBedNotation(unittest.TestCase):
    """BED with chrX or NC_* vs sumstats CHR=23 (numeric)."""

    def setUp(self):
        self._tmpdir = tempfile.mkdtemp(prefix="gwaslab_sexchr_bed_")

    def tearDown(self):
        shutil.rmtree(self._tmpdir, ignore_errors=True)

    def test_filter_bed_chrx_vs_numeric_23(self):
        bed = os.path.join(self._tmpdir, "r.bed")
        with open(bed, "w", encoding="utf-8") as f:
            f.write("chrX\t0\t500\n")
        df = pd.DataFrame({"CHR": pd.Series([23], dtype="Int64"), "POS": [100]})
        out = _filter_bed(df, path=bed, chrom="CHR", pos="POS", keep=True, verbose=False)
        self.assertEqual(len(out), 1)

    def test_filter_bed_nc_x_vs_numeric_23(self):
        bed = os.path.join(self._tmpdir, "n.bed")
        with open(bed, "w", encoding="utf-8") as f:
            f.write(f"{NC_HUMAN_38_X}\t0\t500\n")
        df = pd.DataFrame({"CHR": pd.Series([23], dtype="Int64"), "POS": [100]})
        out = _filter_bed(df, path=bed, chrom="CHR", pos="POS", keep=True, verbose=False)
        self.assertEqual(len(out), 1)


class TestSexchrMtCheckRef(unittest.TestCase):
    """Single MT row: check_ref against a tiny chrM FASTA."""

    def test_check_ref_mt_numeric_chr(self):
        tmp = tempfile.mkdtemp(prefix="gwaslab_sexchr_mt_")
        try:
            fa = os.path.join(tmp, "m.fa")
            with open(fa, "w", encoding="utf-8") as f:
                f.write(">chrM\n")
                f.write("ATCGATCG\n")
            df = pd.DataFrame(
                {
                    "SNPID": ["m1"],
                    "CHR": [25],
                    "POS": [2],
                    "EA": ["C"],
                    "NEA": ["T"],
                    "EAF": [0.05],
                    "BETA": [0.1],
                    "SE": [0.2],
                    "P": [0.3],
                    "N": [100],
                    "STATUS": [3800000],
                }
            )
            gl = Sumstats(df, fmt="gwaslab", build="38", verbose=False)
            gl.fix_chr(verbose=False)
            gl.fix_pos(verbose=False)
            gl.check_ref(fa, verbose=False)
            d6 = gl.data["STATUS"].map(_ref_align_digit)
            self.assertEqual((d6 == 8).sum(), 0, "variant should align to reference")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


@unittest.skipUnless(_HAVE_SUMSTATS_LIFTOVER, "sumstats_liftover package required")
class TestSexchrLiftoverMock(unittest.TestCase):
    """Mocked liftover path for numeric chr 23 (no real chain file)."""

    def test_liftover_numeric_chr23(self):
        from gwaslab.hm.hm_liftover_v2 import _liftover_variant

        data = pd.DataFrame(
            {
                "SNPID": ["x1"],
                "CHR": [23],
                "POS": [1000],
                "EA": ["A"],
                "NEA": ["G"],
                "BETA": [0.1],
                "SE": [0.01],
                "P": [0.01],
                "STATUS": [1900000],
            }
        )
        s = Sumstats(data, fmt="gwaslab", build="19", verbose=False)
        s.data["STATUS"] = 1900000

        def mock_lo(df, **kwargs):
            df = df.copy()
            df["CHR_LIFT"] = 23
            df["POS_LIFT"] = df["POS"] + 50
            df["STRAND_LIFT"] = "+"
            return df

        with patch("gwaslab.hm.hm_liftover_v2.get_chain", return_value="dummy_chain"), patch(
            "gwaslab.hm.hm_liftover_v2.liftover_df", side_effect=mock_lo
        ), patch("gwaslab.hm.hm_liftover_v2._process_build", return_value="38"):
            result = _liftover_variant(
                s.data,
                chrom="CHR",
                pos="POS",
                from_build="19",
                to_build="38",
                status="STATUS",
                chain_path="dummy.chain",
                remove=False,
                filter_by_status=True,
                verbose=False,
                log=Log(),
            )
        self.assertGreater(len(result), 0)
        self.assertTrue((result["POS"] == 1050).any())


if __name__ == "__main__":
    unittest.main()
